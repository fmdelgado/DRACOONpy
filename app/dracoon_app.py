import pandas as pd
import streamlit as st
import base64
from dracoon import dracoon
from utils import plot_dracoon_network, plot_volcano
import json
import pickle
import os
import requests
from PIL import Image
from io import BytesIO
protected_indices = ["example"]


def store_data_as_json(file_name, data):
    with open(file_name, 'w') as file:
        json.dump(data, file)


def get_json(filename):
    with open(filename, 'r') as f:
        data = json.load(f)
    return data


def write_file(filename, content):
    with open(filename, 'wb') as file:
        file.write(content)


def read_file(filename):
    with open(filename, 'rb') as file:
        return file.read()


def get_index_for_files(expression_file, condition_file, structure_file):
    expr_data = pd.read_csv(expression_file, index_col=0)
    cond_data = pd.read_csv(condition_file, index_col=0)
    str_data = pd.read_csv(structure_file, index_col=0)

    dracodict = {'expr_data': expr_data, 'cond_data': cond_data, 'str_data': str_data}
    return dracodict


def store_index_in_db(index, name):
    write_file(f"{name}.pkl", pickle.dumps(index))


def load_index_from_db(index_name):
    index = pickle.loads(read_file(f"{index_name}.pkl"))
    return index


def running_page():


    st.title("Run DRaCOoN")
    pass

    try:
        config = get_json("config.json")
        index_name = config["index"]

    except:
        st.info("No dataset found. Please configure one!")
        st.stop()

    index = load_index_from_db(index_name)
    # Set default parameter values
    default_biom_data = index['expr_data']
    default_cond_data = index['cond_data']
    default_significance = 0.01
    default_associations_df = index['str_data']
    default_iters = 1000
    default_verbose = True
    default_matrixform = True
    default_diffmetric = 'shift'


    # Create sliders or input fields for each parameter
    significance = st.slider("significance", min_value=0.0, max_value=1.0, value=default_significance,step=0.00001,
                                   format="%.5f")

    association_measure = st.selectbox("association_measure", ['entropy', 'pearson', 'spearman'], index=0)
    differential_metric = st.selectbox("differential_metric", ['shift', 'absdiff', 'both'], index=0)

    pvalue_adjustment_method = st.selectbox("pvalue_adjustment_method", ['fdr_bh', 'bonferroni', 'sidak', 'holm-sidak',
                                                                         'holm', 'simes-hochberg', 'hommel',  'fdr_by',
                                                                         'fdr_tsbh', 'fdr_tsbky'], index=0)
    dracoon_program = st.selectbox("dracoon_program", ['DR', 'other_program'], index=0)
    pval_method = st.selectbox("pval_method", ['fitted_background', 'background'], index=0)
    iters = st.slider("iters", min_value=100, max_value=10000, value=default_iters, step=100)



    if st.button("Run"):
        dracoon_net = dracoon(biom_data=default_biom_data,
                                  cond_data=default_cond_data,
                                  significance=significance,
                                  association_measure=association_measure,
                                  pvalue_adjustment_method=pvalue_adjustment_method,
                                  dracoon_program=dracoon_program,
                                  associations_df=default_associations_df,
                                  pval_method=pval_method,
                                  iters=iters,
                                  verbose=default_verbose,
                                  matrixform=default_matrixform)
        dracoon_net.run()
        st.success('Done!')
        volcano_fig = plot_volcano(dracoon_net)
        # Display the Plotly figure in Streamlit
        st.plotly_chart(volcano_fig)
        # Apply metric threshold
        if differential_metric == 'both':
            dracoon_net.res_padj = dracoon_net.res_padj[dracoon_net.res_padj['sig'] == 'absdiff_shift']
        elif differential_metric == 'shift':
            dracoon_net.res_padj = dracoon_net.res_padj[dracoon_net.res_padj['sig'].str.contains('shift')]
        elif differential_metric == 'absdiff':
            dracoon_net.res_padj = dracoon_net.res_padj[dracoon_net.res_padj['sig'].str.contains('absdiff')]
        # Display the DataFrame
        # st.dataframe(dracoon_net.res_padj)
        # Display the network
        network_graph = plot_dracoon_network(dracoon_net.res_padj)

        network_graph.show('result.html')
        # Read the HTML file and render in Streamlit
        with open('result.html', 'r') as f:
            html_string = f.read()
        st.components.v1.html(html_string, width=800, height=600)

        # Add a download button for the DataFrame 'dracoon_net.res_padj'
        csv = dracoon_net.res_padj.to_csv(index=False)
        b64 = base64.b64encode(csv.encode()).decode()
        href = f'<a href="data:file/csv;base64,{b64}" download="dracoon_net_{index_name}_res_padj.csv">Download network</a>'
        st.markdown(href, unsafe_allow_html=True)



def upload_page():

    try:
        indices = get_json("index-list.json")
        if indices == None:
            indices = []
    except:
        indices = []
    st.title("Select input data")

    st.markdown("""
       ### Data Requirements
       To run the DRaCOoN algorithm, you need an expression dataset typically obtained through microarray or RNA-Seq. This dataset should contain multiple samples \( M \) screened over two conditions \( A \) and \( B \).
    
        #### Please upload:
        
        - **Expression dataset**: a matrix of shape \( M \times N \) where \( M \) is the number of genes and \( N \) is the number of samples. The matrix should be in the form of a CSV file with genes as columns and samples as rows. The first column should contain the gene sample. The first row should contain the gene names. The remaining cells should contain the expression values.
        - **Condition dataset**: a matrix with two columns in the form of a CSV file with samples as rows and one column named 'condition' that can take the two condition values subjected to differential networking. The first column should contain the sample names, which must coincide with those at the expression dataset.
        - **Structure dataset**: a table of the gene-gene interactions to be evaluated with the columns 'source' and 'target'. Gene names must coincide with those at the expression dataset.
       
       To proceed, please upload your datasets or select one of the available ones.
        """)






    dip = indices + ["Create New"]
    select_index = st.selectbox("Select dataset", options=dip)
    protected_indices = ["example"]

    if select_index == "Create New":
        with st.form(key="index"):
            st.write("##### Upload your data")
            expr_file = st.file_uploader("Step 1 - Upload your expression dataset", type="csv",
                                         accept_multiple_files=False)
            cond_file = st.file_uploader("Step 2 - Upload your condition dataframe", type="csv",
                                         accept_multiple_files=False)
            structure_file = st.file_uploader("Step 3 - Upload your structure file", type="csv",
                                              accept_multiple_files=False)
            name = st.text_input("Step 4 - Choose a name for DRaCOoN experiment")

            button = st.form_submit_button("Upload")

            if button:
                with st.spinner('uploading data'):
                    index = get_index_for_files(expr_file, cond_file, structure_file)
                    index_name = "exp_" + name
                    store_index_in_db(index, name=index_name)
                    indices = indices + [index_name]
                    store_data_as_json("index-list.json", indices)
                st.success("Finished uploading new data")
                st.experimental_rerun()
    else:
        delete = st.button("Delete")
        config = {}
        config["index"] = select_index

        store_data_as_json("config.json", config)

        if delete:
            # Only delete the files if the index is not protected
            if select_index not in protected_indices:
                # Delete the files
                pkl_file = f"{select_index}.pkl"
                if os.path.exists(pkl_file):
                    os.remove(pkl_file)

                indices.remove(select_index)
                store_data_as_json("index-list.json", indices)

            else:

                st.warning("Dataset is protected and cannot be deleted.")

            st.experimental_rerun()


def about_page():

    st.header("Overview")

    st.markdown("""
        _DRaCOoN_ (Differential Regulation and CO-expression Networks) is a data-driven tool optimized for effectively retrieving differential relationships between genes across two distinct conditions. The tool is designed to handle large datasets and offers various working modes for comprehensive analysis. It produces a differential network involving gene pairs whose association changes between two conditions.

        An overview of the algorithm is shown in the figure below:
        """)

    # Download and display the graphical abstract image
    abstract_url = "https://github.com/fmdelgado/DRACOONpy/blob/master/img/graphical_abstract.jpg?raw=true"
    abstract_response = requests.get(abstract_url)
    abstract_image = Image.open(BytesIO(abstract_response.content))
    st.image(abstract_image, caption="Graphical Abstract")

    # Download and display the GitHub logo
    github_logo_url = "https://cdn-icons-png.flaticon.com/512/25/25231.png"
    github_logo_response = requests.get(github_logo_url)
    github_logo_image = Image.open(BytesIO(github_logo_response.content))

    st.markdown("### Code and Documentation")
    col4, col5 = st.columns([1, 9])
    with col4:
        st.image(github_logo_image, width=30)  # Using higher width for better resolution
    with col5:
        st.markdown("[Visit our GitHub Repository](https://github.com/fmdelgado/DRACOONpy)")


def main():
    # Default content displayed when the app is first loaded
    st.markdown(
        f'<p align="center"> <img src="https://github.com/fmdelgado/DRACOONpy/raw/master/img/dracoon_logo.png" width="300"/> </p>',
        unsafe_allow_html=True,
    )
    st.title("DRaCOoN")
    st.subheader("Differential Regulation and CO-expression Networks")
    page = st.sidebar.selectbox(
        "Choose a page",
        ["Run DRaCOoN", "Upload Data", "About"]
    )

    if page == 'Run DRaCOoN':
        running_page()
    elif page == "Upload Data":
        upload_page()
    elif page == "About":
        about_page()

if __name__ == "__main__":
    main()
