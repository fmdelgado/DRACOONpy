from matplotlib import pyplot as plt
from pyvis.network import Network
import matplotlib.colors as mcolors
import pandas as pd
import networkx as nx
import seaborn as sns
import math
import plotly.express as px
from statsmodels.stats.multitest import multipletests


def plot_conditional_distribution(gene1, gene2, exprdata, condata, ra, rb, r, absdiff, pabsdiff, shift, pshift, ax,
                                  conditions_dict):
    conditions = sorted(list(set(conditions_dict.keys()).intersection(set(condata['condition'].astype(str)))))

    data_a = exprdata.loc[condata[condata['condition'] == conditions[0]].index]
    data_b = exprdata.loc[condata[condata['condition'] == conditions[1]].index]

    subdataset_control = data_a[[gene1, gene2]].copy()
    subdataset_control.loc[:, 'condition'] = conditions_dict[conditions[0]]

    subdataset_case = data_b[[gene1, gene2]].copy()
    subdataset_case.loc[:, 'condition'] = conditions_dict[conditions[1]]

    result = pd.concat([subdataset_control, subdataset_case])
    scatterplot = sns.scatterplot(x=gene1, y=gene2, hue="condition", data=result, ax=ax)

    # Remove the legend
    legend_colors = [handle.get_facecolor()[0] for handle in scatterplot.legend_.legendHandles]
    scatterplot.legend_.remove()

    # Add the info inside the plot area in the top left corner
    info_labels = [
        f"r_{conditions_dict[conditions[0]]}: {ra:.3f}\nr_{conditions_dict[conditions[1]]}: {rb:.3f}\nr: {r:.3f}",
        f"absdiff: {absdiff:.3f}, {pabsdiff:.3f}",
        f"shift: {shift:.3f}, {pshift:.3f}"
    ]

    info_text = "\n".join(info_labels)
    ax.text(0.05, 0.95, info_text, fontsize=9, transform=ax.transAxes, verticalalignment='top')

    # Create a custom legend with condition colors
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', label=condition,
                                  markerfacecolor=color, markersize=8)
                       for condition, color in zip(conditions, legend_colors)]

    ax.legend(handles=legend_elements, loc='lower right', title='Condition')


def plot_dracoon_network(subnet):
    # Sort the dataframe
    try:
        subnet = subnet.sort_values(by=['padj_absdiff', 'padj_shift'], ascending=False)
    except:
        try:
            subnet = subnet.sort_values(by=['padj_shift'], ascending=True)
        except:
            subnet = subnet.sort_values(by=['padj_absdiff'], ascending=False)

    # Create a gradient scale from red to yellow (ylOrRd)
    cmap = plt.get_cmap('YlOrRd')
    norm = mcolors.Normalize(vmin=0, vmax=len(subnet) - 1)

    # Assign colors to edges based on their row position
    subnet['color'] = [mcolors.to_hex(cmap(norm(i))) for i in range(len(subnet))]

    # Remove duplicate edges
    subnet = subnet[~subnet[['source', 'target']].duplicated(keep='first')]

    # Initialize the network
    g = Network(notebook=True, directed=False)

    # Add nodes
    nodes_involved = pd.DataFrame({'node': pd.unique(subnet[['source', 'target']].values.ravel('K'))})

    # Calculate betweenness centrality using networkx
    G_nx = nx.from_pandas_edgelist(subnet, source='source', target='target')
    betweenness_centrality = nx.betweenness_centrality(G_nx)
    # print(betweenness_centrality)
    nodes_involved['betweenness'] = nodes_involved['node'].map(betweenness_centrality)

    # Normalize betweenness centrality values for setting node sizes
    min_size, max_size = 10, 30
    norm_betweenness = mcolors.Normalize(vmin=nodes_involved['betweenness'].min(),
                                         vmax=nodes_involved['betweenness'].max())
    nodes_involved['size'] = [min_size + (max_size - min_size) * norm_betweenness(val) for val in
                              nodes_involved['betweenness']]

    g.add_nodes(nodes_involved.node.to_list(), label=nodes_involved.node.astype(str),
                value=nodes_involved['size'].tolist())

    # Add edges with the assigned colors and hover titles
    for src, dest, color, padj_absdiff, padj_shift in zip(subnet['source'], subnet['target'], subnet['color'],
                                                          subnet['padj_absdiff'], subnet['padj_shift']):
        title = f"padj_absdiff: {padj_absdiff:.4f} | padj_shift: {padj_shift:.4f}"

        # Set the edge style based on 'assoc_type' value
        g.add_edge(src, dest, color=color, width=5, style='solid', title=title)

    # Show the network
    return g


def plot_differetial_relationships(dracoonobj, conditions_dict=None):
    # Iterate over the top 10 rows and plot the conditional distributions on the grid
    df = dracoonobj.res_padj.sort_values(by=['padj_absdiff', 'padj_shift'], ascending=True)
    df = df[(df['padj_shift'] < 0.1)]
    if not conditions_dict:
        conditions_dict = {str(dracoonobj.conditions[0]): str(dracoonobj.conditions[0]),
                           str(dracoonobj.conditions[1]): str(dracoonobj.conditions[1])}

    if not df.empty:
        # Calculate the grid size based on the length of df
        grid_size = int(math.ceil(math.sqrt(len(df))))
        fig, axes = plt.subplots(grid_size, grid_size, figsize=(grid_size * 3, grid_size * 3))

        for i, (index, row) in enumerate(df.iterrows()):
            ax = axes[i // grid_size, i % grid_size]
            plot_conditional_distribution(gene1=row['source'], gene2=row['target'], exprdata=dracoonobj.biom_data,
                                          condata=dracoonobj.cond_data,
                                          ra=row[f'r_{str(dracoonobj.conditions[0])}'],
                                          rb=row[f'r_{str(dracoonobj.conditions[1])}'],
                                          r=row['r'], absdiff=row['absdiff'],
                                          pabsdiff=row['padj_absdiff'],
                                          shift=row['shift'], pshift=row['padj_shift'], ax=ax,
                                          conditions_dict=conditions_dict)
        # Remove empty subplots if any
        for j in range(i + 1, grid_size * grid_size):
            fig.delaxes(axes[j // grid_size, j % grid_size])

        # Adjust the layout and show the figure
        plt.tight_layout()
        fig = plt.gcf()  # get the current figure
        fig.set_dpi(300)  # Set the DPI to 300
        plt.show()


def plot_volcano(dracobj):
    thres = dracobj.significance
    finalnet = dracobj.res_p.copy()
    finalnet['padj_absdiff'] = multipletests(finalnet['p_absdiff'], method=dracobj.pvalue_adjustment_method)[1]
    finalnet['padj_shift'] = multipletests(finalnet['p_shift'], method=dracobj.pvalue_adjustment_method)[1]

    def assign_color(row):
        if row['padj_shift'] < thres and row['padj_absdiff'] < thres:
            return 'both'
        elif row['padj_shift'] < thres:
            return 'padj_shift'
        elif row['padj_absdiff'] < thres:
            return 'padj_absdiff'
        else:
            return 'none'

    finalnet['significance_group'] = finalnet.apply(assign_color, axis=1)
    fig = px.scatter(finalnet, x='absdiff', y='shift',
                     color='significance_group', hover_data=['source', 'target', 'padj_shift', 'padj_absdiff'])
    fig.update_layout(title='Overview of differential associations')  # Add title here

    return fig
