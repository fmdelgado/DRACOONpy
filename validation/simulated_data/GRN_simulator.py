import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
import seaborn as sns
import scipy.stats as ss
import networkx as nx
from pyvis.network import Network
import os
import math


class GRN_simulator:
    def __init__(self, n_genes, n_tfs, n_samples,
                 n_shutdowns, n_inversions, n_shifts, n_uncorrelateds, proportion_group_case=0.5,
                 simdata_noise=0.1, perturbation_noise=0.1, mean_expression=2, shift_alpha=2,
                 ingroup_perturbed_proportion=1, zscore_normalize=False):

        # Structure info
        self.structure = None
        self.genes = None
        self.TFs = None

        # Perturbation inform.
        self.shutdown_bioms = None
        self.inverted_bioms = None
        self.shifted_bioms = None
        self.loc_bioms = None
        self.perturbed_bioms = None
        self.shift_alpha = None

        # Dataset info
        self.simulated_data = None
        self.condition_data = None
        self.control_df = None
        self.case_df = None
        self.adjacency = None
        self.refnet = None
        self.zscore_normalize = zscore_normalize

        # Details biomolecules
        self.n_genes = n_genes
        self.n_tfs = n_tfs  # genes considered TFs
        # Details samples
        self.proportion_group_case = proportion_group_case
        self.n_samples = n_samples
        # Details simulation
        self.simdata_noise = simdata_noise  # from 0 to 1, corresponds to SD
        self.mean_expression = mean_expression

        # Details perturbation
        self.ingroup_perturbed_proportion = ingroup_perturbed_proportion
        self.perturbation_noise = perturbation_noise  # from 0 to 1, corresponds to SD
        self.n_shutdowns = n_shutdowns  # genes to shut down
        self.n_inversions = n_inversions  # genes to invert
        self.n_shifts = n_shifts
        self.n_uncorrelateds = n_uncorrelateds
        self.shift_alpha = shift_alpha

        if self.n_samples < 10:
            exit('ERROR: Number of samples should be over 10.')

    def generate_DRN_graphsim(self):
        # print('Estimating control dataset...')
        self.structure = self.recreate_tree_structure()
        # We add information about activation (1) or repression (-1). In principle, they are all activating.
        grn_structure_path = '~/Documents/Research/DraCooN/checkups/data_simulator/graphsim_structure.csv'
        self.structure.to_csv(grn_structure_path, index=False)
        correlation_level = 1
        outpath = '/Users/fernando/Documents/Research/DraCooN/checkups/graphsim/graphsim_simdata.csv'
        os.system(
            f"~/local/bin/Rscript /Users/fernando/Documents/Research/DraCooN/algorithm/graphism_grnsimulator.R {grn_structure_path} {self.n_samples} {outpath} {correlation_level} {self.mean_expression}")
        simdata = pd.read_csv(outpath, index_col=0)
        # print('Estimating perturbed dataset...')
        self.simulated_data, self.condition_data, self.shutdown_bioms, self.inverted_bioms, self.perturbed_bioms = self.split_and_perturb_simdata(
            simdata=simdata)
        self.simulated_data.index = list(map(lambda each: each[7:], self.simulated_data.index.to_list()))
        self.condition_data.index = list(map(lambda each: each[7:], self.condition_data.index.to_list()))
        # print('Estimating reference differential network...')
        self.refnet = self.get_reference_DN()

    def generate_DRN(self):
        # print('Estimating control dataset...')
        self.recreate_tree_structure()
        simdata = self.simulate_expression_from_net()

        # print('Estimating perturbed dataset...')
        self.simulated_data = self.split_and_perturb_simdata(simdata=simdata)
        # print('Estimating reference differential network...')
        self.refnet = self.get_reference_DN()
        # print('Done!')

    def simulate_expression_from_adjacency(self, adjamat=None, n_samples=None, simdata_noise=None):
        if adjamat is not None:
            self.adjacency = adjamat
        if n_samples is not None:
            self.n_samples = n_samples
        if simdata_noise is not None:
            self.simdata_noise = simdata_noise

        # noise = 0.5 # .1 strong correlation, 3 weak correlation
        # random initialize data matrix
        n_rows = self.adjacency.shape[0]  # number of genes
        mu, sigma = self.mean_expression, 1  # mean and standard deviation
        dataset = np.random.normal(mu, sigma, self.n_samples * n_rows).reshape(n_rows, self.n_samples)

        already_modified = np.repeat(0, n_rows)  # N=n1+n2 genes
        already_modified[0] = 1  # leave the first gene alone, base case

        for i in range(self.adjacency.shape[1]):
            for j in range(i + 1, self.adjacency.shape[1]):
                # print(f'Considering row: {i}, column: {j} of A')
                if (self.adjacency[i, j] == 1) & (already_modified[j] == 0):
                    dataset[j, :] = dataset[i, :] + np.random.normal(0, self.simdata_noise, dataset.shape[1])
                    already_modified[j] = 1
                elif (self.adjacency[i, j] == -1) & (already_modified[j] == 0):
                    dataset[j, :] = -dataset[i, :] + np.random.normal(0, self.simdata_noise, dataset.shape[1])
                    already_modified[j] = 1
                elif (already_modified[j] == 1) & (already_modified[i] == 0):
                    # if j is  modified, we accordingly modify i. If i has been modified, we do nothing.
                    dataset[i, :] = dataset[j, :] + np.random.normal(0, self.simdata_noise, dataset.shape[1])
                    already_modified[i] = 1
        return dataset

    def recreate_tree_structure(self):
        self.genes = ['g' + str(s) for s in list(range(self.n_genes))]
        # print('Estimating obtaining TF list...')
        self.TFs = random.sample(self.genes, self.n_tfs)
        targets = list(set(self.genes) - set(self.TFs))
        targets_lists = list(self.split(targets, self.n_tfs))

        tree = nx.random_tree(n=len(self.TFs), create_using=nx.DiGraph)
        # print(nx.forest_str(tree))
        tftfdf = nx.to_pandas_edgelist(tree)
        tf_mapping = dict(zip(range(len(self.TFs)), self.TFs))
        tftfdf = tftfdf.applymap(tf_mapping.get)
        tftfdf['type'] = 'TF-TF'

        df = tftfdf
        # Create a dictionary to store the level of each node
        levels = {}

        # Initialize the level of the root node to 0
        root_node = list(set(df['source']) - set(df['target']))[0]
        levels[root_node] = 0

        # Define a function to set the level of each node
        def set_levels(node):
            children = df[df['source'] == node]['target']
            for child in children:
                if child not in levels:
                    levels[child] = levels[node] + 1
                    set_levels(child)

        # Set the levels of all nodes
        set_levels(root_node)

        # Add the levels to the DataFrame
        df['level'] = df['source'].map(levels)

        # Sort the DataFrame based on the levels
        tftfdf = df.sort_values('level')

        source = []
        target = []

        for tf, target_list in zip(self.TFs, targets_lists):
            for thetarget in target_list:
                source.append(tf)
                target.append(thetarget)
                # print(tf, thetarget)
        tftargetdf = pd.DataFrame({'source': source, 'target': target})
        tftargetdf['type'] = 'TF-gene'

        self.structure = pd.concat(
            [tftfdf[['source', 'target', 'type']], tftargetdf])  # we join tf-tf rels with tf-gene rels
        # self.structure['weight'] = np.random.choice(a=[1, -1], size=(self.structure.shape[0]), p=[p, 1 - p])
        self.structure['weight'] = 1
        self.structure['color'] = self.structure['weight'].map(dict({1: 'blue', -1: 'red'}))
        self.structure.reset_index(inplace=True, drop=True)

        # nx.draw(nx.from_pandas_edgelist(self.structure, edge_attr=["weight", "color"],  create_using=nx.DiGraph), with_labels=True)
        # plt.show()
        return self.structure

    def simulate_expression_from_net(self, structure_net=None, n_samples=None, simdata_noise=None):
        if n_samples is not None:
            self.n_samples = n_samples
        if simdata_noise is not None:
            self.simdata_noise = simdata_noise
        if structure_net is not None:
            self.structure = structure_net

        # Initialize biomolecuels to simulate with random parameters
        n_rows = len(self.genes)
        mu, sigma = self.mean_expression, 1  # mean and standard deviation
        dataset = pd.DataFrame(np.random.normal(mu, sigma, self.n_samples * n_rows).reshape(n_rows, self.n_samples))
        dataset.index = self.genes
        # sns.heatmap(dataset.T.corr(),  vmin=-1, vmax=1, annot=False, cmap='RdBu')
        # plt.show()
        # dataset = self.simulated_data.T

        for source, target, rel in zip(self.structure.source, self.structure.target, self.structure.weight):
            # print(source, target, rel)
            dataset.loc[target] = dataset.loc[source].to_numpy() + np.random.normal(0, self.simdata_noise,
                                                                                    dataset.shape[1])

            r, _ = ss.pearsonr(dataset.loc[source], dataset.loc[target])
        return dataset

    @staticmethod
    def split(a, n):
        k, m = divmod(len(a), n)
        return (a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))

    def split_and_perturb_simdata(self, simdata, n_shutdowns=None, n_inversions=None, n_shifts=None,
                                  perturbation_noise=None, proportion_group_case=None):
        if n_shutdowns is not None:
            self.n_shutdowns = n_shutdowns
        if n_inversions is not None:
            self.n_inversions = n_inversions
        if n_shifts is not None:
            self.n_shifts = n_shifts
        if perturbation_noise is not None:
            self.perturbation_noise = perturbation_noise
        if proportion_group_case is not None:
            self.proportion_group_case = proportion_group_case

        groupsize = int(self.proportion_group_case * simdata.shape[1])
        case = simdata.iloc[:, :groupsize]
        control = simdata.iloc[:, groupsize:]
        ingroup_perturbed = int(self.ingroup_perturbed_proportion * groupsize)
        # print(ingroup_perturbed)

        # We avoid repetition so inverted genes cannot be shut down
        self.perturbed_bioms = np.random.choice(simdata.index.to_list(),
                                                self.n_shutdowns + self.n_inversions + self.n_shifts + self.n_uncorrelateds,
                                                replace=False)
        self.shutdown_bioms, self.inverted_bioms, self.shifted_bioms, self.loc_bioms, _ = np.split(self.perturbed_bioms,
                                                                                                   np.cumsum([
                                                                                                                 self.n_shutdowns,
                                                                                                                 self.n_inversions,
                                                                                                                 self.n_shifts,
                                                                                                                 self.n_uncorrelateds]))

        case = case.T
        # We now shut down genes
        if self.shutdown_bioms.size != 0:
            for biom in self.shutdown_bioms:
                # biom = self.shifted_bioms[0]
                zeroed = self.estimate_correlating_vector_given_rho(x1=case[biom].to_numpy(), rho=0) + np.random.normal(
                    0, self.perturbation_noise, groupsize)
                # zeroed = np.random.normal(0, self.perturbation_noise, groupsize).reshape(groupsize, )
                mask = np.array([1] * ingroup_perturbed + [0] * (groupsize - ingroup_perturbed)).astype(bool)
                np.random.shuffle(mask)
                case[biom][mask] = zeroed[mask]

        # We now invert genes
        if self.inverted_bioms.size != 0:
            for biom in self.inverted_bioms:
                inverted = self.estimate_correlating_vector_given_rho(x1=case[biom].to_numpy(),
                                                                      rho=-0.99) + np.random.normal(
                    self.mean_expression, self.perturbation_noise, groupsize)
                # inverted = -case.loc[biom] + np.random.normal(0, self.perturbation_noise, groupsize)
                mask = np.array([1] * ingroup_perturbed + [0] * (groupsize - ingroup_perturbed)).astype(bool)
                np.random.shuffle(mask)
                case[biom][mask] = inverted[mask]

        # We now shift some genes
        if self.shifted_bioms.size != 0:
            for biom in self.shifted_bioms:
                shifted = case[biom] + np.random.normal(self.shift_alpha * (np.mean(case[biom])),
                                                        self.perturbation_noise, groupsize)
                mask = np.array([1] * ingroup_perturbed + [0] * (groupsize - ingroup_perturbed)).astype(bool)
                np.random.shuffle(mask)
                case[biom][mask] = shifted[mask]

        # We now shuffle some genes
        if self.loc_bioms.size != 0:
            for biom in self.loc_bioms:
                uncorrelated = np.random.permutation(case[biom].to_numpy()) + np.random.normal(0,
                                                                                               self.perturbation_noise,
                                                                                               groupsize)
                mask = np.array([1] * ingroup_perturbed + [0] * (groupsize - ingroup_perturbed)).astype(bool)
                np.random.shuffle(mask)
                case[biom][mask] = uncorrelated[mask]

        case = case.T
        self.simulated_data = pd.concat([control, case], axis=1).T
        self.simulated_data.index = 'sample_' + self.simulated_data.index.astype(str)
        convector = [[1] * control.shape[1], [2] * case.shape[1]]
        convector = [x for xs in convector for x in xs]

        self.condition_data = pd.DataFrame({'sample': self.simulated_data.index, 'condition': convector})
        self.condition_data.set_index('sample', inplace=True)

        if self.zscore_normalize:
            # self.case_df = ss.zscore( self.case_df)
            # self.control_df = ss.zscore(self.control_df)
            self.simulated_data = ss.zscore(self.simulated_data.T).T

        self.case_df = self.simulated_data.loc[
            self.condition_data.loc[self.condition_data.condition == 2].index.to_list()]
        self.control_df = self.simulated_data.loc[
            self.condition_data.loc[self.condition_data.condition == 1].index.to_list()]
        '''test = self.simulated_data
        test[test.columns].plot(kind='kde', alpha=0.5)
        plt.show()'''
        return self.simulated_data

    def get_reference_DN(self):
        df1 = self.structure[self.structure['target'].isin(self.perturbed_bioms)]
        df2 = self.structure[self.structure['source'].isin(self.perturbed_bioms)]
        df = pd.concat([df1, df2], ignore_index=True)
        self.refnet = df.drop_duplicates(subset=['source', 'target'], keep='last')
        return self.refnet

    def check_type(self, value):
        if value in self.TFs:
            return 'TF'
        elif value in self.genes:
            return 'gene'

    @staticmethod
    def ground_truth_from_adjamat(adjamat, biom_names):
        gene_a = []
        gene_b = []
        rels = []
        for i in range(adjamat.shape[1]):
            for j in range(i + 1, adjamat.shape[1]):
                if adjamat[i][j] != 0:
                    # print(i, j, adjamat[i, j])
                    gene_a.append(i)
                    gene_b.append(j)
                    rels.append(np.where(adjamat[i, j] > 0, "[+]", "[-]"))

        pairs = pd.DataFrame({'source': gene_a, 'target': gene_b, 'correlation': rels})
        mygenedict = dict(zip(range(adjamat.shape[1]), biom_names))
        pairs.source = pairs.source.map(mygenedict)
        pairs.target = pairs.target.map(mygenedict)
        return pairs

    # plottings
    def plot_pyvis(self, relationships=None):
        if relationships is None:
            relationships = self.structure

        # input_genes = example.expression_data
        # Tfs = example.TFs
        # relationships = example.net_case
        g = Network(notebook=True, directed=True)

        nodes_involved = pd.DataFrame({'node': pd.unique(relationships[['source', 'target']].values.ravel('K'))})
        # Adding all the genes that
        # took part
        nodes_involved['biom_type'] = nodes_involved['node'].apply(self.check_type)  # Identify the Tfs
        conditions = [nodes_involved['biom_type'] == 'TF',
                      nodes_involved['biom_type'] == 'gene',
                      nodes_involved['biom_type'] == 'cpg']
        values = ['diamond', 'circle', 'triangle']

        # create a list of the values we want to assign for each condition
        nodes_involved['shape'] = np.select(conditions, values)
        nodes_involved['perturbed'] = np.where(nodes_involved['node'].isin(self.perturbed_bioms), 'red', 'skyblue')
        # print(nodes_involved)

        # print(nodes_involved)
        g.add_nodes(nodes_involved.node.to_list(), title=nodes_involved.biom_type.to_list(),
                    shape=nodes_involved['shape'].to_list(), color=nodes_involved.perturbed.to_list())  # Add the nodes

        for src, dest in zip(relationships['source'], relationships['target']):
            # print(src, dest)
            g.add_edge(src, dest)

        g.show('mytftarget.html')
        return g

    def plot_clustermap(self, cluster=True):
        my_palette = dict(zip(self.condition_data.condition.unique(), ["blue", "orange"]))
        row_colors = self.condition_data.condition.map(my_palette)
        sns.clustermap(self.simulated_data, cmap="coolwarm", row_colors=row_colors, figsize=(15, 10), center=0.00,
                       row_cluster=False,
                       col_cluster=cluster)
        plt.show()

    def plot_conditional_distributions(self, gene1, gene2, dataset_a=None, dataset_b=None):
        if dataset_a is None:
            dataset_a = self.control_df
        if dataset_b is None:
            dataset_b = self.case_df

        subdataset_control = dataset_a[[gene1, gene2]]
        subdataset_control['condition'] = 'control'
        subdataset_case = dataset_b[[gene1, gene2]]
        subdataset_case['condition'] = 'case'
        result = pd.concat([subdataset_control, subdataset_case])
        sns.jointplot(x=gene1, y=gene2, hue="condition", data=result)
        plt.show()

    def estimate_correlating_vector_given_rho(self, x1, rho):
        n = len(x1)  # length of vector
        theta = math.acos(rho)  # corresponding angle

        x2 = np.random.normal(self.mean_expression, self.simdata_noise, n)
        X = np.vstack((x1, x2)).T

        Xctr = ss.zscore(X)  # centered columns (mean 0)
        Id = np.diag(np.ones(n))  # identity matrix
        Q = np.linalg.qr(Xctr)[0][:, 0]  # QR-decomposition, just matrix Q
        P = Q.reshape(-1, 1) @ Q.reshape(1, -1)  # projection onto space defined by x1
        x2o = (Id - P) @ Xctr[:, 1]  # x2ctr made orthogonal to x1ctr
        Xc2 = np.vstack((Xctr[:, 0], x2o)).T  # bind to matrix
        Y = Xc2 @ np.diag(1 / np.sum(Xc2 ** 2, axis=0) ** 0.5)  # scale columns to length 1
        x = Y[:, 1] + (1 / math.tan(theta)) * Y[:, 0]  # final new vector
        return x


# Working example

'''

---------------------------------------------------------------------------------------------------

n_iters = 1000
num_genes = 50
num_shuts = int(0.2 * num_genes) #20% of genes will be shut down
num_TFs = int(0.3 * num_genes) #30% of genes will be TFs

affected = int(np.round(num_shuts/4))

simulation = GRN_simulator(n_genes=num_genes, n_tfs=num_TFs, n_samples=300,
                           n_shutdowns=affected, n_shifts =affected, 
                           n_inversions=affected, n_uncorrelateds=affected, 
                           proportion_group_case=0.5,
                           simdata_noise=2,
                           perturbation_noise=2)
simulation.generate_DRN()
#simulation.simulated_data = ss.zscore(simulation.simulated_data)

simulation.plot_clustermap(cluster=True)
print(simulation.perturbed_bioms)



biom2perturbation = {**{elem:'shifted' for elem in simulation.shifted_bioms},
                     **{elem:'shutdown' for elem in simulation.shutdown_bioms},
                    **{elem:'inverted' for elem in simulation.inverted_bioms},
                    **{elem:'loc' for elem in simulation.loc_bioms}}

unperturbed_genes = set(simulation.simulated_data.columns) - set(biom2perturbation.keys())
biom2perturbation = {**biom2perturbation,
                     **{elem:'unperturbed' for elem in unperturbed_genes}}

df = simulation.refnet.copy()

conditions_dict = {'1' : 'control',  '2':'case'}
size = 50

if (not df.empty):
    if len(df) > size:
        df = df.head(size)
        print(f"plotting first {size}")

    # Calculate the grid size based on the length of df
    grid_size = int(math.ceil(math.sqrt(len(df))))
    fig, axes = plt.subplots(grid_size, grid_size, figsize=(grid_size * 3, grid_size * 3))

    for i, (index, row) in enumerate(df.iterrows()):
        ax = axes[i // grid_size, i % grid_size]

        gene1 = row['source']
        gene2 = row['target']
        expr_data  = simulation.simulated_data.copy()
        cond_data = simulation.condition_data.copy()


        conditions = sorted(list(set(conditions_dict.keys()).intersection(set(cond_data['condition'].astype(str)))))
        #print(conditions)


        data_A = expr_data.loc[cond_data[cond_data['condition'] == int(conditions[0])].index]
        data_B = expr_data.loc[cond_data[cond_data['condition'] == int(conditions[1])].index]

        subdataset_control = data_A[[gene1, gene2]].copy()
        subdataset_control.loc[:, 'condition'] = conditions_dict[conditions[0]]

        subdataset_case = data_B[[gene1, gene2]].copy()
        subdataset_case.loc[:, 'condition'] = conditions_dict[conditions[1]]

        result = pd.concat([subdataset_control, subdataset_case])
        result = result.merge(cond_data, left_index=True, right_index=True)
        # print(result)
        palette = {1: 'tab:blue', 2: 'tab:orange'}


        scatterplot = sns.scatterplot(x=gene1, y=gene2, hue="condition_y", data=result, ax=ax, palette=palette,)
        ax.set(xlabel=f'{gene1} {biom2perturbation[gene1]}', ylabel=f'{gene2} {biom2perturbation[gene2]}')
        ax.legend(loc='lower right', title='Condition')
        ax.get_legend().remove()


    # Remove empty subplots if any
    for j in range(i+1, grid_size * grid_size):
        fig.delaxes(axes[j // grid_size, j % grid_size])


    # After the for loop, add a general legend to the figure
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center',  bbox_to_anchor=(0.9, 0.2))

    # Adjust the layout and show the figure
    plt.tight_layout()
    fig = plt.gcf()  # get the current figure
    fig.set_dpi(300)  # Set the DPI to 300
    plt.show()

self = simulation
'''
