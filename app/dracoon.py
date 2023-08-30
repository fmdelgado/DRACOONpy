import sys
import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from statsmodels.stats.multitest import multipletests
from fitter import Fitter, get_common_distributions, get_distributions
import scipy.stats as ss
import scipy.special as sc
import re
from itertools import product
# time measures
from tqdm import tqdm
import time
import datetime
# parallelization
import functools
# optimizer
TRY_USE_NUMBA = True


def jit_if_available(func):
    # default "do nothing" decorator with the numba-like interface
    def decorated(*args, **kwargs):
        return func(*args, **kwargs)
    return decorated


if TRY_USE_NUMBA:
    try:
        import numba as nb  # as jit_if_available

        jit_if_available = nb.njit()
        print("Numba is available. Optimization on.")
    except:
        print("Numba is not available. Install numba for a bit faster calculations")


# @nb.njit
def pearsoncorr(gene1, gene2):
    r = np.corrcoef(gene1, gene2)[0, 1]
    return 0 if np.isnan(r) else r


# @nb.njit
def numpyspearmancorr(gene1, gene2):
    ranked_gene1 = np.empty_like(gene1)
    ranked_gene1[np.argsort(gene1)] = np.arange(len(gene1))
    ranked_gene1 += 1

    ranked_gene2 = np.empty_like(gene2)
    ranked_gene2[np.argsort(gene2)] = np.arange(len(gene2))
    ranked_gene2 += 1

    r = pearsoncorr(ranked_gene1, ranked_gene2)

    return r


# @nb.njit
def test_significance(observed_value, simulated_results, two_sided=True):
    greater = np.sum(simulated_results >= observed_value)
    smaller = np.sum(simulated_results <= observed_value)
    if two_sided:
        p_value = min(greater, smaller) / (len(simulated_results) - 1)
        return 2 * p_value
    else:
        return greater / (len(simulated_results) - 1)


# @nb.njit
def entropy_between_2_genes(gene_a_vals, gene_b_vals, epsilon=np.finfo(float).eps):
    rho = pearsoncorr(gene_a_vals, gene_b_vals)
    if rho == 1:  # remedies error arising from  divide by zero encountered in log2
        rho = rho - epsilon
    elif rho == -1:
        rho = rho + epsilon
    entropy = -((1 + abs(rho)) / 2 * np.log((1 + abs(rho)) / 2) + (1 - abs(rho)) / 2 * np.log((1 - abs(rho)) / 2))
    return entropy


# @nb.njit
def simulate_values_for_differential_params(gene_i_values_a, gene_i_values_b, gene_j_values_a, gene_j_values_b,
                                            association_func):
    shuff_gene_i_values_a = np.random.permutation(gene_i_values_a)
    shuff_gene_i_values_b = np.random.permutation(gene_i_values_b)
    genei = np.concatenate((shuff_gene_i_values_a, shuff_gene_i_values_b), axis=0)
    np.random.shuffle(genei)

    shuff_gene_j_values_a = np.random.permutation(gene_j_values_a)
    shuff_gene_j_values_b = np.random.permutation(gene_j_values_b)
    genej = np.concatenate((shuff_gene_j_values_a, shuff_gene_j_values_b), axis=0)
    np.random.shuffle(genej)

    sim_r_a = association_func(shuff_gene_i_values_a, shuff_gene_j_values_a)
    sim_r_b = association_func(shuff_gene_i_values_b, shuff_gene_j_values_b)
    sim_r = association_func(genei, genej)

    sim_abs_dif = abs(sim_r_a - sim_r_b)
    sim_shift = (sim_r_a + sim_r_b) / 2 - sim_r
    return sim_abs_dif, sim_shift


# @nb.njit
def njitn_random_genepairs(numpydata_a, numpydata_b, am, niters):
    background_distribution_absdiff = []
    background_distribution_shift = []

    ngenes = numpydata_a.shape[1]
    # pbar = tqdm(total = niters)
    # while len(background_distribution_absdiff) < niters:
    for _ in range(niters):
        try:
            gene_i, gene_j = np.random.choice(ngenes, 2, replace=False)

            simul_abs_dif, simul_shift = simulate_values_for_differential_params(
                gene_i_values_a=numpydata_a[:, gene_i],
                gene_i_values_b=numpydata_b[:, gene_i],
                gene_j_values_a=numpydata_a[:, gene_j],
                gene_j_values_b=numpydata_b[:, gene_j],
                association_func=am)

            background_distribution_absdiff.append(simul_abs_dif)
            background_distribution_shift.append(simul_shift)
            # pbar.update(1)
        except:
            pass
    return np.asarray(background_distribution_absdiff), np.asarray(background_distribution_shift)


# @nb.njit
def run_pearson(kombi, data_a, data_b, biom_data):
    ra = pearsoncorr(data_a[:, kombi[0]], data_a[:, kombi[1]])
    rb = pearsoncorr(data_b[:, kombi[0]], data_b[:, kombi[1]])
    rall = pearsoncorr(biom_data[:, kombi[0]], biom_data[:, kombi[1]])
    return ra, rb, rall


# @nb.njit
def run_spearman(kombi, data_a, data_b, biom_data):
    ra = numpyspearmancorr(data_a[:, kombi[0]], data_a[:, kombi[1]])
    rb = numpyspearmancorr(data_b[:, kombi[0]], data_b[:, kombi[1]])
    rall = numpyspearmancorr(biom_data[:, kombi[0]], biom_data[:, kombi[1]])
    return ra, rb, rall


# @nb.njit
def run_entropy(kombi, data_a, data_b, biom_data):
    ra = entropy_between_2_genes(data_a[:, kombi[0]], data_a[:, kombi[1]])
    rb = entropy_between_2_genes(data_b[:, kombi[0]], data_b[:, kombi[1]])
    rall = entropy_between_2_genes(biom_data[:, kombi[0]], biom_data[:, kombi[1]])
    return ra, rb, rall


# @nb.njit
def print_correlation_line(kombi, data_a, data_b, biom_data, bg_abs, bg_shift, alpha_shift=2, fun=run_pearson):
    ra, rb, rall = fun(kombi, data_a, data_b, biom_data)
    # Compute differential metrics
    abs_dif = abs(ra - rb)
    shift = (ra + rb) / alpha_shift - rall
    # Compute simulated values for the differential parameters n_iter times, 10000 by default
    abs_dif_pval = test_significance(simulated_results=bg_abs, observed_value=abs_dif, two_sided=False)
    shift_pval = test_significance(simulated_results=bg_shift, observed_value=shift, two_sided=True)
    return float(kombi[0]), float(kombi[1]), ra, rb, rall, abs_dif, abs_dif_pval, shift, shift_pval


# @nb.njit(parallel=True)
def prangeprint_correlations(k1, k2, data_a, data_b, biom_data, bg_abs, bg_shift, am, alpha_shift=2):
    arr = np.empty((len(k1), 9), dtype=nb.f4)
    for i in nb.prange(len(k1)):
        # print(k1[i])
        arr[i, :] = print_correlation_line(kombi=(k1[i], k2[i]), data_a=data_a, data_b=data_b, biom_data=biom_data,
                                           bg_abs=bg_abs, bg_shift=bg_shift, alpha_shift=alpha_shift, fun=am)
    return arr


# @nb.njit
def print_correlation_line_nopval(kombi, data_A, data_B, biom_data, alpha_shift=2, fun=run_pearson):
    ra, rb, rall = fun(kombi, data_A, data_B, biom_data)
    # Compute differential metrics
    abs_dif = abs(ra - rb)
    shift = (ra + rb) / alpha_shift - rall
    # Compute simulated values for the differential parameters n_iter times, 10000 by default
    return float(kombi[0]), float(kombi[1]), ra, rb, rall, abs_dif, shift


# @nb.njit(parallel=True)
def prangeprint_correlations_nopval(k1, k2, data_A, data_B, biom_data, am, alpha_shift=2):
    arr = np.empty((len(k1), 7), dtype=nb.f4)
    for i in nb.prange(len(k1)):
        # print(k1[i])
        arr[i, :] = print_correlation_line_nopval(kombi=(k1[i], k2[i]), data_A=data_A, data_B=data_B,
                                                  biom_data=biom_data, alpha_shift=alpha_shift, fun=am)
    return arr


# @nb.njit
def random_dist_genepair(gene_i_values_a, gene_i_values_b, gene_j_values_a, gene_j_values_b, am, niters):
    background_distribution_absdiff = []
    background_distribution_shift = []

    for _ in nb.prange(niters):
        try:
            simul_abs_dif, simul_shift = simulate_values_for_differential_params(
                gene_i_values_a=gene_i_values_a,
                gene_i_values_b=gene_i_values_b,
                gene_j_values_a=gene_j_values_a,
                gene_j_values_b=gene_j_values_b,
                association_func=am)

            background_distribution_absdiff.append(simul_abs_dif)
            background_distribution_shift.append(simul_shift)
            # pbar.update(1)
        except:
            pass
    return np.asarray(background_distribution_absdiff), np.asarray(background_distribution_shift)


class DraCooN_GRN:
    def __init__(self, biom_data, cond_data, dracoon_program, significance=0.05,
                 association_measure='pearson', pval_method='background', pvalue_adjustment_method='fdr_bh',
                 associations_df=None, association_pvalue_filter=None, iters=None, timeout_fitter=None,
                 zscore_normalize=False, log2_transform=False,
                 verbose=True, matrixform=True):
        # Input Data
        self.cond_data = cond_data
        self.biom_data = biom_data
        self.associations_df = associations_df
        self.source_targets = None
        self.assessed_interactions = None
        self.not_assessed_interactions = None
        self.conditions = None
        self.position_to_biom = None
        self.biom_to_position = None
        self.targets = None
        self.sources = None
        self.data_A = None
        self.data_B = None
        self.timeout_fitter = None

        # Algorithm parameters
        self.verbose = verbose
        self.dracoon_program = dracoon_program  # DR: differential regulation or DC: differential co-expression
        self.iters = iters
        self.significance = significance
        self.association_measure = association_measure
        self.alpha_shift = 2
        self.association_pvalue_filter = association_pvalue_filter
        self.pvalue_adjustment_method = pvalue_adjustment_method
        self.log2_transform = log2_transform
        self.zscore_normalize = zscore_normalize
        self.matrixform = matrixform

        # Running mode
        self.pval_method = pval_method

        # Results
        self.background_distribution_absdiff = None
        self.background_distribution_shift = None
        self.res_unthresholded = None
        self.res_p = None
        self.res_padj = None

        if self.iters is None:
            ngenes = self.biom_data.shape[1]
            pairwises = int(ngenes * (ngenes - 1) / 2)
            self.iters = int(10 ** np.log10(pairwises).round()) * 10
            if self.verbose:
                print(f"Number of permutations not specified, running using {self.iters}")

        if self.timeout_fitter is None:
            self.timeout_fitter = 60

        if self.dracoon_program not in ['DC', 'DR', 'DCR']:
            exit("ERROR: invalid Dracoon program!")

        asoc_list = ['spearman', 'pearson', 'entropy']
        if self.association_measure not in asoc_list:
            exit("ERROR: invalid correlation measure!")

        pval_methods = ['permutation', 'background', 'fitted_background']
        if self.pval_method not in pval_methods:
            exit('ERROR: invalid p-value obtention method!')

        pvaladj_methods = ['bonferroni', 'sidak', 'holm-sidak', 'holm', 'simes-hochberg', 'hommel', 'fdr_bh', 'fdr_by',
                           'fdr_tsbh', 'fdr_tsbky']
        if self.pvalue_adjustment_method not in pvaladj_methods:
            exit('ERROR: invalid p-value adjustment method!')

    def run(self):
        self.preprocessing()
        # Check if background distribution is needed
        if self.pval_method != 'permutation':
            self.estimate_background_model()

        t = time.time()
        if self.pval_method == 'permutation':
            self.calculate_correlations_perm()

        elif self.matrixform:
            self.calculate_correlations_smallinput()
        else:
            if self.pval_method == 'background':
                self.calculate_correlations_bg()
            elif self.pval_method == 'fitted_background':
                self.calculate_correlations_fittedbg()

        # Show process
        if self.verbose:
            print(f'Time for correlations estimation: {time.time() - t}')

        self.threshold_results()

    def preprocessing(self):
        # We assume columns correspond to genes and rows correspond to samples.
        # We anyway retrieve all the combinatios that could be picked for the background model generation
        self.position_to_biom = dict(zip(range(self.biom_data.shape[1]), self.biom_data.columns))
        self.biom_to_position = dict(zip(self.biom_data.columns, range(self.biom_data.shape[1])))

        if self.dracoon_program == 'DC':
            combs = np.vectorize(self.position_to_biom.get)(np.triu_indices(self.biom_data.shape[1], 1))
            self.source_targets = list(zip(combs[0], combs[1]))

        elif self.dracoon_program == 'DR':
            # Check which of the interactions meant to study actually contain bioms in the quantitative matrix df
            self.assessed_interactions = self.associations_df[
                (self.associations_df.source.isin(self.biom_data.columns)) & (
                    self.associations_df.target.isin(self.biom_data.columns))]
            # Check which of the interactions meant to study don't contain bioms in the quantitative matrix df
            self.not_assessed_interactions = self.associations_df[
                (~self.associations_df.source.isin(self.biom_data.columns)) | (
                    ~self.associations_df.target.isin(self.biom_data.columns))]
            # Check which biomolecules are exclusively sources and which are exclusively targets
            self.sources = set(self.assessed_interactions.source)
            self.targets = set(self.assessed_interactions.target).difference(self.sources)
            # boi = self.targets.union(self.sources)  # get unique biomolecules
            # Check the dataset contains only bioms in associations -> NO! It affects background model
            # self.biom_data = self.biom_data[self.biom_data.columns.intersection(boi)]
            self.source_targets = list(zip(self.assessed_interactions.source, self.assessed_interactions.target))

        elif self.dracoon_program == 'DCR':
            # Check which of the interactions meant to study actually contain bioms in the quantitative matrix df
            self.assessed_interactions = self.associations_df[
                (self.associations_df.source.isin(self.biom_data.columns))]
            # Check which of the interactions meant to study don't contain bioms in the quantitative matrix df
            self.not_assessed_interactions = self.associations_df[
                ~self.associations_df.source.isin(self.biom_data.columns)]
            self.sources = set(self.assessed_interactions.source)
            self.targets = set(self.biom_data.columns)
            self.source_targets = list(product(self.sources, self.targets))
            # Remove self loops
            self.source_targets[:] = [x for x in self.source_targets if x[0] != x[1]]

            # Divide biomolecules data to the condition data frame and perform some transformations

        if 'condition' in self.cond_data.columns:
            self.conditions = list(self.cond_data['condition'].unique())
            self.conditions.sort()
        else:
            exit('ERROR: Condition column not found')

        if len(self.conditions) > 2:
            exit('ERROR: more than 2 conditions are not supported!')

        if self.log2_transform:
            if self.verbose:
                print('log2(n+1) transforming input data...')
            self.biom_data = np.log2(self.biom_data + 1)

        if self.zscore_normalize:
            if self.verbose:
                print('z-score normalizing input data...')
            self.biom_data = ss.zscore(self.biom_data, axis=0)
            '''
            print('conditional z-score normalizing input data...')
            r1 = ss.zscore(self.biom_data.loc[self.cond_data[self.cond_data['condition'] == self.conditions[0]].index])
            r2 = ss.zscore(self.biom_data.loc[self.cond_data[self.cond_data['condition'] == self.conditions[1]].index])
            self.biom_data = pd.concat([r1, r2], axis=0)
            '''
        # control
        common_samples = set(self.biom_data.index).intersection(set(self.cond_data.index))

        if len(common_samples) != len(self.biom_data.index):
            print(
                f"WARNING: No condition found for the following {len(set(self.biom_data.index).difference(common_samples))} samples:\n"
                f"{set(self.biom_data.index).difference(common_samples)}\n"
                f"Estimation will be performed for correctly labeled samples present in dataset.")
        if len(common_samples) != len(self.cond_data.index):
            print(
                f"WARNING: The following {len(set(self.cond_data.index).difference(common_samples))} samples in the condition table were not found in the dataset:\n"
                f"{set(self.cond_data.index).difference(common_samples)}\n"
                f"Estimation will be performed for correctly labeled samples present in dataset.")

        self.cond_data = self.cond_data[self.cond_data.index.isin(common_samples)]
        self.data_A = self.biom_data.loc[
            self.cond_data.index[self.cond_data['condition'] == self.conditions[0]].tolist()]
        # case
        self.data_B = self.biom_data.loc[
            self.cond_data.index[self.cond_data['condition'] == self.conditions[1]].tolist()]

        if self.verbose:  # Show process
            print(f'Condition {self.conditions[0]}: ' + str(self.data_A.shape[0]) + ' samples')
            print(f'Condition {self.conditions[1]}: ' + str(self.data_B.shape[0]) + ' samples')
            print(f"{self.dracoon_program} Relationships to estimate: {len(self.source_targets)}")

    def estimate_background_model(self):
        metric_to_func = {'pearson': pearsoncorr, 'spearman': numpyspearmancorr, 'entropy': entropy_between_2_genes}
        bgabsdiff, bgshift = njitn_random_genepairs(numpydata_a=np.array(self.data_A),
                                                    numpydata_b=np.array(self.data_B),
                                                    niters=self.iters,
                                                    am=metric_to_func[self.association_measure])
        self.background_distribution_absdiff = np.sort(bgabsdiff)
        self.background_distribution_shift = np.sort(bgshift)

    def calculate_correlations_smallinput(self, background_absdiff=None, background_shift=None):
        if background_absdiff is not None:
            self.background_distribution_absdiff = background_absdiff
        if background_shift is not None:
            self.background_distribution_shift = background_shift
        if self.verbose:
            # self =result_DR
            print(f'Estimating associations using {self.association_measure}')

        ra, pa = self.corrcoef_pval(nparray=np.array(self.data_A))
        rb, pb = self.corrcoef_pval(nparray=np.array(self.data_B))
        r, p = self.corrcoef_pval(nparray=np.array(self.biom_data))
        absdiff = abs(ra - rb)
        shift = (ra + rb) / self.alpha_shift - r

        if self.pval_method == 'background':
            absdiffpval = list(map(functools.partial(test_significance,
                                                     simulated_results=self.background_distribution_absdiff,
                                                     two_sided=False), absdiff))
            shiftpval = list(map(functools.partial(test_significance,
                                                   simulated_results=self.background_distribution_shift,
                                                   two_sided=True), shift))
        elif self.pval_method == 'fitted_background':
            best_dists_absdiff = ['expon', 'logistic', 'rayleigh', 'norm', 'gumbel_r', 'laplace', 'kstwobign',
                                  'moyal', 'beta']
            best_dists_shift = ['logistic', 'norm', 'laplace', 'gumbel_l', 'gumbel_r', 'uniform', 'expon', 'rayleigh',
                                'hypsecant']
            fitdist_absdiff, params_absdiff = self.estimate_bestdist_params(
                distribution=self.background_distribution_absdiff, possible_distributions=best_dists_absdiff,
                timeout=self.timeout_fitter, dist_name='absdiff')
            fitdist_shift, params_shift = self.estimate_bestdist_params(distribution=self.background_distribution_shift,
                                                                        possible_distributions=best_dists_shift,
                                                                        timeout=self.timeout_fitter, dist_name='shift')
            absdiffpval = list(map(functools.partial(self.estimate_pval_fitdist, dist=fitdist_absdiff,
                                                     param=params_absdiff, two_sided=False), absdiff))
            shiftpval = list(map(functools.partial(self.estimate_pval_fitdist, dist=fitdist_shift,
                                                   param=params_shift, two_sided=True), shift))

        res = np.column_stack((ra, pa))
        res = np.column_stack((res, rb, pb))
        res = np.column_stack((res, r, p))
        res = np.column_stack((res, absdiff, absdiffpval))
        res = np.column_stack((res, shift, shiftpval))

        if self.dracoon_program == 'DR':
            # res_thr = pd.DataFrame(np.column_stack((self.assessed_interactions[['source', 'target']], res)))
            res_thr = pd.concat(
                [self.assessed_interactions[['source', 'target']].reset_index(drop=True), pd.DataFrame(res)], axis=1,
                ignore_index=True)
        else:
            # res_thr = pd.DataFrame(np.column_stack((self.source_targets, res)))
            res_thr = pd.concat([pd.DataFrame(self.source_targets), pd.DataFrame(res)], axis=1)

        res_thr.columns = ["source", "target", "r_A", "p_A", "r_B", "p_B", "r", "p_r", "absdiff", "p_absdiff", 'shift',
                           'p_shift']
        res_thr.dropna(inplace=True)
        res_thr.reset_index(drop=True)
        # print("Thresholding to get significant correlations")
        if self.verbose:
            print(f'\nCorrecting for multiple testing using {self.pvalue_adjustment_method}')
        res_thr['padj_absdiff'] = multipletests(res_thr['p_absdiff'], method=self.pvalue_adjustment_method)[1]
        res_thr['padj_shift'] = multipletests(res_thr['p_shift'], method=self.pvalue_adjustment_method)[1]

        res_thr['padj_A'] = multipletests(res_thr['p_A'], method=self.pvalue_adjustment_method)[1]
        res_thr['padj_B'] = multipletests(res_thr['p_B'], method=self.pvalue_adjustment_method)[1]
        res_thr['padj_r'] = multipletests(res_thr['p_r'], method=self.pvalue_adjustment_method)[1]

        column_order = ['source', 'target', 'r_A', 'padj_A', 'r_B', 'padj_B', 'r', 'padj_r', 'absdiff', 'p_absdiff',
                        'padj_absdiff', 'shift', 'p_shift', 'padj_shift']
        self.res_unthresholded = res_thr[column_order]  # Unfiltered
        return self.res_unthresholded

    def calculate_correlations_bg(self, combis=None):
        if combis is None:
            if self.dracoon_program == 'DC':
                combis = self.source_targets
            elif self.dracoon_program in ['DR', 'DCR']:
                combis = self.source_targets

        metric_to_funcs = {'pearson': run_pearson, 'spearman': run_spearman, 'entropy': run_entropy}
        corrfunc = metric_to_funcs[self.association_measure]
        n = 100000  # chunksize
        now = str(datetime.datetime.now()).replace(' ', '')
        tempfilename = f'tmp_dracoon_perm_{now}.tsv'

        with open(tempfilename, 'wb') as f:
            pvals = np.empty((0, 5), np.float64)

            for i in tqdm(range(0, len(combis), n)):
                kombis = [(self.biom_to_position[combi[0]], self.biom_to_position[combi[1]]) for combi in
                          combis[i:i + n]]
                k1 = [combi[0] for combi in kombis]
                k2 = [combi[1] for combi in kombis]
                resk1 = prangeprint_correlations(k1, k2, data_a=np.array(self.data_A), data_b=np.array(self.data_B),
                                                 biom_data=np.array(self.biom_data),
                                                 bg_abs=self.background_distribution_absdiff,
                                                 bg_shift=self.background_distribution_shift,
                                                 alpha_shift=self.alpha_shift,
                                                 am=corrfunc)

                # resk1[:,[0,1]] = resk1[:,[0,1]].astype(int)
                pvals = np.vstack((pvals, resk1[:, [2, 3, 4, 6, 8]]))
                np.savetxt(f, resk1, delimiter='\t')

        pvals[:, 0] = self.estimate_corrpval(r=pvals[:, 0], corrpoints=self.data_A.shape[0])
        pvals[:, 1] = self.estimate_corrpval(r=pvals[:, 1], corrpoints=self.data_B.shape[0])
        pvals[:, 2] = self.estimate_corrpval(r=pvals[:, 2], corrpoints=self.biom_data.shape[0])

        if self.verbose:
            print(f'\nCorrecting for multiple testing using: {self.pvalue_adjustment_method}')
        padjusteds = np.empty(pvals.shape, np.float64)

        for col in range(pvals.shape[1]):
            padjusteds[:, col] = multipletests(pvals[:, col], method=self.pvalue_adjustment_method)[1]

        sign_lines = np.where(np.sum(padjusteds[:, 3:5] < self.significance, axis=1) > 0)[0]
        # np.sum(padjusteds < signif, axis=0)
        array = []
        with open(tempfilename, 'r') as f:
            lines = f.readlines()
            for myline in sign_lines:
                array.append([float(x) for x in lines[myline].strip('\n').split('\t')])

        results_p = pd.concat([pd.DataFrame(array), pd.DataFrame(padjusteds[sign_lines, :])], axis=1)
        results_p.columns = ["source", "target", "r_A", "r_B", "r", "absdiff", "p_absdiff", 'shift',
                             'p_shift', 'padj_A', 'padj_B', 'padj_r', 'padj_absdiff', 'padj_shift']
        results_p['source'] = results_p['source'].astype(int).map(self.position_to_biom)
        results_p['target'] = results_p['target'].astype(int).map(self.position_to_biom)

        # print("Thresholding to get significant correlations")
        results_p.dropna(inplace=True)
        results_p.reset_index(drop=True)
        # In regular DC we remove self differential associations, unless these are specified in the source - target file
        if self.dracoon_program == 'DC':
            results_p = results_p[results_p['source'] != results_p['target']]
        self.res_unthresholded = results_p[["source", "target", "r_A", "padj_A", "r_B", "padj_B", "r", "padj_r",
                                            "absdiff", "p_absdiff", 'padj_absdiff', 'shift', 'p_shift', 'padj_shift']]
        os.system(f'rm -rf {tempfilename}')
        return self.res_unthresholded

    def calculate_correlations_fittedbg(self, combis=None):
        if combis is None:
            if self.dracoon_program == 'DC':
                combis = self.source_targets
            elif self.dracoon_program in ['DR', 'DCR']:
                combis = self.source_targets

        # Estimating the fitted distribution
        best_dists_absdiff = ['expon', 'logistic', 'rayleigh', 'norm', 'gumbel_r', 'pareto', 'laplace', 'kstwobign',
                              'moyal', 'halfnorm']
        best_dists_shift = ['logistic', 'norm', 'laplace', 'gumbel_l', 'gumbel_r', 'uniform', 'expon', 'rayleigh',
                            'hypsecant']
        fitdist_absdiff, params_absdiff = self.estimate_bestdist_params(
            distribution=self.background_distribution_absdiff,
            possible_distributions=best_dists_absdiff,
            timeout=self.timeout_fitter, dist_name='absdiff')
        fitdist_shift, params_shift = self.estimate_bestdist_params(distribution=self.background_distribution_shift,
                                                                    possible_distributions=best_dists_shift,
                                                                    timeout=self.timeout_fitter, dist_name='shift')
        # Saving results per iteration
        metric_to_funcs = {'pearson': run_pearson, 'spearman': run_spearman, 'entropy': run_entropy}
        corrfunc = metric_to_funcs[self.association_measure]
        n = 100000  # chunksize
        now = str(datetime.datetime.now()).replace(' ', '')
        tempfilename = f'tmp_dracoon_perm_{now}.tsv'
        with open(tempfilename, 'wb') as f:
            pvals = np.empty((0, 5), np.float64)

            for i in tqdm(range(0, len(combis), n)):
                kombis = [(self.biom_to_position[combi[0]], self.biom_to_position[combi[1]]) for combi in
                          combis[i:i + n]]
                k1 = [combi[0] for combi in kombis]
                k2 = [combi[1] for combi in kombis]
                resk1 = prangeprint_correlations_nopval(k1, k2, data_A=np.array(self.data_A),
                                                        data_B=np.array(self.data_B),
                                                        biom_data=np.array(self.biom_data),
                                                        alpha_shift=self.alpha_shift,
                                                        am=corrfunc)

                # resk1[:,[0,1]] = resk1[:,[0,1]].astype(int)
                # ra, rb, rall, absdiff, shift
                pvals = np.vstack((pvals, resk1[:, [2, 3, 4, 5, 6]]))
                np.savetxt(f, resk1, delimiter='\t')

        pvals[:, 0] = self.estimate_corrpval(r=pvals[:, 0], corrpoints=self.data_A.shape[0])
        pvals[:, 1] = self.estimate_corrpval(r=pvals[:, 1], corrpoints=self.data_B.shape[0])
        pvals[:, 2] = self.estimate_corrpval(r=pvals[:, 2], corrpoints=self.biom_data.shape[0])
        pvals[:, 3] = list(map(functools.partial(self.estimate_pval_fitdist, dist=fitdist_absdiff,
                                                 param=params_absdiff, two_sided=False), pvals[:, 3]))
        pvals[:, 4] = list(map(functools.partial(self.estimate_pval_fitdist, dist=fitdist_shift,
                                                 param=params_shift, two_sided=True), pvals[:, 4]))

        if self.verbose:
            print(f'\nCorrecting for multiple testing using: {self.pvalue_adjustment_method}')
        padjusteds = np.empty(pvals.shape, np.float64)

        for col in range(pvals.shape[1]):
            padjusteds[:, col] = multipletests(pvals[:, col], method=self.pvalue_adjustment_method)[1]

        sign_lines = np.where(np.sum(padjusteds[:, 3:5] < self.significance, axis=1) > 0)[0]
        # np.sum(padjusteds < signif, axis=0)
        array = []
        with open(tempfilename, 'r') as f:
            lines = f.readlines()
            for myline in sign_lines:
                array.append([float(x) for x in lines[myline].strip('\n').split('\t')])

        results_p = pd.concat(
            [pd.DataFrame(array), pd.DataFrame(pvals[sign_lines, :]), pd.DataFrame(padjusteds[sign_lines, :])], axis=1)
        results_p.columns = ["source", "target", "r_A", "r_B", "r", "absdiff", 'shift', 'p_A', 'p_B', 'p_r',
                             'p_absdiff', 'p_shift', 'padj_A', 'padj_B', 'padj_r', 'padj_absdiff', 'padj_shift']
        results_p['source'] = results_p['source'].astype(int).map(self.position_to_biom)
        results_p['target'] = results_p['target'].astype(int).map(self.position_to_biom)

        # print("Thresholding to get significant correlations")
        results_p.dropna(inplace=True)
        results_p.reset_index(drop=True)
        # In regular DC we remove self differential associations, unless these are specified in the source - target file
        if self.dracoon_program == 'DC':
            results_p = results_p[results_p['source'] != results_p['target']]

        colorder = ["source", "target", "r_A", "p_A", "padj_A", "r_B", "p_B", "padj_B", "r", "p_r", "padj_r", "absdiff",
                    "p_absdiff", "padj_absdiff", "shift", "p_shift", "padj_shift"]

        self.res_unthresholded = results_p[colorder]
        os.system(f'rm -rf {tempfilename}')
        return self.res_unthresholded

    def calculate_correlations_perm(self, combis=None):
        if combis is None:
            if self.dracoon_program == 'DC':
                combis = self.source_targets
            elif self.dracoon_program in ['DR', 'DCR']:
                combis = self.source_targets

        metric_to_funcs = {'pearson': run_pearson, 'spearman': run_spearman, 'entropy': run_entropy}
        corrfunc = metric_to_funcs[self.association_measure]
        metric_to_funcs = {'pearson': pearsoncorr, 'spearman': numpyspearmancorr, 'entropy': entropy_between_2_genes}
        minifunc = metric_to_funcs[self.association_measure]

        n = 100000  # chunksize
        now = str(datetime.datetime.now()).replace(' ', '')
        tempfilename = f'tmp_dracoon_perm_{now}.tsv'
        with open(tempfilename, 'wb') as f:
            pvals = np.empty((0, 5), np.float64)

            for i in range(0, len(combis), n):
                if self.verbose:
                    print(f'Running on chunk {i + 1}')
                kombis = [(self.biom_to_position[combi[0]], self.biom_to_position[combi[1]]) for combi in
                          combis[i:i + n]]
                arr = np.empty((len(kombis), 9))

                for j in tqdm(range(len(kombis))):
                    arr[j, :] = self.print_correlation_inlinepermutation(kombi=kombis[j],
                                                                         data_A=np.array(self.data_A),
                                                                         data_B=np.array(self.data_B),
                                                                         biom_data=np.array(self.biom_data),
                                                                         fun=corrfunc,
                                                                         minifun=minifunc)

                pvals = np.vstack((pvals, arr[:, [2, 3, 4, 6, 8]]))
                np.savetxt(f, arr, delimiter='\t')

        pvals[:, 0] = self.estimate_corrpval(r=pvals[:, 0], corrpoints=self.data_A.shape[0])
        pvals[:, 1] = self.estimate_corrpval(r=pvals[:, 1], corrpoints=self.data_B.shape[0])
        pvals[:, 2] = self.estimate_corrpval(r=pvals[:, 2], corrpoints=self.biom_data.shape[0])

        if self.verbose:
            print(f'\nCorrecting for multiple testing using {self.pvalue_adjustment_method}')
        padjusteds = np.empty(pvals.shape, np.float64)

        for col in range(pvals.shape[1]):
            padjusteds[:, col] = multipletests(pvals[:, col], method=self.pvalue_adjustment_method)[1]

        sign_lines = np.where(np.sum(padjusteds < self.significance, axis=1) > 0)[0]
        # np.sum(padjusteds < signif, axis=0)
        array = []

        with open(tempfilename, 'r') as f:
            lines = f.readlines()
            for myline in sign_lines:
                array.append([float(x) for x in lines[myline].strip('\n').split('\t')])

        results_p = pd.concat([pd.DataFrame(array), pd.DataFrame(padjusteds[sign_lines, :])], axis=1)
        results_p.columns = ["source", "target", "r_A", "r_B", "r", "absdiff", "p_absdiff", 'shift',
                             'p_shift', 'padj_A', 'padj_B', 'padj_r', 'padj_absdiff', 'padj_shift']
        results_p['source'] = results_p['source'].astype(int).map(self.position_to_biom)
        results_p['target'] = results_p['target'].astype(int).map(self.position_to_biom)

        # print("Thresholding to get significant correlations")
        results_p.dropna(inplace=True)
        results_p.reset_index(drop=True)
        # In regular DC we remove self differential associations, unless these are specified in the source - target file
        if self.dracoon_program == 'DC':
            results_p = results_p[results_p['source'] != results_p['target']]
        self.res_unthresholded = results_p[["source", "target", "r_A", "padj_A", "r_B", "padj_B", "r", "padj_r",
                                            "absdiff", "p_absdiff", 'padj_absdiff', 'shift', 'p_shift', 'padj_shift']]
        os.system(f'rm -rf {tempfilename}')
        return self.res_unthresholded

    def threshold_results(self, df=None, corr_pval_filter=None):
        if df is not None:
            self.res_unthresholded = df
        if corr_pval_filter is not None:
            self.association_pvalue_filter = corr_pval_filter
        # Filtered by unadjusted q val
        self.res_p = self.threshold_pvalues(dataframe=self.res_unthresholded, column_prefix='p_', significance=1.1)
        # Filtered by adjusted q val
        self.res_padj = self.threshold_pvalues(dataframe=self.res_unthresholded, column_prefix='padj_', significance=self.significance)

        if self.association_pvalue_filter:
            # self.res_p = self.res_p[ (self.res_p.p_A < self.association_pvalue_filter)
            # | (self.res_p.p_B < self.association_pvalue_filter)]
            self.res_padj = self.res_padj[(self.res_padj.padj_A < self.association_pvalue_filter) | (
                    self.res_padj.padj_B < self.association_pvalue_filter)]

        self.res_unthresholded = self.res_unthresholded.rename(
            columns=lambda x: re.sub('_A$', '_' + str(self.conditions[0]), x))
        self.res_unthresholded = self.res_unthresholded.rename(
            columns=lambda x: re.sub('_B$', '_' + str(self.conditions[1]), x))

        self.res_p = self.res_p.rename(columns=lambda x: re.sub('_A$', '_' + str(self.conditions[0]), x))
        self.res_p = self.res_p.rename(columns=lambda x: re.sub('_B$', '_' + str(self.conditions[1]), x))

        self.res_padj = self.res_padj.rename(columns=lambda x: re.sub('_A$', '_' + str(self.conditions[0]), x))
        self.res_padj = self.res_padj.rename(columns=lambda x: re.sub('_B$', '_' + str(self.conditions[1]), x))
        return self.res_p, self.res_padj


    @staticmethod
    def threshold_pvalues(dataframe, column_prefix, significance=0.05):
        absdiff = column_prefix + "absdiff"
        shift = column_prefix + "shift"
        # create a list of our conditions
        conditioned = dataframe.iloc[:, dataframe.columns.isin([absdiff, shift])] < significance

        for column in conditioned.columns:
            conditioned[column] = np.where(conditioned[column] == True, column.split('_')[1] + '_', '')
        conditioned['sig'] = conditioned.apply(lambda row: ''.join(row.values.astype(str)), axis=1).str[:-1]
        conditioned['sig'].replace('', 'not_sig', inplace=True)
        # create a list of the values we want to assign for each condition
        dataframe = pd.concat([dataframe, conditioned['sig']], axis=1)
        dataframe = dataframe[dataframe['sig'] != 'not_sig']
        column_order = ['source', 'target', 'r_A', 'padj_A', 'r_B', 'padj_B', 'r', 'padj_r', 'absdiff', absdiff,
                        'shift', shift, 'sig']
        dataframe = dataframe[column_order]
        return dataframe

    def corrcoef_pval(self, nparray, epsilon=np.finfo(float).eps):
        if self.association_measure == 'spearman':
            for col in range(nparray.shape[1]):
                ranked_gene = np.empty_like(nparray[:, col])
                ranked_gene[np.argsort(nparray[:, col])] = np.arange(len(nparray[:, col]))
                ranked_gene += 1
                nparray[:, col] = ranked_gene

        r = np.corrcoef(nparray.T)
        if self.association_measure == 'entropy':
            r -= epsilon
            r = -((1 + abs(r)) / 2 * np.log((1 + abs(r)) / 2) + (1 - abs(r)) / 2 * np.log((1 - abs(r)) / 2))

        # r[np.diag_indices(r.shape[0])] = 0
        rf = r[np.triu_indices(r.shape[0], 1)]
        df = nparray.shape[0] - 2
        ts = rf * rf * (df / (1 - rf * rf))
        pf = sc.betainc(0.5 * df, 0.5, df / (df + ts))
        # Then we estimate p-values which indicate how different is the correlation coefficient from 0.
        p = np.zeros(shape=r.shape)
        p[np.triu_indices(p.shape[0], 1)] = pf
        p[np.tril_indices(p.shape[0], -1)] = p.T[np.tril_indices(p.shape[0], -1)]
        p[np.diag_indices(p.shape[0])] = np.ones(p.shape[0])

        if self.dracoon_program in ['DR', 'DCR']:
            DR_indices = [(self.biom_to_position[k[0]], self.biom_to_position[k[1]]) for k in self.source_targets]
            rf = np.array(list(map(functools.partial(self.findpos, array=r), DR_indices)))
            pf = np.array(list(map(functools.partial(self.findpos, array=p), DR_indices)))

        return rf, pf

    def print_correlation_inlinepermutation(self, kombi, data_A, data_B, biom_data, minifun, fun):
        ra, rb, rall = fun(kombi, data_A, data_B, biom_data)
        # Compute differential metrics
        abs_dif = abs(ra - rb)
        shift = (ra + rb) / self.alpha_shift - rall
        # Compute simulated values for the differential parameters n_iter times, 10000 by default
        bg_abs, bg_shift = random_dist_genepair(gene_i_values_a=data_A[:, kombi[0]],
                                                gene_i_values_b=data_B[:, kombi[0]],
                                                gene_j_values_a=data_A[:, kombi[1]],
                                                gene_j_values_b=data_B[:, kombi[1]],
                                                am=minifun, niters=self.iters)
        bg_abs = np.sort(bg_abs)
        bg_shift = np.sort(bg_shift)

        abs_dif_pval = test_significance(simulated_results=bg_abs, observed_value=abs_dif, two_sided=False)
        shift_pval = test_significance(simulated_results=bg_shift, observed_value=shift, two_sided=True)

        return float(kombi[0]), float(kombi[1]), ra, rb, rall, abs_dif, abs_dif_pval, shift, shift_pval

    @staticmethod
    def estimate_pval_fitdist(value, dist, param, two_sided=False):
        cd = dist.cdf(value, *param)  # Cumulative distribution value for this x.
        if two_sided:
            return min(cd, 1 - cd) * 2
        else:
            return 1 - cd

    def estimate_bestdist_params(self, distribution, possible_distributions, dist_name, timeout=60):
        fit_dist = Fitter(distribution, distributions=possible_distributions, timeout=timeout)
        fit_dist.fit()
        # fit_dist.summary()
        if self.verbose:
            fit_dist.summary(plot=True)
            plt.title('Fitting for ' + dist_name + ' metric background distribution')
            plt.show()
        best_fit = fit_dist.summary(plot=False).index.to_list()[0]
        params = fit_dist.fitted_param[best_fit]
        return getattr(ss, best_fit), params

    @staticmethod
    def findpos(position, array):
        return array[position]

    @staticmethod
    def estimate_corrpval(r, corrpoints, epsilon=np.finfo(float).eps):
        rf = np.array(r) - epsilon
        df = corrpoints - 2
        ts = rf * rf * (df / (1 - rf * rf))
        pf = sc.betainc(0.5 * df, 0.5, df / (df + ts))
        return pf


'''
expression_data = pd.read_csv('/Users/fernando/Documents/Research/DRACOONpy/app/expression_example.csv', index_col=0)
condition_data = pd.read_csv('/Users/fernando/Documents/Research/DRACOONpy/demo/condition_example.csv', index_col=0)
structure = pd.read_csv('/Users/fernando/Documents/Research/DRACOONpy/demo/structure_example.csv', index_col=0)


draconet = DraCooN_GRN(biom_data=expression_data,
                      cond_data=condition_data,
                      significance=0.01,
                      association_measure='entropy',
                      pvalue_adjustment_method='fdr_bh',
                      dracoon_program='DR',
                      associations_df=structure,
                      association_pvalue_filter=None,
                      pval_method='fitted_background',
                      iters=1000,
                      verbose=True,
                      matrixform=True)
draconet.run()                    
self = draconet
'''
