from typing import Union, Tuple, Dict, List
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
import statsmodels.stats.multitest as multi
import warnings
import subprocess
import tempfile
import gseapy as gp
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from mlxtend.evaluate import mcnemar
from ..settings import init_logger

HERE = Path(__file__).parent.parent

warnings.simplefilter('ignore', ConvergenceWarning)
warnings.simplefilter('ignore', RuntimeWarning)
logger = init_logger(__name__)


def wilcoxonTest(x: pd.Series, y: pd.Series):

    x1 = x[y == 0].dropna()
    x2 = x[y == 1].dropna()
    try:
        _, p =  stats.mannwhitneyu(x1, x2, method="exact")
    except ValueError:
        return np.nan, np.nan, 'Wilcoxon'
    delta = (x2.mean() - x1.mean()) / x.std()

    return delta, p, 'Wilcoxon'


def contigencyTableTest(x: pd.Series, y: pd.Series,paired:Union[pd.Series,None]=None,method:str='Fisher_exact'):
    if isinstance(paired,pd.Series):
        x1 = x[y==0]
        x1.index = x1.index.map(paired.to_dict())
        x2 = x[y==1]
        x2.index = x2.index.map(paired.to_dict())
        x1 = x1[x2.index] # make sure the pt ids is the same
        table = pd.crosstab(x1.values,x2.values,dropna=False).values
    else:
        table = pd.crosstab(x, y,dropna=False).values
    if table.shape != (2, 2):
        return np.nan, np.nan, 'fisher_exact'
    if method =='Fisher_exact':
        oddsr, p = stats.fisher_exact(table, alternative='two-sided')
    elif method =='McNemar_exact':
        # exact if sample size < 25
        exact = (table.sum().sum()<=30)
        oddsr, p =mcnemar(table,exact=exact, corrected=True)
        oddsr = f'{(x1==True).sum()}:{(x2==True).sum()}'
        if x.name == 'ERCC4':
            print(table)
    else:
        raise KeyError(f'Unknown {method}')

    return oddsr, p, method

def fdrAdjust(p_value, method='fdr_bh'):
    """Adjust p value
    Parameters
    ----------
    p_value : pd.Series
        Original p value
    method : str, optional
        Method used for testing and adjustment of pvalues. Default: "fdr_bh"
        Can be either the full name or initial letters. Available methods are:
        bonferroni : one-step correction
        sidak : one-step correction
        holm-sidak : step down method using Sidak adjustments
        holm : step-down method using Bonferroni adjustments
        simes-hochberg : step-up method (independent)
        hommel : closed method based on Simes tests (non-negative)
        fdr_bh : Benjamini/Hochberg (non-negative)
        fdr_by : Benjamini/Yekutieli (negative)
        fdr_tsbh : two stage fdr correction (non-negative)
        fdr_tsbky : two stage fdr correction (non-negative)
    Returns
    -------
    pd.Series
        p-values corrected for multiple tests
    """
    adjust_p = multi.multipletests(p_value, method=method)[1]
    return adjust_p


def glmFit(
    x: pd.Series,
    y: pd.Series,
    covars: Union[pd.Series, pd.DataFrame, None] = None,
    family: sm.families = sm.families.Binomial,
) -> Tuple[float, float]:
    """
    Using GLM regression to access the association of ``x`` with ``y`` while adjusting
    effects from ``covaris``

    Parameters
    ----------
    x : pd.Series
        dependent variables
    y : pd.Series
        independent variables
    covars : Union[pd.Series,pd.DataFrame,None], optional
        confounders, default: None.
    family : sm.families, optional
       Distribution of y. Each family can take a link instance as an argument.
       See statsmodels.family.family for more information., by default sm.families.Binomial()

    Returns
    -------
    Tuple[float,float]

        ``(tvalue,pvalue)`` for the variable x.

        - tvalue: Z score of x from the regression.
        - pvalue: P-value of x from the regression.
    """
    if covars is not None:
        X = pd.concat([x, covars], axis=1)
    else:
        X = x
    X = sm.add_constant(X)
    glm_m = sm.GLM(y, X, family=family(), missing='drop')
    method_name = f'GLM-{family.__name__}'
    try:
        result = glm_m.fit(disp=False)
    except:
        return np.nan, np.nan, method_name

    return result.tvalues[x.name], result.pvalues[x.name], method_name


def getDEGs(
    count: pd.DataFrame,
    tpm: pd.DataFrame,
    y: pd.Series,
    covars: Union[pd.Series, pd.DataFrame, None] = None,
) -> pd.DataFrame:
    """
    Perform DE analysis via Wilcoxon test, edgeR, DESeq2

    Parameters
    ----------
    count: pd.DataFrame
        gene x sample, raw count gene expression quantification matrix
    tpm: pd.DataFrame
        gene x sample, TPM gene expression quantification matrix
    y : pd.Series
        binary variable indicates the control (0) and experiment group (1).
    covars : Union[pd.Series,pd.DataFrame]
        confounding factors

    Returns
    -------
    pd.DataFrame
        DE analysis result, with columns: ``Gene, LFC, Pvalue, FDR, Method``
    """
    results = []
    method_dict = {
        "DESeq2": "d",
    }
    y = y[count.columns]
    tpm = tpm[count.columns]
    tpm.columns = y.values
    gene_mean = tpm.groupby(level=0, axis=1).mean()
    gene_mean.columns = gene_mean.columns.map(lambda x: f'AveTPM_Group{x:0.0f}')
    gene_mean.index.name = 'Gene'

    with tempfile.TemporaryDirectory() as tmpdirname:
        in_file = str(Path(tmpdirname) / 'gene.readCount.matrix.tsv')
        condition = str(Path(tmpdirname) / 'Conditions.txt')
        output = str(Path(tmpdirname) / 'result.tsv')

        # write count and condition files
        count.to_csv(in_file, sep='\t')
        y.to_frame().T.to_csv(condition, index=False, sep='\t', header=None)

        cmd = " ".join(
            [
                'Rscript',
                f"{HERE}/Rscripts/DEGs.R",
                f"-g={in_file}",
                f"-c={condition}",
                f"-o={output}",
                "-f=1",  # keep all DEG result for sake of the vocano plot
            ]
        )
        for method, abbr in method_dict.items():
            the_cmd = cmd + f' -s={abbr}'
            logger.info(f'Call: {the_cmd}')
            subprocess.call(
                the_cmd,
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            try:
                result = pd.read_csv(output, sep='\t')
            except pd.errors.EmptyDataError:
                continue

            result.rename(
                columns={
                    'baseMean': 'AveExpr',  # DESeq2
                    'log2FoldChange': 'log2FC',
                    'logFC': 'log2FC',
                    'P.Value': 'pvalue',
                    'PValue': 'pvalue',
                    'p-value': 'pvalue',
                    'adjPval': 'padj',
                    'adj.P.Val': 'padj',
                    'FDR': 'padj',  # edgeR
                    'rawPval': 'pvalue',  # dearseq
                },
                inplace=True,
            )

            if not 'log2FC' in result.columns:  # dearseq and wilcoxon-test
                result['log2FC'] = np.log2(
                    (gene_mean['AveTPM_Group1'] + 1e-5)
                    / (gene_mean['AveTPM_Group0'] + 1e-5)
                )[result.Gene].values

            if 'prob' in result.columns:  # NOISeq
                result['padj'] = 1 - result['prob']
                result['pvalue'] = np.nan
            used_columns = ['Gene', 'log2FC', 'pvalue', 'padj']
            result = result[used_columns].merge(
                gene_mean.reset_index(), on='Gene', how='left'
            )
            result['Method'] = method
            results.append(result)

    return pd.concat(results, axis=0).set_index('Gene').reset_index()


def getComparisonXAndY(
    X: pd.DataFrame,
    Y: pd.DataFrame,
    compare_column: str,
    compare_column_rename: Union[Dict[str, int], None] = None,
    subset_string: Union[Dict[str, List], None] = None,
    confounders: Union[List[str], None] = None,
    min_n_samples: int = 4,
    min_feature_variance: Union[float, None] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Generate X and Y variables for compare the features in X based on labels on Y

    Parameters
    ----------
    X : pd.DataFrame
        feature matrix want to compare (Sample x Feature).
    Y : pd.DataFrame
        metadata information data frame containing both response and confounding variables.
    compare_column:str
        name of the response column or the column wants to compare the X on.
    compare_column_rename:Dict[str,int]
        a dict to encode the str variable in the ``compare_column`` in to binary values.
    subset_string : Union[Dict[str, List],None], optional
        a dict to subset datasets. e.g. only include samples has ``value`` in the ``key`` column of ``Y``, by default None.
    confounders : Union[List[str], None], optional
        the confounding factors should be considered when performing the comparison. the column names in the ``Y``, by default None.
    min_n_samples : int, optional
        number of minimum required samples to perform the analysis, by default 4
    min_feature_variance: Union[float,None], optional
        features with at least the variance to be included in the analysis, by default None
    Returns
    -------
    Tuple[pd.DataFrame,pd.DataFrame]
        preprocessed X and Y variables.
    """
    samples = Y.index
    # avoid overrwrite
    Y = Y.copy()
    X = X.copy()
    assert compare_column in Y.columns, f"{compare_column} cannot be found in the Y."
    if confounders is not None:
        miss_cols = pd.Index(confounders).difference(Y.columns)
        assert miss_cols.size == 0, f"{','.join(miss_cols)} cannot be found in the Y."
    # --------------------Filtering out samples
    if subset_string is not None:
        miss_cols = pd.Index(subset_string.keys()).difference(Y.columns)
        assert miss_cols.size == 0, f"{','.join(miss_cols)} cannot be found in the Y."
        ## subsetting samples based on the subset_string
        for col, selected_groups in subset_string.items():
            samples = samples.intersection(Y.index[Y[col].isin(selected_groups)])
    ## convert the str variables to numeric
    if compare_column_rename is not None:
        logger.info(f"Rename {compare_column}: {compare_column_rename}")
        Y[compare_column] = Y[compare_column].map(compare_column_rename)

    ## remove samples with NAs in slected columns
    involve_cols = [compare_column]
    if confounders is not None:
        involve_cols += confounders
    Y_test = Y.loc[samples, involve_cols]
    # check NAs in the comparing and confounding variables
    na_stats = Y_test.isna().sum()
    if na_stats.sum() > 0:
        logger.warn(f'Remove rows with NAs in columns:\n {na_stats}')
    Y_test.dropna(inplace=True)
    assert (
        Y_test.shape[0] > min_n_samples
    ), f"Remaining number of samples < {min_n_samples}"
    # make the compared column categorical
    Y[compare_column] = pd.Categorical(Y[compare_column])
    classes = Y_test[compare_column].value_counts()
    assert classes.size > 1, f"Only one class found: {classes.size}."
    logger.info(
        f"Size of classes from {subset_string}: {','.join([ str(k)+':'+str(v) for k,v in classes.items()])}"
    )
    # --------------------Filtering out testing variables
    X_test = X.loc[Y_test.index, :]
    if min_feature_variance is not None:
        X_test = X_test.loc[:, X_test.var(axis=0) > min_feature_variance]
        assert (
            X_test.shape[1] > 0
        ), f"No feature with variance > {min_feature_variance}."
        logger.info(
            f"{X_test.shape[1]} features with variance > {min_feature_variance} left for downstream analysis."
        )

    return X_test, Y_test


class bulkAnalysis:
    def __init__(
        self,
        feature_mtrx: pd.DataFrame,
        metadata: pd.DataFrame,
        identity_col: Union[str, None] ='Patient',
        confounders: Union[List[str], None] = ['facets_purity', 'Nonsynonymous'],
        arms_of_interested: List[str] = [
            'Chemo->ICI',
            'ICI->Chemo',
            'Chemo->ICI,ICI->Chemo',
        ],
        response_of_interested: List[str] = ['0-I', 'II-III', '0-I,II-III'],
        timepoints_of_interested: List[str] = ['Baseline', 'W3D1', 'W7D1', 'Surg+AC'],
        response_map_dict: Dict[str, int] = {
            '0-I': 1,
            'II-III': 0,
        },
        min_n_samples: int = 4,
        significance_column: str = 'padj',
        significance_threshold: float = 0.05,
        statistics_column: str = 'Zscore',
        statistics_threshold: float = 0.0,
        min_feature_variance: Union[float, None] = None,
    ):
        """
        Class for performing downstream analysis for bulk sequencing data.

        Parameters
        ----------
        feature_mtrx : pd.DataFrame
            raw count, somatic mutation, or CNV profile. (sample x feature)
        metadata : pd.DataFrame
            index by the sample name. the sample metadata with at lease three columns ``['BestResponse', 'Timepoint','Treatment_Exposure']``
            - ``BestResponse``: should contains response information, where ``CR/PR`` represents response and ``SD/PD`` represents non-response
            - ``Timepoint``: indicates when the sample was captured.
            - ``Treatment_Exposure``: indicates the treatment the sample being exposed to.
        confounders : Union[List[str],None], optional
            confounders adjusted on the association analysis, by default ``['facets_purity','Nonsynonymous']``
        """
        ol_samples = feature_mtrx.index.intersection(metadata.index)
        self.feature_mtrx = feature_mtrx.loc[ol_samples, :]
        self.metadata = metadata.loc[ol_samples, :]
        self.confounders = confounders if confounders is not None else []
        self.timepoints_of_interested = timepoints_of_interested
        self.arms_of_interested = arms_of_interested
        self.response_of_interested = response_of_interested
        self.response_map_dict = response_map_dict
        self._check_data_format()
        self.min_n_samples = min_n_samples
        self.min_feature_variance = min_feature_variance
        self.significance_column = significance_column
        self.significance_threshold = significance_threshold
        self.statistics_column = statistics_column
        self.statistics_threshold = statistics_threshold
        self.identity_col = identity_col

    def _check_data_format(self):
        for c in ['BestResponse', 'Timepoint', 'Treatment_Arm']:
            assert (
                c in self.metadata.columns
            ), f"required variable {c} not found in the metadata."
        for c in self.confounders:
            assert (
                c in self.metadata.columns
            ), f"confounding variable {c} not found in the metadata."
        for c in self.timepoints_of_interested:
            assert (
                c in self.metadata['Timepoint'].unique().tolist()
            ), f"Timepoint {c} not found in the metadata."

        self.samples = self.feature_mtrx.index

    def test_function(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        covars: Union[pd.Series, pd.DataFrame, None],
        paired: Union[pd.Series,None]=None,
    ) -> pd.DataFrame:
        """
        Find the features in the X that associates with the

        Parameters
        ----------
        X : pd.DataFrame
            independent variables (sample,feature)
        y : pd.Series
            dependent variables (sample,)
        covars : Union[pd.Series,pd.DataFrame,None]
            confounding variables (sample,confounding variable)
        paired: Union[pd.Series,None], defualt: None
            indicate the paired test.
        Returns
        -------
        pd.DataFrame
            association analysis results with ``"Zscore","Pvalue","FDR"`` for each feature
        """
        if isinstance(paired,pd.Series):
            result = X.apply(lambda c: contigencyTableTest(x=c, y=y,method='McNemar_exact',paired=paired), axis=0).T
        elif (covars is None) and X[X.columns[0]].dtype == 'category':
            result = X.apply(lambda c: contigencyTableTest(x=c, y=y,method='Fisher_exact'), axis=0).T
        else:
            result = X.apply(lambda c: glmFit(x=c, y=y, covars=covars), axis=0).T
        result = result.apply(pd.Series)
        result.columns = ['Zscore', 'Pvalue', 'Method']
        result['padj'] = np.nan
        try:
            result.loc[~result.Pvalue.isna(), 'padj'] = fdrAdjust(
                p_value=result['Pvalue'][~result.Pvalue.isna()]
            )
        except:
            result['padj'] = np.nan

        return result.reset_index()

    def get_clinical_benefits_features(
        self,
        timepoint: Union[str, None] = 'Baseline',
        treatment_arm: Union[str, None] = 'Pembro->Abraxane',
        *args,
        **kwargs,
    ) -> pd.DataFrame:
        """
        Calculate the clinical benefits of genes harboring mutations

         Parameters
         ----------
         timepoint : Union[str,None], optional
             perform the analysis on samples captured at the timepoint , by default ``'Baseline'``
         treatment_arm: Union[str,None], optional
             perform the analysis on samples on the treatment arm, by default ``'Pembro->Abraxane'``
         min_n_samples : int, optional
             requires at least ``min_n_samples`` samples with the annotation, by default 4
         Returns
         -------
         pd.DataFrame
             testing result with four columns ``[category,Zscore,Pvalue,FDR]``
        """
        compare_column = "BestResponse"
        # construct the subset string
        subset_string = {}
        if timepoint is not None:
            subset_string.update({"Timepoint": [timepoint]})
        if treatment_arm is not None:
            subset_string.update({"Treatment_Arm": treatment_arm.split(',')})
        if len(subset_string) == 0:
            subset_string = None
        # construct the final X, Y
        X, Y = getComparisonXAndY(
            X=self.feature_mtrx,
            Y=self.metadata,
            compare_column=compare_column,
            compare_column_rename=self.response_map_dict,
            subset_string=subset_string,
            confounders=self.confounders,
            min_n_samples=self.min_n_samples,
            min_feature_variance=self.min_feature_variance,
        )
        y = Y[compare_column]
        covars = Y[self.confounders] if self.confounders is not None else None
        result = self.test_function(X=X, y=y, covars=covars, *args, **kwargs)
        result['Timepoint'] = timepoint
        result['Treatment_Arm'] = treatment_arm
        result['Comparison'] = ','.join(
            [f'{k}:{v}' for k, v in self.response_map_dict.items()]
        )
        result['N_Sample'] = ','.join([f'{k}:{v}' for k, v in y.value_counts().items()])
        return result

    def get_treatment_impacted_features(
        self,
        response: Union[str, None] = '0-I',
        treatment_arm: Union[str, None] = 'Pembro->Abraxane',
        *args,
        **kwargs,
    ) -> pd.DataFrame:
        """
        Calculate the features that might be impacted by the treatment exposure

        Parameters
        ----------
        response : Union[str,None], optional
            perform the analysis on samples that from responders/non-responders , by default ``'CR/PR'``
        treatment_arm: Union[str,None], optional
            perform the analysis on samples on the treatment arm, by default ``'Pembro->Abraxane'``
        Returns
        -------
        pd.DataFrame
            testing result with four columns ``[category,Zscore,Pvalue,FDR]``
        """
        compare_column = "Timepoint"
        results = []
        for i in range(len(self.timepoints_of_interested) - 1):
            for j in range(i + 1, len(self.timepoints_of_interested)):
                compare_column_rename = {
                    self.timepoints_of_interested[i]: 0,
                    self.timepoints_of_interested[j]: 1,
                }
                # construct the subset string
                subset_string = {"Timepoint": list(compare_column_rename.keys())}
                if response is not None:
                    subset_string.update({"BestResponse": response.split(',')})
                if treatment_arm is not None:
                    subset_string.update({"Treatment_Arm": treatment_arm.split(',')})
                if len(subset_string) == 0:
                    subset_string = None
                    
                if self.identity_col is not None:
                    logger.info(f'Perform paired analysis based on the {self.identity_col}.')
                    n_sample = self.metadata.loc[self.metadata[compare_column].isin(compare_column_rename.keys()),
                                               self.identity_col].value_counts()
                    overlap_ids = n_sample[n_sample>1].index # focus on overlapped samples
                    Y_select = self.metadata.loc[self.metadata[self.identity_col].isin(overlap_ids),:].copy()
                else:
                    Y_select = self.metadata.copy()
                try:
                    # construct the final X, Y
                    X, Y = getComparisonXAndY(
                        X=self.feature_mtrx,
                        Y=Y_select,
                        compare_column=compare_column,
                        compare_column_rename=compare_column_rename,
                        subset_string=subset_string,
                        confounders=self.confounders,
                        min_n_samples=self.min_n_samples,
                        min_feature_variance=self.min_feature_variance,
                    )
                    y = Y[compare_column]
                    covars = (
                        Y[self.confounders] if self.confounders is not None else None
                    )
                    result = self.test_function(
                        X=X, y=y, covars=covars,paired= Y_select[self.identity_col],*args, **kwargs
                    )
                except AssertionError as e:
                    logger.warn(e)
                    continue
                result['Comparison'] = ','.join(
                    [f'{k}:{v}' for k, v in compare_column_rename.items()]
                )
                result['N_Sample'] = ','.join(
                    [f'{k}:{v}' for k, v in y.value_counts().items()]
                )
                result['BestResponse'] = response
                result['Treatment_Arm'] = treatment_arm

                results.append(result)
        if len(results) == 0:
            return pd.DataFrame()
        return pd.concat(results, axis=0, ignore_index=True)

    def analyze_arm_response_difference(
        self,
        pathway_analysis=False,
        output_prefix: Union[str, None] = None,
        *args,
        **kwargs,
    ):
        clin_benefits_features = [pd.DataFrame(), pd.DataFrame()]
        clin_benefits_pathways = [pd.DataFrame(), pd.DataFrame()]

        for arm in self.arms_of_interested:
            for timepoint in self.timepoints_of_interested:
                logger.info(f'Analyzing {arm} {timepoint}.')
                # single feature analysis
                try:
                    result = self.get_clinical_benefits_features(
                        timepoint=timepoint, treatment_arm=arm, *args, **kwargs
                    )
                except AssertionError as e:
                    logger.warn(e)
                    continue
                clin_benefits_features.append(result)

                if not pathway_analysis:
                    continue
                # enrichment analysis
                try:
                    pathway = self.get_enriched_pathways(
                        df=result, test_name=f"{timepoint}_{arm}"
                    )
                    pathway["Treatment_Arm"] = arm
                    pathway['Timepoint'] = timepoint

                except AssertionError as e:
                    logger.warn(e)
                    continue

                clin_benefits_pathways.append(pathway)
        if output_prefix is None:
            return pd.concat(
                clin_benefits_features, axis=0, ignore_index=True
            ), pd.concat(clin_benefits_pathways, axis=0, ignore_index=True)
        else:
            pd.concat(clin_benefits_features, axis=0, ignore_index=True).to_csv(
                f"{output_prefix}_clinical_benefits_features.csv", index=False
            )
            pd.concat(clin_benefits_pathways, axis=0, ignore_index=True).to_csv(
                f"{output_prefix}_clinical_benefits_pathways.csv", index=False
            )

    def analyze_treatment_impact(
        self,
        pathway_analysis=False,
        output_prefix: Union[str, None] = None,
        *args,
        **kwargs,
    ):
        impacts_features = [pd.DataFrame(), pd.DataFrame()]
        impacts_pathways = [pd.DataFrame(), pd.DataFrame()]
        for arm in self.arms_of_interested:
            for response in self.response_of_interested:
                logger.info(f'Analyzing {arm} {response}.')
                # single feature analysis
                result = self.get_treatment_impacted_features(
                    response=response, treatment_arm=arm, *args, **kwargs
                )
                impacts_features.append(result)

                # enrichment analysis
                if not pathway_analysis:
                    continue
                try:
                    pathway = self.get_enriched_pathways(
                        df=result, test_name=f"{arm}_{response}"
                    )
                    pathway["Treatment_Arm"] = arm
                    pathway['BestResponse'] = response

                except AssertionError as e:
                    logger.warn(e)
                    continue
                impacts_pathways.append(pathway)

        if output_prefix is None:
            return pd.concat(impacts_features, axis=0, ignore_index=True), pd.concat(
                impacts_pathways, axis=0, ignore_index=True
            )
        else:
            pd.concat(impacts_features, axis=0, ignore_index=True).to_csv(
                f'{output_prefix}_treatment_impact_features.csv', index=False
            )
            pd.concat(impacts_pathways, axis=0, ignore_index=True).to_csv(
                f'{output_prefix}_treatment_impact_pathways.csv', index=False
            )

    def get_enriched_pathways(self, df: pd.DataFrame, test_name: str) -> pd.DataFrame:
        assert "Gene" in df.columns, "Cannot find the Gene name column."
        results = [pd.DataFrame(), pd.DataFrame()]
        logger.warn(
            f'used {self.significance_column} < {self.significance_threshold} for DEGs detetion.'
        )
        # filtered_df = df.loc[df[self.significance_column]<self.significance_threshold,:]
        filtered_df = df.loc[
            (df[self.significance_column] < self.significance_threshold)
            & (  # significance
                df[self.statistics_column].abs() > self.statistics_threshold
            ),
            :,  # minimum changes
        ].copy()
        filtered_df['Direction'] = 'Up'
        filtered_df.loc[filtered_df[self.statistics_column] < 0, 'Direction'] = 'Down'
        for (comparison, method, direction), sub_df in filtered_df.groupby(
            ['Comparison', 'Method', 'Direction']
        ):
            gene_list = sub_df['Gene'].tolist()
            logger.info(
                f'[{test_name}]Start {comparison} for {len(gene_list)} detected by {method}.'
            )
            if len(gene_list) < 10:
                continue

            try:
                enr = gp.enrichr(
                    gene_list=gene_list,
                    gene_sets=['MSigDB_Hallmark_2020'],
                    organism='Human',
                    outdir='reports/enrichr_kegg/{comparison}_{method}_{test_name}',
                    no_plot=False,
                    cutoff=0.5,  # This argument control the terms (e.g FDR < 0.05) that will be shown on figures, not the result table output.
                )
                print(enr.results.head())
            except Exception as e:
                logger.warn(e)
                continue
            logger.info(f'Finish {comparison} for {method}: {test_name}')
            result = enr.results
            # result=pre_res.res2d
            result['Method'] = method
            result['Test_name'] = test_name
            result['Direction'] = direction
            result['Comparison'] = comparison
            results.append(result)

        return pd.concat(results, axis=0, ignore_index=True)