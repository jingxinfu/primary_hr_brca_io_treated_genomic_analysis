# -*- coding: utf-8 -*-
from typing import Tuple
import os
import click
import glob
from pathlib import Path
import numpy as np
import pandas as pd
from ..settings import init_logger,HERE
from . import Bulk_Correct_Sample_Name

logger = init_logger(__name__)
def pam50_normalization(expr:pd.DataFrame,ihc_classes:pd.DataFrame)->pd.DataFrame:
    """
    normalize the expression data for the PAM50 subtyping

    Parameters
    ----------
    expr : pd.DataFrame
        log normalized gene expression matrix (n_genes, n_samples)
    ihc_classes : pd.Series
        IHC classes Immunohistochemistry subgroup of breast cancer to be subtyped (i.e., ERpos_HER2neg, HER2pos_ERneg, HER2pos_ERpos, TNBC)(n_samples,)
        (n_samples,)

    Returns
    -------
    pd.DataFrame
        PAM50 normalized expression data (n_genes, n_samples)
    """
    centered = []
    PAM50_sigma = pd.read_csv(HERE/'external/SIGMAS_UNC232_v4.0.txt',sep='\t',index_col=0)
    # For each class in the test set get the samples and apply the quantiles from the train set
    for class_ in ihc_classes['IHC_class'].unique():
        samples = ihc_classes.index[ihc_classes['IHC_class'] == class_]
        gene_expr_class = expr[samples]
        percentile_group = PAM50_sigma[class_]
        centered.append(quantile_centering(gene_expr_class, percentile_group))
    result = pd.concat(centered, sort=False, axis=1)
    return result

def quantile_centering(expr_matrix, gene_quantile):
    """Do row centering based on the quantile and IHC group.
    Reference: https://unclineberger.org/peroulab/wp-content/uploads/sites/1008/2022/02/JCO-2020-Subgroup-specific-gene-centering-method-AFM.zip
    
    :param expr_matrix: pandas.DataFrame where row are genes and columns are samples
    :param gene_quantile: pandas.DataFrame or Series containig all the genes in the
      first parameter and the value of the quantile to be used, i.e. .5 if one wants
      to do row centering using the mean.

    :return: Subgroup-specific Centered dataframe.
    """
    gene_renames = {
         "CDCA1":"NUF2",
         "KNTC2": "NDC80",
         "ORC6L":"ORC6"}
    res = []
    for name in gene_quantile.index:
        q = gene_quantile.loc[name]
        if name in gene_renames:
            name = gene_renames[name]
        q_value = expr_matrix.loc[name].quantile(q)
        res.append(expr_matrix.loc[name] - q_value)
    res = pd.concat(res,axis=1,sort=False).T
    res.index = gene_quantile.index
    return res

def mergeTIN(folder_path:Path)-> pd.DataFrame:

    """Merge single sample's rseqQC TIN files into cohort level gene by sample data frame

    Parameters
    ----------
    folder_path : str
        Path to the cohort folder, where all salmon quantification output files locate

    Returns
    -------
    pd.DataFrame
        1. gene (Hugo Symbol) by sample Raw Count data frame
        2. gene (Hugo Symbol) by sample TPM data frame
    """
    tin_suffix = 'tsv'
    tins = []
    for file_path in glob.glob(str(folder_path / "*")):
        if not file_path.endswith(tin_suffix):
            continue
        sample_name = '_'.join(os.path.basename(file_path).split('_')[:4])
        # Correct the bulk sample name.
        if sample_name in Bulk_Correct_Sample_Name:
            sample_name = Bulk_Correct_Sample_Name[sample_name]
        df = pd.read_csv(file_path, sep="\t")
        tin = df['TIN(median)']
        tin.index = [sample_name]
        tins.append(tin)
    tins = pd.concat(tins, axis=0).to_frame()
    tins.index.name = 'Sample'
    return tins
        
def mergeSalmonGeneOutput(folder_path: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Merge single sample's salmon gene quantification files into cohort level gene by sample data frame

    Parameters
    ----------
    folder_path : str
        Path to the cohort folder, where all salmon quantification output files locate

    Returns
    -------
    (pd.DataFrame,pd.DataFrame)
        1. gene (Hugo Symbol) by sample Raw Count data frame
        2. gene (Hugo Symbol) by sample TPM data frame

    """
    salmon_gene_suffix = ".sf.gz"
    tpms = []
    counts = []
    batch_infos = []
    for file_path in glob.glob(str(folder_path / "*")):
        if not file_path.endswith(salmon_gene_suffix):
            continue
        sample_name = '_'.join(os.path.basename(file_path).split('_')[:4])
        # Correct the bulk sample name.
        if sample_name in Bulk_Correct_Sample_Name:
            sample_name = Bulk_Correct_Sample_Name[sample_name]
        df = pd.read_csv(file_path, sep="\t", index_col=0)
        tpm = df['TPM']
        tpm.name = sample_name
        tpms.append(tpm)

        count = df["NumReads"]
        count.name = sample_name
        counts.append(count.astype(int))

        batch_info = pd.Series(
            [os.path.basename(file_path).split('_')[5]],
            index=[sample_name],
            name='Sequence_Batch',
        )
        batch_infos.append(batch_info)

    return (
        pd.concat(counts, axis=1),
        pd.concat(tpms, axis=1),
        pd.concat(batch_infos, axis=0),
    )


@click.command()
@click.argument('input_filepath', type=click.Path(exists=True))
@click.argument('output_filepath', type=click.Path())
def main(input_filepath, output_filepath):
    """Runs data processing scripts to turn raw data from (../raw) into
    cleaned data ready to be analyzed (saved in ../processed).
    """
    input_filepath = Path(input_filepath)
    output_filepath = Path(output_filepath)
    tins = mergeTIN(folder_path=input_filepath/'rseQC_tin')
    tins.to_csv(output_filepath/'bulkRNA_medTIN.csv')

    counts, tpms, batchs = mergeSalmonGeneOutput(
        folder_path=input_filepath/ 'salmon_gene'
    )
    count_out, tpm_out, batch_out = (
        output_filepath / 'bulkRNA_Count.csv',
        output_filepath / 'bulkRNA_TPM.csv',
        output_filepath / 'bulkRNA_Sequence_Batch.csv',
    )
    logger.info(
        f"Store cohort level TPM , raw count matrix, and batch information into{(count_out, tpm_out,batch_out)}"
    )
    counts.to_csv(count_out)
    tpms.to_csv(tpm_out)
    batchs.to_csv(batch_out)

    logger.info('-------------------- PAM50 normalization --------------------')
    ihc_classes = pd.DataFrame(index=tpms.columns)
    ihc_classes['IHC_class'] = 'ERpos_HER2neg'
    pam50_normalized = pam50_normalization(expr=np.log2(tpms+1),ihc_classes=ihc_classes)
    pam50_normalized_out = output_filepath / 'bulkRNA_PAM50_Normalized_Expression.csv'
    logger.info(
        f"Store cohort level PAM50 normalized expression into{pam50_normalized_out}"
    )
    pam50_normalized.to_csv(pam50_normalized_out)

    logger.info('-------------------- Downstream analysis --------------------')
    logger.info('#-Perform CytoSig analysis.')
    import CytoSig

    signature = os.path.join(
        "src/external/CytoSig/CytoSig", 'signature.centroid'
    )  # load cytokine response signature installed in your python system path
    signature = pd.read_csv(signature, sep='\t', index_col=0)
    beta, std, zscore, pvalue = CytoSig.ridge_significance_test(
        signature,
        np.log(1 + tpms),
        alpha=1e4,
        alternative="two-sided",
        nrand=1000,
        cnt_thres=10,
        flag_normalize=True,
        verbose=True,
    )
    zscore.index = 'Cytokine_' + zscore.index
    zscore.T.to_csv(output_filepath / 'bulkRNA_Cytokine.csv')

    import gseapy as gp

    gene_sets = ["MSigDB_Hallmark_2020",
                 "data/external/gene_sets_for_antigen_presentation.gmt",
                 "data/external/gene_sets_for_estrogen.gmt"]
    ssgsea_result = pd.DataFrame()
    for gene_set in gene_sets:
        logger.info(
            f'#-Perform single sample GSEA analysis on {gene_set}--------------------'
        )
        geneset_result = gp.ssgsea(
            data=tpms,
            gene_sets=gene_set,
            outdir=None,
            sample_norm_method='rank',
            no_plot=True,
            processes=4,
        )
        ssgsea_result = pd.concat([ssgsea_result, geneset_result.res2d.pivot(index='Name', columns='Term', values='NES')], axis=1)
    ssgsea_result.columns = ssgsea_result.columns+'_ssGSEA'
    ssgsea_out = output_filepath / "bulkRNA_ssGSEA.csv"
    logger.info(f"Store cohort level ssGSEA result into {ssgsea_out}")
    ssgsea_result.to_csv(ssgsea_out)

    logger.info(
        '#-Perform immune cell deconvolution analysis using six public methods--------------------'
    )
    immunedeconv_out = output_filepath / 'bulkRNA_TIMER2.csv'
    logger.warn(
        f'Manually upload the {tpm_out} to the TIMER2 website and store the estimation result into {immunedeconv_out} '
    )


if __name__ == '__main__':
    main()
