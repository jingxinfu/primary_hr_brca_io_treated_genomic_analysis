# -*- coding: utf-8 -*-
from typing import List, Tuple
import os
import click
from comut import fileparsers
from pathlib import Path
import pandas as pd
from ..settings import init_logger,HERE
from . import Bulk_Correct_Sample_Name

logger = init_logger(__name__)

def summarize_WES_features(
    pair_table: pd.DataFrame, used_features: List[str], output_folder: Path
) -> None:
    """
    Extract needed information from Terra Pair Table

    Parameters
    ----------
    pair_table : pd.DataFrame
        Terra Pair Table
    used_features : List[str]
        feature names of interests
    output_folder : Path
        the output directory to store the file
    """
    df = pair_table[used_features].copy()

    df.index = (
        df.index.to_series()
        .map(lambda x: '_'.join(x.split('_')[:4]))
        .replace(Bulk_Correct_Sample_Name)
    )
    df.index.name = 'sample'
    df.to_csv(output_folder / 'CGA_WES_features.csv')


def summarize_somatic_mutation(
    maf_out_dir: Path,
    maf_for_tmb_out_dir: Path,
    base_coverage: pd.Series,
    output_folder: Path,
) -> None:
    """
    Summary mutation statistics from MAF files stored on ``maf_out_dir``

    Parameters
    ----------
    maf_out_dir : Path
        Folder path where stores MAF files
    maf_for_tmb_out_dir : Path
        Folder path where stores MAF files with CCF information.
    base_coverage: pd.Series
        the base coverage information
    output_folder: Path
        the output directory to store all three files

    Returns
    -------
    None

        - ``mutation_data``:mutation_data, with columns ``['sample','category','value']``
        - ``mutation_burden``: mutation_burden,with columns ``['sample','Nonsynonymous','Synonymous','Clonal','Subclonal']``
    """

    def get_tmb(
        maf_path: Path,
    ) -> pd.DataFrame:
        """
        Get Tumor mutation burden

        Parameters
        ----------
        maf_path : Path
            MAF file contains somatic mutation calling from the single sample
        Returns
        -------
        pd.DataFrame
            mutation burden table
        """
        # Synonymous_Classification = ['Silent']
        Nonsynonymous_Classification = [
            'Frame_Shift_Del',
            'Frame_Shift_Ins',
            'In_Frame_Del',
            'In_Frame_Ins',
            'Missense_Mutation',
            'Nonsense_Mutation',
            'Splice_Site',
            'Nonstop_Mutation',
        ]
        selected_maf_cols = [
            'Hugo_Symbol',
            'Entrez_Gene_Id',
            'Chromosome',
            'Start_position',
            'End_position',
            'Strand',
            'cDNA_Change',
            'Codon_Change',
            'Protein_Change',
            'Reference_Allele',
            'Tumor_Seq_Allele1',
            'Tumor_Seq_Allele2',
            'Variant_Type',
            'Variant_Classification',
            'tumor_f',
            'CCF_lower',  #  lower bound of the 95% confidence interval
            'clonal.ix',
            'subclonal.ix',
            'ccf_hat'
        ]
        mutation_data = pd.read_csv(
            maf_path,
            comment='#',
            sep='\t',
            encoding='unicode_escape',
        )
        mutation_data = mutation_data.loc[
            :, mutation_data.columns.intersection(selected_maf_cols)
        ]
        sample_name = os.path.basename(maf_path).split('.')[0]
        mutation_data['Tumor_Sample_Barcode'] = sample_name
        # NOTE: must have reads support the variant, ref: Brendan
        ## prior papers from the lab have calculated TMB using the force called MAF,
        ## you just have to make sure that you are  removing variants with zero supporting reads for each sample
        mutation_data = mutation_data.loc[mutation_data['tumor_f'] > 0, :]
        mutation_data['substitution'] = (
            mutation_data['Variant_Classification']
            .isin(Nonsynonymous_Classification)
            .astype(int)
            .map({1: 'Nonsynonymous', 0: 'Synonymous'})
        )
        nonsyn_mutation = (
            mutation_data.loc[mutation_data['substitution'] == 'Nonsynonymous', :]
            .dropna(subset=['ccf_hat'])
            .copy()
        )
        nonsyn_mutation['clonality'] = ''
        nonsyn_mutation.loc[nonsyn_mutation['clonal.ix']==1.0, 'clonality'] = 'Clonal'
        nonsyn_mutation.loc[nonsyn_mutation['subclonal.ix']==1.0, 'clonality'] = 'Subclonal'

        if all(
            nonsyn_mutation.loc[:, 'clonality'] == ''
        ):  # for MAF that does NOT contain the CCF information
            clonality = pd.DataFrame({'Tumor_Sample_Barcode': [sample_name]})
        else:
            clonality = pd.crosstab(
                nonsyn_mutation['Tumor_Sample_Barcode'], nonsyn_mutation['clonality']
            ).reset_index()
        clonality['All_Somatic'] = nonsyn_mutation.shape[0]

        for c in ['Clonal', 'Subclonal']:
            if not c in clonality.columns:
                clonality[c] = 0

        mutation_burden = pd.crosstab(
            mutation_data['Tumor_Sample_Barcode'], mutation_data['substitution']
        ).reset_index()

        mutation_burden = pd.merge(
            clonality, mutation_burden, on='Tumor_Sample_Barcode', how='right'
        )

        mutation_burden.rename(columns={'Tumor_Sample_Barcode': 'sample'}, inplace=True)
        return mutation_burden

    def get_mutation_data_and_reformat_maf(
        maf_path: Path,
    ) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        reformat MAF file and get mutation data

        Parameters
        ----------
        maf_path : Path
            MAF file contains somatic mutation calling from the single sample
        Returns
        -------
        Tuple[pd.DataFrame,pd.DataFrame,pd.DataFrame]
            comut parsed format table, mutation burden table, MAF format table

        """
        # Synonymous_Classification = ['Silent']
        Nonsynonymous_Classification = [
            'Frame_Shift_Del',
            'Frame_Shift_Ins',
            'In_Frame_Del',
            'In_Frame_Ins',
            'Missense_Mutation',
            'Nonsense_Mutation',
            'Splice_Site',
            'Nonstop_Mutation',
        ]
        selected_maf_cols = [
            'Hugo_Symbol',
            'Entrez_Gene_Id',
            'Center',
            'NCBI_Build',
            'Chromosome',
            'Start_position',
            'End_position',
            'Strand',
            "Variant_Classification",
            "Variant_Type",
            "Reference_Allele",
            "Tumor_Seq_Allele1",
            "Tumor_Seq_Allele2",
            "dbSNP_RS",
            "dbSNP_Val_Status",
            "Tumor_Sample_Barcode", # the first 16 columns are required for sigProfiler Assignment
            'cDNA_Change',
            'Codon_Change',
            'Protein_Change',
            'tumor_f',
            'CCF_lower',  #  lower bound of the 95% confidence interval
            'subclonal.ix',
        ]
        mutation_data = pd.read_csv(
            maf_path,
            comment='#',
            sep='\t',
            encoding='unicode_escape',
        )
        mutation_data = mutation_data.loc[
            :, mutation_data.columns.intersection(selected_maf_cols)
        ]
        sample_name = os.path.basename(maf_path).split('.')[0]
        mutation_data['Tumor_Sample_Barcode'] = sample_name
        # NOTE: must have reads support the variant, ref: Brendan
        ## prior papers from the lab have calculated TMB using the force called MAF,
        ## you just have to make sure that you are  removing variants with zero supporting reads for each sample
        if 'tumor_f' in mutation_data.columns:
            mutation_data = mutation_data.loc[mutation_data['tumor_f'] > 0, :]
        mutation_df = fileparsers.parse_maf(mutation_data) # for comut plot
        mutation_data['substitution'] = (
            mutation_data['Variant_Classification']
            .isin(Nonsynonymous_Classification)
            .astype(int)
            .map({1: 'Nonsynonymous', 0: 'Synonymous'})
        )

        return mutation_df, mutation_data

    df_list = [
        get_mutation_data_and_reformat_maf(maf_path)
        for maf_path in maf_out_dir.glob('*.maf')
    ]
    tmb_list = [get_tmb(maf_path) for maf_path in maf_for_tmb_out_dir.glob('*.maf')]
    mutation_df = pd.concat([x[0] for x in df_list])
    ## store reformed maf file
    output_maf_folder = output_folder / 'reformat_maf'
    if not output_maf_folder.is_dir():
        output_maf_folder.mkdir()
    for _,formated_maf in df_list:
        sample_name = formated_maf.Tumor_Sample_Barcode.map(
                lambda x: '_'.join(x.split('_')[:4])
            ).replace(Bulk_Correct_Sample_Name).unique()[0]
        formated_maf.Tumor_Sample_Barcode = sample_name
        formated_maf.to_csv(output_maf_folder/ f"{sample_name}.maf",sep='\t',index=False)
    ## perform sigProfiler analysis
    merged_maf = pd.concat([x[1] for x in df_list])
    mutation_burden = pd.concat(tmb_list).set_index('sample')
    heterogenity = mutation_burden['Subclonal'] / mutation_burden['All_Somatic']

    mutation_burden.columns = mutation_burden.columns + '_TMB'

    mutation_burden = (
        mutation_burden.T * 1e6 / base_coverage[mutation_burden.index]
    ).T
    mutation_burden['Heterogeneity'] = heterogenity
    mutation_burden = mutation_burden.reset_index()

    ## Replace sample name
    merged_maf['Tumor_Sample_Barcode'] = (
        merged_maf['Tumor_Sample_Barcode']
        .map(lambda x: '_'.join(x.split('_')[:4]))
        .replace(Bulk_Correct_Sample_Name)
    )
    mutation_df['sample'] = (
        mutation_df['sample']
        .map(lambda x: '_'.join(x.split('_')[:4]))
        .replace(Bulk_Correct_Sample_Name)
    )
    mutation_burden['sample'] = (
        mutation_burden['sample']
        .map(lambda x: '_'.join(x.split('_')[:4]))
        .replace(Bulk_Correct_Sample_Name)
    )

    for df, name, sep in [
        (mutation_burden, 'somatic_mutation_burden.csv', ','),
        (mutation_df, 'somatic_mutation_data.csv', ','),
        (merged_maf, 'somatic_mutation_merged.maf', '\t'),
    ]:
        logger.info(f"Write {name.split('.')[0]} into {output_folder/name}")
        df.to_csv(output_folder / name, index=False, sep=sep)


@click.command()
@click.argument('input_filepath', type=click.Path(exists=True))
@click.argument('output_filepath', type=click.Path())
def main(input_filepath, output_filepath):
    """Runs data processing scripts to turn raw data from (../raw) into
    cleaned data ready to be analyzed (saved in ../processed).
    """
    used_features = [
        'combined_fracContam',
        'combined_TiN',
        'absolute_ploidy',
        'absolute_purity',
        'facets_wgd_bool'
    ]
    used_somatic_maf_col = 'combined_mutation_validator_oncotated_maf'
    used_somatic_maf_for_tmb_col = 'absolute_annotated_maf_capture'
    used_base_coverage = 'combined_mutect1_basescovered'
    output_filepath = Path(output_filepath)
    input_filepath = Path(input_filepath)
    pair_table = pd.read_csv(f"{input_filepath}/pair.tsv",sep='\t',index_col=0)
   
    
    logger.info(f"-Summarizing WES calling features:{used_features}")

    summarize_WES_features(
        pair_table=pair_table,
        used_features=used_features,
        output_folder=output_filepath,
    )

    logger.info('-Summarizing somatic mutation and generating files:')
    summarize_somatic_mutation(
        maf_out_dir=input_filepath / used_somatic_maf_col,
        maf_for_tmb_out_dir=input_filepath / used_somatic_maf_for_tmb_col,
        base_coverage=pair_table[used_base_coverage],
        output_folder=output_filepath,
    )
    exclude_signature_subgroups =None
    logger.info(f'-Perform mutational signature analysis [exclude_signature:{exclude_signature_subgroups}]:')
    
    logger.info('-----Using COSMIC Version 2:')
    from SigProfilerAssignment import Analyzer as Analyze
    Analyze.cosmic_fit(samples=f"{output_filepath}/reformat_maf/", output=f"{output_filepath}/SigProfilerAssignment_COSMIC_v2", input_type="vcf", context_type="96",
                   collapse_to_SBS96=False, cosmic_version=2,exome=True,
                   genome_build="GRCh37", signature_database=None,
                   exclude_signature_subgroups=exclude_signature_subgroups, export_probabilities=True,
                   export_probabilities_per_mutation=False, make_plots=True,sample_reconstruction_plots='png',
                   verbose=False)
    ## the aetiology reference was extracted from maftools
    aetiology_v2 = pd.read_csv( HERE / 'external/cosmic_v2_aetiology.csv',index_col=0)
    aetiology_v2.index = aetiology_v2.index.str.replace('COSMIC','Signature')
    aetiology_v2 = aetiology_v2['aetiology'].str.capitalize().to_dict()
    sigscore = pd.read_csv(output_filepath / 'SigProfilerAssignment_COSMIC_v2/Assignment_Solution/Activities/Assignment_Solution_Activities.txt',
                        sep='\t',index_col=0)
    sigs_n_sample = sigscore.sum()
    nonzero_sigs = sigs_n_sample.index[sigs_n_sample>0]
    sigscore = sigscore[nonzero_sigs].T
    sigscore = (sigscore / sigscore.sum()).T
    ## store dominant signature
    sigmap = sigscore.idxmax(axis=1).rename('Sig').to_frame()
    sigmap['aetiology'] = sigmap.Sig.map(aetiology_v2)
    sigmap.to_csv(output_filepath/'SigProfilerAssignment_COSMIC_v2/dominant_mutation_signature.csv')
    ## store signature score
    sigscore.to_csv(output_filepath/'SigProfilerAssignment_COSMIC_v2/mutation_signature.csv')
    
    logger.info('-----Using COSMIC Version 3:')
    Analyze.cosmic_fit(samples=f"{output_filepath}/reformat_maf/", output=f"{output_filepath}/SigProfilerAssignment_COSMIC_v3", input_type="vcf", context_type="96",
                   collapse_to_SBS96=False, cosmic_version=3.4,exome=True,
                   genome_build="GRCh37", signature_database=None,
                   exclude_signature_subgroups=exclude_signature_subgroups, export_probabilities=True,
                   export_probabilities_per_mutation=False, make_plots=True,sample_reconstruction_plots='png',
                   verbose=False)
    ## the aetiology reference was extracted from maftools
    aetiology_v3 = pd.read_csv( HERE / 'external/cosmic_v3_aetiology.csv',index_col=0)['aetiology'].str.capitalize().to_dict()
    sigscore = pd.read_csv(output_filepath / 'SigProfilerAssignment_COSMIC_v3/Assignment_Solution/Activities/Assignment_Solution_Activities.txt',
                        sep='\t',index_col=0)
    sigs_n_sample = sigscore.sum()
    nonzero_sigs = sigs_n_sample.index[sigs_n_sample>0]
    sigscore = sigscore[nonzero_sigs].T
    sigscore = (sigscore / sigscore.sum()).T
    ## store dominant signature
    sigmap = sigscore.idxmax(axis=1).rename('Sig').to_frame()
    sigmap['aetiology'] = sigmap.Sig.map(aetiology_v3)
    sigmap.to_csv(output_filepath/'SigProfilerAssignment_COSMIC_v3/dominant_mutation_signature.csv')
    ## store signature score
    sigscore.to_csv(output_filepath/'SigProfilerAssignment_COSMIC_v3/mutation_signature.csv')


if __name__ == '__main__':
    main()
