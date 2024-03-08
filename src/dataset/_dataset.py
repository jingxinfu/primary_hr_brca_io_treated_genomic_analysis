from pathlib import Path
import pandas as pd
from ..settings import init_logger

logger = init_logger(__name__)


class Dataset:
    def __init__(self, files_name_map_dict, data_processed_folder: Path):
        self._files_name_map_dict = files_name_map_dict
        if isinstance(data_processed_folder, str):
            data_processed_folder = Path(data_processed_folder)
        assert data_processed_folder.is_dir(), f"Cannot find {data_processed_folder}."
        self.folder = data_processed_folder
        self._parse_available_files()

    @property
    def files_name_map_dict(self):
        return self._files_name_map_dict

    @files_name_map_dict.setter
    def files_name_map_dict(self, new_file_dict):
        self._files_name_map_dict = new_file_dict
        self._parse_available_files()

    def _parse_available_files(self):
        self.available_files = {
            name: self.folder / rel_path
            for name, rel_path in self.files_name_map_dict.items()
            if (self.folder / rel_path).is_file()
        }

    def load(self, name, *args, **kwargs):
        assert (
            name in self.available_files
        ), f"File {name} was not defined in the dataset."
        file_path = self.available_files[name]
        suffix = str(file_path).split('.')[-1]
        if suffix == 'csv':
            df = pd.read_csv(file_path, *args, **kwargs)
        elif suffix == 'tsv' or suffix == 'txt':
            df = pd.read_csv(file_path, sep='\t', *args, **kwargs)
        else:
            raise ValueError(f'unknown file type: {suffix}.')

        ## Manually operation on selected dataset.
        if name == 'TIMER2':
            logger.warn('Transpose the TIMER2 result.')
            df = df.T

        return df


class bulkDataset(Dataset):
    def __init__(self, data_processed_folder):
        # Default setting
        self.Default_Metadata_Loading = {
            ## Clinical and sequencing QC information
            "clinical_sequencing": {
                # change it to the manually corrected clinical metadata with IHC result
                "clinical_and_sequencing_metadata": "clin_IHC.csv",# "clinical_and_sequencing_metadata.csv",
            },
            ## WES Features
            "WES": {
                "mutation_signature": "SigProfilerAssignment_COSMIC_v3/mutation_signature.csv",
                "somatic_mutation_burden": 'somatic_mutation_burden.csv',
                "cga_wes_features": 'CGA_WES_features.csv',
            },
            ## bulk RNA-Seq Features
            'bulkRNA': {
                "medTIN":"bulkRNA_medTIN.csv", # rseqQC medTIN metrics
                "TIMER2": "bulkRNA_TIMER2.csv",  # immune cells infiltration
                "ssGSEA": "bulkRNA_ssGSEA.csv",  # KEGG pathway enrichment scores
                "RNA_Seq_Batch": "bulkRNA_Sequence_Batch.csv",  # sequence batch
                "PAM50": "bulkRNA_PAM50.csv", # PAM50 estimation from genefu function
                "Cytokine": "bulkRNA_Cytokine.csv",  # Cytokine estimation from CytoSig
            },
        }
        _files_name_map_dict = {
            # -------------------- WES Pipeline output
            "somatic_mutation": "somatic_mutation_data.csv",
            ## Somatic mutation calling statistics
            "MutSigCV": "MutSigCV2_sigGenes_with_PatientID/RP-2423_16466.tsv",
            # -------------------- RNA-Seq pipeline output
            "TPM": "bulkRNA_TPM.csv",
            "gene_counts": "bulkRNA_Count.csv",
            # -------------------- Downstream Analysis Result
            ## Somatic mutation
            "mutation_Response": "mutation_clinical_benefits_features.csv",
            ## RNA-Seq
            "DEGs_Response": "bulkRNA_clinical_benefits_features.csv",
            "Pathways_Response": "bulkRNA_clinical_benefits_pathways.csv",

            ## RNA-Seq without AC
            "NoAC_DEGs_Response": "bulkRNA_NoAC_clinical_benefits_features.csv",
            "NoAC_Pathways_Response": "bulkRNA_NoAC_clinical_benefits_pathways.csv",
            ## all genomic and transcriptomic features
            "Feature_Response": "DownstreamFeatures_clinical_benefits_features.csv",
        }
        for v in self.Default_Metadata_Loading.values():
            _files_name_map_dict.update(v)

        super().__init__(
            data_processed_folder=data_processed_folder,
            files_name_map_dict=_files_name_map_dict
        )
        self.metadata = self._getMetadata()

    def _getMetadata(self):
        metadata = []
        for tech, features in self.Default_Metadata_Loading.items():
            for feature_name in features.keys():
                logger.info(
                    f'Appending {feature_name} feature estimated from {tech} sequencing into metadata'
                )
                feature_df = self.load(feature_name, index_col=0)

                if tech != 'clinical_sequencing':
                    logger.info(
                        f'Attaching {tech} prefix to feature names in {feature_name}'
                    )
                    feature_df.columns = tech + '_' + feature_df.columns

                metadata.append(feature_df)

        logger.info("Merging metadata: Outer concat.")
        metadata = pd.concat(metadata, axis=1, join='outer')
        her2 = (metadata.Patient=='P25')
        logger.warn(
            f"Remove {her2.sum()} samples from P25 which is HER+"
            )
        metadata = metadata.loc[~her2,:]
        # non_rcbs = metadata['RCB'].isna()
        # if non_rcbs.sum() > 0:
        #     logger.warn(
        #         f"Remove {non_rcbs.sum()} samples without RCB information: {','.join(metadata.index[non_rcbs].tolist())}"
        #     )
        #     metadata = metadata.loc[~non_rcbs, :]
        logger.info(f"Fill NAs in series with object dtype by N/A.")
        for c in metadata.columns:
            if metadata[c].dtype == 'O':
                metadata[c].fillna('N/A', inplace=True)

        return metadata
