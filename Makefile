.PHONY: clean data lint

#################################################################################
# GLOBALS                                                                       #
#################################################################################
PROJECT_NAME = hr_brca_16-466
SHELL=/bin/bash
CONDAROOT = $$(conda info --base)
#################################################################################
# COMMANDS                                                                      #
#################################################################################

## Delete all compiled Python files
clean:
	find . -type f -name "*.py[co]" -delete
	find . -type d -name "__pycache__" -delete

develop: clean
	pip install -e '.[dev,doc,test]'

## Lint using black
lint:
	black src

## Set up python interpreter environment
create_environment:
	@echo ">>> Detected conda, creating conda environment."
	mamba env create -f conda-env.yml
	@echo ">>> New conda env created. Activate with:\nsource activate $(PROJECT_NAME)"
	source $(CONDAROOT)/bin/activate $(PROJECT_NAME) && \
	pip install -r requirements.txt && \
	Rscript INSTALL.R

#################################################################################
# PROJECT RULES                                                                 #
#################################################################################

## Make WES Dataset
wes_data: #develop
	make_wes_dataset data/raw data/processed
	echo "Perform mutational signature analysis by deconstructSig on COSMIC Version 3"
	Rscript src/Rscripts/deconstructSig.R \
			--maf_path data/processed/somatic_mutation_merged.maf \
			--signature_path src/external/COSMIC_v3.4_SBS_GRCh37.txt \
			--output_folder data/processed/deconstructSig_COSMIC_V3



## Make bulk RNA-Seq Dataset
rna_data: #develop
	make_rna_dataset data/raw data/processed
	Rscript src/Rscripts/pam50.R

## Make all data
data: wes_data rna_data

## Analyze mutation data from the processed data folder
analyze_mutation:
	analyze_mutation data/processed

## Generate Comut plot
gen_comut:
	analyze_landscape data/processed 

## Analyze bulk RNA-seq data from the processed data folder
analyze_rna:
	analyze_rna data/processed
## Analyze genomic and transcriptomic features from the processed data folder
analyze_feature:
	analyze_feature data/processed

#################################################################################
# Self Documenting Commands                                                     #
#################################################################################

.DEFAULT_GOAL := help

# Inspired by <http://marmelab.com/blog/2016/02/29/auto-documented-makefile.html>
# sed script explained:
# /^##/:
# 	* save line in hold space
# 	* purge line
# 	* Loop:
# 		* append newline + line to hold space
# 		* go to next line
# 		* if line starts with doc comment, strip comment character off and loop
# 	* remove target prerequisites
# 	* append hold space (+ newline) to line
# 	* replace newline plus comments by `---`
# 	* print line
# Separate expressions are necessary because labels cannot be delimited by
# semicolon; see <http://stackoverflow.com/a/11799865/1968>
.PHONY: help
help:
	@echo "$$(tput bold)Available rules:$$(tput sgr0)"
	@echo
	@sed -n -e "/^## / { \
		h; \
		s/.*//; \
		:doc" \
		-e "H; \
		n; \
		s/^## //; \
		t doc" \
		-e "s/:.*//; \
		G; \
		s/\\n## /---/; \
		s/\\n/ /g; \
		p; \
	}" ${MAKEFILE_LIST} \
	| LC_ALL='C' sort --ignore-case \
	| awk -F '---' \
		-v ncol=$$(tput cols) \
		-v indent=19 \
		-v col_on="$$(tput setaf 6)" \
		-v col_off="$$(tput sgr0)" \
	'{ \
		printf "%s%*s%s ", col_on, -indent, $$1, col_off; \
		n = split($$2, words, " "); \
		line_length = ncol - indent; \
		for (i = 1; i <= n; i++) { \
			line_length -= length(words[i]) + 1; \
			if (line_length <= 0) { \
				line_length = ncol - indent - length(words[i]) - 1; \
				printf "\n%*s ", -indent, " "; \
			} \
			printf "%s ", words[i]; \
		} \
		printf "\n"; \
	}' \
	| more $(shell test $(shell uname) = Darwin && echo '--no-init --raw-control-chars')
