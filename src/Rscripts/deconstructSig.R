library(argparse, quiet = T)
library(deconstructSigs, quiet = T)
library(plyr, quiet = T)
library(data.table)
DESCRIPTION <- "Mutational Signature Analysis by DeconstructSig"
MutSigEstimate <- function(maf_path, signature_path) {
    ref <- data.frame(t(read.table(signature_path, sep = "\t", row.names = 1, header = T, check.names = FALSE)), check.names = FALSE)
    snp.maf <- data.frame(fread(maf_path, sep = "\t", header = T)) # read.table cannot handle table with too many NAs. Don't know how it filsl values, but it mixup lines.
    snp.maf[, "Chromosome"] <- as.character(snp.maf[, "Chromosome"])
    snp.maf[, "Chromosome"] <- sub("^([0-9XY])", "chr\\1", snp.maf[, "Chromosome"])
    snp.maf[, "Chromosome"] <- sub("^MT", "chrM", snp.maf[, "Chromosome"])
    snp.maf[, "Chromosome"] <- sub("^(GL[0-9]+).[0-9]", "chrUn_\\L\\1", snp.maf[, "Chromosome"], perl = T)

    sigs.input <- deconstructSigs::mut.to.sigs.input(
        mut.ref = snp.maf,
        sample.id = "Tumor_Sample_Barcode",
        chr = "Chromosome",
        pos = "Start_position",
        ref = "Reference_Allele",
        alt = "Tumor_Seq_Allele2"
    )
    ref <- ref[, colnames(sigs.input)] # align column names
    sig.contrib.list <- lapply(row.names(sigs.input), function(x) {
        row.df <- deconstructSigs::whichSignatures(
            tumor.ref = sigs.input,
            sample.id = x,
            contexts.needed = TRUE, # normalization on tumor.ref by tri.counts.method
            tri.counts.method = "exome2genome", # normalization is performed to reflect the absolute frequency of each trinucleotide context as it would occur across the whole genome
            signatures.ref = ref,
            signature.cutoff = 0.08
        )$weights
        row.sum <- rowSums(row.df)
        row.df <- row.df / row.sum
        rownames(row.df) <- x
        return(as.data.frame(row.df))
    })

    deconstructSigs.results <- plyr::rbind.fill(sig.contrib.list)
    rownames(deconstructSigs.results) <- rownames(sigs.input)

    return(deconstructSigs.results)
}

main <- function() {
    parser <- ArgumentParser(description = DESCRIPTION)
    parser$add_argument("--maf_path", help = "MAF file path")
    parser$add_argument("--signature_path", help = "signature reference path")
    parser$add_argument("--output_folder", help = "path to store output data")
    args <- parser$parse_args()
    result <- MutSigEstimate(maf_path = args$maf_path, signature_path = args$signature_path)

    # store results
    dir.create(args$output_folder, showWarnings = FALSE)
    write.csv(result, file.path(args$output_folder, "mutation_signature.csv"), quote = F)
}

if (!interactive()) {
    main()
}