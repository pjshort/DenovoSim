# generate simulation data based on trinucleotide rate
# this can be fed to denovoTF to generate simulation data to compare against experimental data

# INPUT:
# --n_snps -> number of snps to simulate
# --regions -> regions in which SNPs should be simulated
# tab delimited df with chr, start, stop
# --n_probands -> number of patients to simulate

# OUTPUT:
# de novo output file with columns "unique_id", "chr", "pos", "ref", "alt", "tf_name", "jaspar_internal", "ref_score", "alt_score" 
# with ONE ROW PER TF binding event. the output file will likely have more rows than the input file (many more if score threshold is low)

# documentation notes:
# a triple hash (### description xyz) denotes a 'section header' in the code while (# comments..) denotes a more simple comment

### dependencies
library(optparse)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TFBSTools)
library(JASPAR2014)

source("../R/core.R")
source("../R/simulate.R")
source("../R/mutation_null_model.R")

### command line options
option_list <- list(
  make_option("--n_snps", default=0,
              help="Pass the genomic regions that should be annotated with predicted TF binding sites."),
  make_option("--n_probands", default=421,
              help="Number of probands which should be simulated (this can be used to simulate diagnosed/undiagnosed effects)."),
  make_option("--iterations", default=10, help="Set the number of simulation outputs to generate."),
  make_option("--iteration_start", default=1, help="Set where to start the iteration counting from - helpful if running in batches."),
  make_option("--n_chunks", default=2, help = "Number of smaller files to split simulated data into (to reduce memory overhead 
  and allow parallel processing)"),
  make_option("--base_name", default="../data/simulated_dn", help = "Directory to save the chunks. Passing /path/to/chunk 
  will yield /path/to/chunk.1.txt, path/to/chunk.2.txt, etc."),
  make_option("--regions", default="../data/DDD_well_cov_regions.txt",
              help="Pass set of genomic regions for simulation."),
  make_option("--verbose", action="store_true", default=FALSE,
              help="Print extra output advising the user of progression through the code.")
)

args <- parse_args(OptionParser(option_list=option_list))

# load in regions file with required columns: chr, start, stop
regions <- read.table(args$regions, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
regions$region_id <- paste(regions$chr, regions$start, regions$stop, sep = ".")
regions$seq = as.character(get_sequence(regions$chr, regions$start, regions$stop))

if ( args$verbose ) { write("Computing per-base mutation probability for all of the genomic regions that were passed. This may take a little while...", stderr()) }

# get probability of mutation per region and probability of selecting each region (prop of total probability)
if (args$n_snps != 0) {
  # only needed for supervised simulation
  regions$p_snp_null <- sapply(regions$seq, p_sequence)
  regions$p_relative <- regions$p_snp_null/sum(regions$p_snp_null)
}

if ( args$verbose ) { write("Simulating de novos drawn from sequence-context-specific null distribution over regions provided...", stderr()) }

# create large data frame with columns for id, chr, pos, ref, alt, iteration
if (args$n_snps != 0){
  # supervised - condition on number of probands and number of SNPs
  
  seq_probabilities_normalized = relative_haploid_seq_probabilities(regions) # run these two lines to regenerate
  #save(seq_probabilities, file = "../data/sequence_probabilities.out")
  #attach("../data/sequence_probabilities.out")
  
  sim_out = lapply(seq(args$iteration_start, args$iteration_start + args$iterations), function(i) simulate_de_novos(regions, seq_probabilities_normalized, args$n_snps, args$n_probands, i))
  sim_df = do.call(rbind, sim_out)
  sim_df = sim_df[,c("person_stable_id", "chr", "pos", "ref", "alt", "iteration")]
} else { 
  if ( args$verbose ) { write("Simulating de novos without conditioning on the number of SNPs.", stderr()) }
  
  # do not condition on the number of SNPs
  
  seq_probabilities_absolute = absolute_haploid_seq_probabilities(regions) # run these two lines to regenerate
  #save(seq_probabilities, file = "../data/sequence_probabilities.out")
  #attach("../data/sequence_probabilities.out")
  
  sim_out = lapply(seq(args$iteration_start, args$iteration_start + args$iterations), function(i) unsupervised_sim(regions, seq_probabilities_absolute, args$n_probands, i))
  sim_df = do.call(rbind, sim_out)
}

sim_df$iteration = sim_df$iteration + args$iteration_start - 1 # shift iterations up depending on where they are supposed to start
bkp = seq(args$iteration_start - 1, args$iterations, length.out = args$n_chunks + 1)

if ( args$verbose ) { write("Saving simulation files in chunks - feed these to denovoTF to annotate with TF binding predictions.", stderr()) }

if (args$n_chunks != 1){
  for (i in seq(1, args$n_chunks)){
    fname = sprintf("%s.%i.txt", args$base_name, i)
    write.table(sim_df[sim_df$iteration > bkp[i] & sim_df$iteration <= bkp[i+1],], file = fname, col.names = TRUE, row.names = FALSE, sep="\t", quote = FALSE)
    prev_chunk = i
  }
} else {
  write.table(sim_df, file = args$base_name, col.names = TRUE, row.names = FALSE, sep="\t", quote = FALSE)
}
