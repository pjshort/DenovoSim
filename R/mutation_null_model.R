# Null model for de novo regulatory mutations

# Probability of mutation in middle base of each trinucleotide is calculated from empirical data.
# The probabilities are assumed poisson, independent and summed across each sequence

mut_rates <- read.table("../data/forSanger_1KG_mutation_rate_table.txt", header=TRUE)

### SNP MODEL - based on sequence context ####

# the indel null model is inferred directly from global snp mutation rate and data - see rupit_core.R for generate_indel_null function

p_base <- function(from){
  
  # probability of mutation from base in position 2 to any of three other possible bases
  
  p = mut_rates$mu_snp[c(mut_rates$from == from)]
  return(sum(p))
}

p_position <- function(sequence, normalize = FALSE){
  
  # return vector length nchar(sequence) - 2 with probability of mutation at each base
  
  sequence = as.character(sequence)
  p = sapply(seq(1:nchar(sequence)-2), function(z) p_base(substr(sequence, z, z+2)))
  if (normalize == FALSE){
    return(p)
  } else {
    return(p/sum(p))
  }
}

p_sequence <- function(sequence){
  
  # sum p_all across each trinucleotide sliding along sequence
  
  p = p_position(sequence)
  return(sum(as.numeric(p)))
  
}

sample_alt <- function(ref_tri) {
  
  # given a ref sequence and mutated position, sample a new alt based on null model (middle base changes)
  
  m = mut_rates[which(mut_rates$from == ref_tri), c("to", "mu_snp")]
  tri_alt = sample(m$to, size=1, prob=m$mu_snp)
  return(as.character(substr(tri_alt,2,2)))
}

