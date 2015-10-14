# create de novo files by simulation

relative_haploid_seq_probabilities <- function(regions){
  
  # get probability of mutation at each point in sequence (based on sequence context) and normalize
  seq_probabilities = sapply(as.character(regions$seq), function(s) p_position(s, normalize = TRUE), USE.NAMES = FALSE)
  names(seq_probabilities) = regions$region_id
  return(seq_probabilities)
}

absolute_haploid_seq_probabilities <- function(regions){
  
  # get probability of mutation at each point in sequence (based on sequence context) - does not normalize
  seq_probabilities = sapply(as.character(regions$seq), function(s) p_position(s, normalize = FALSE), USE.NAMES = FALSE)
  names(seq_probabilities) = regions$region_id
  return(seq_probabilities)
  
}

simulate_de_novos <- function(regions, seq_probabilities, n_snps, n_probands, iteration) {
  
  # takes vector of region_ids to consider, list matching region_id to sequence relative probability, and total snps to simulate
  # samples n_snps region ids (with replacement). for each region, specific relative position will be randomly sampled
  
  region_ids = sample(regions$region_id, n_snps, replace = TRUE, prob = regions$p_relative)
  
  if (n_probands < n_snps) { # TODO: think if this is the best way to sample... seems reasonable to automatically assign one de novo per prband because if they are observed,
    # then they must have at least one de novo... sampling naively creates a distribution that does not match the observed one where n_probands slightly less than n_snps
    proband_ids1 = paste0("sim_person_id_", sample(seq(n_probands), n_probands, replace = FALSE)) # assign one snp to each proband
    proband_ids2 = paste0("sim_person_id_", sample(seq(n_probands), n_snps - n_probands, replace = TRUE)) # assign remaining snps
    proband_ids = c(proband_ids1, proband_ids2)
  } else {
    proband_ids = paste0("sim_person_id_", sample(seq(n_probands), n_snps, replace = TRUE)) # assign remaining snps
  }
  
  rel_pos = sapply(as.character(region_ids), function(id) sample(seq(1,length(seq_probabilities[id][[1]])), 1, prob = seq_probabilities[id][[1]]))
  coords = do.call(rbind, strsplit(region_ids, "\\.")) # chr, start, stop
  pos = as.integer(coords[,2]) + as.integer(rel_pos)
  
  # once the position within the sequence has been chosen, want to choose the mutated (alt) base based on null model
  ref_tri = as.character(mapply(function(s, p) substr(s, p-1, p+1), regions[match(region_ids, regions$region_id), "seq"], rel_pos))
  
  # replace any di nucleotides if simulation fell on edge TODO write test to check this works
  di_idx = as.integer(which(sapply(ref_tri, nchar) < 3))
  if (length(di_idx) > 0) {
    ref_tri[di_idx] = as.character(get_sequence(chr = paste0("chr", as.character(coords[di_idx, 1])), start = pos[di_idx] - 1, stop = pos[di_idx] + 1))
  }
  
  ref = as.character(sapply(ref_tri, function(s) substr(s,2,2)))
  alt = as.character(sapply(ref_tri, sample_alt))
  
  sim_dn = data.frame("person_stable_id" = proband_ids, "chr" = coords[,1], "pos" = pos, "ref" = ref, "alt" = alt, "iteration" = iteration)
  
  return(sim_dn)  # data frame
}

simulate_proband <- function(null_probs, n_bases){
  # generate a vector of random numbers (one random number per genomic position)
  p_vector = runif(n_bases)
  mutated_positions = which(p_vector < null_probs)
  return(mutated_positions)
}


record_snp <- function(mutated_position, proband_id, regions, region_break_points, seq_probabilities){
  
  # take a proband id and position and return one line of vcf
  region_idx = sum(region_break_points < mutated_position) + 1
  rel_pos = as.numeric(mutated_position - region_break_points[region_idx - 1])
  coords = strsplit(names(seq_probabilities)[region_idx], "\\.")[[1]]
  chr = coords[1]
  start = as.numeric(coords[2])
  end = as.numeric(coords[3])
  pos = start + rel_pos
  
  print(paste(chr,start,end,rel_pos,sep = "."))
  ref_tri = substr(regions$seq[region_idx], rel_pos-1, rel_pos+1)
  
  if (nchar(ref_tri) < 3) { # simulation de novo fell on end of sequence
    print("LESS THAN 3")
    ref_tri = as.character(get_sequence(chr = paste0("chr", chr), start = pos - 1, stop = pos + 1)) # from S4 to character
  }
  
  ref = substr(ref_tri,2,2)
  alt = sample_alt(ref_tri)
  
  sim_dn = data.frame("person_stable_id" = proband_id, "chr" = chr, "pos" = pos, "ref" = ref, "alt" = alt)
  return(sim_dn)
}

unsupervised_sim <- function(regions, seq_probabilities, n_probands, iteration){
  # generate null probability of mutation at each base (probability is haploid, so multiply by 2)
  null_probs = as.numeric(unlist(seq_probabilities))*2 # note seq_probabilities should be generated with normalize = FALSE
  n_bases = length(null_probs)
  
  region_break_points = cumsum(sapply(seq_probabilities, function(s) length(s)))

  mutation_coords = replicate(n_probands, simulate_proband(null_probs, n_bases))
  muts_per_proband = sapply(mutation_coords, function(m) length(m))
  mutation_coords = unlist(mutation_coords)
  
  proband_ids = paste0("proband.", seq(n_probands), ".iteration.", iteration)
  mutated_proband_ids = rep(proband_ids, muts_per_proband)
  
  sim_dn = do.call(rbind, mapply(record_snp, mutation_coords, mutated_proband_ids, MoreArgs = list("regions" = regions, "region_break_points" = region_break_points, 
                                                 "seq_probabilities" = seq_probabilities), SIMPLIFY = FALSE))
  sim_dn$iteration = iteration
  
  return(sim_dn)
}
