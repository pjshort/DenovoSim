# create de novo files by simulation

relative_seq_probabilities <- function(regions){
  
  # get probability of mutation at each point in sequence (based on sequence context) and normalize
  seq_probabilities = sapply(as.character(regions$seq), function(s) p_position(s, normalize = TRUE), USE.NAMES = FALSE)
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
