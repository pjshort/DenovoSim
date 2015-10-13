# core routines for TF binding analysis

get_sequence <- function(chr, start, stop, version = "hg18") {
  
  # input: (multiple) chr, start, stop, hg version (defaults to hg19)
  # output: list of sequences as DNAStrings object for each input
  
  if (version == "hg19"){
    library(BSgenome.Hsapiens.UCSC.hg19)
  } else if (version == "hg18"){
    library(BSgenome.Hsapiens.UCSC.hg18) # TODO need to add download of hg18 to build.R
  }
  
  if (!all(grepl(pattern = "^chr", chr))){  # assert that chromosome column have chr in front
    warning("Not all entries in the chromosome column start with \"chr\" - try reformatting this column e.g. \"chrX\" instead of \"X\" with paste0(\"chr\",chr_number")
    chr = paste0("chr", chr)
  }
  
  seqs = getSeq(Hsapiens, chr, start, stop)
  return(seqs)
}


get_alt_sequence <- function(sequence, sub_position, alt) {
  
  # input: sequence, positions where alteration has occured, alteration to substitute in
  # output: new alt_sequence
  
  alt_sequence = sequence
  alt_sequence[sub_position] = alt # TODO: add a test to ensure that sub_position and alt are the same length

  return(alt_sequence)
}

### Annotated sequences with predicted TF binding

single_sequence_coverage <- function(seq, rel_pos, pwm_list, min.score = "95%"){
  
  # input: DNAString sequence, list of PWMs to query, min.score (optional)
  # returns: site
  # returns all regions predicted to have TFB affinity >= min.score
  
  # TODO: write test for this section
  
  # scan full list of PWMs against the sequence provided
  site_seq_list = searchSeq(pwm_list, seq, seqname="ref_sequence", min.score=min.score, strand="*")
  
  # keep only the TFs that have a hit greater than min score
  interval_hits = site_seq_list[which(sapply(site_seq_list, length) > 0)]
  
  # keep only the TFs that have a hit in region that overlaps with the de novo
  overlaps_dn = sapply(interval_hits, function(t) t[(rel_pos >= start(t@views@ranges)) & (rel_pos <= end(t@views@ranges))])
  
  # filter list to remove the empty TFs
  pos_hits = overlaps_dn[which(sapply(overlaps_dn, length) > 0)]

  return(pos_hits) # returns a (possibly empty) list of SiteSet objects
}

scan_regions <- function(sequences, rel_positions, pwm_list, min.score = "95%"){
  
  # input: vector of sequences, vector of relative positions of de novo within sequence, list of PWMs to query, minimum binding score (optional)
  # output: list with one element for each pair of seq, rel_pos that contains predicted de novo binding events (if any)
  
  scan_results = mapply(single_sequence_coverage, sequences, rel_positions, MoreArgs = list("pwm_list" = pwm_list))
  
  return(scan_results)
}

LOBGOB_scan <- function(ref_seq, rel_pos, ref, alt, pwm_list, min.score = "95%"){
  
  # input: single ref sequence, single alt sequence, list of PWMs to query, minimum binding score (optional)
  # returns: SiteSetList of original site that was passed plus any sites that were NOT found with ref (but are found with alt)
  # stand for 'loss of binding gain of binding scan'
  
  # TODO: alter so rel_pos can be a range instead of a point!
  
  if (typeof(ref_seq) != "DNAString") { ref_seq = DNAString(ref_seq)}
  
  alt_seq = get_alt_sequence(ref_seq, rel_pos, alt)
  
  # scan against all PWMs with the reference sequence and alt (after mutation)
  # the only differences between scan results should be as due to a change in binding affinity due to the mutation
  ref_results = single_sequence_coverage(ref_seq, rel_pos, pwm_list, min.score = min.score)
  alt_results = single_sequence_coverage(alt_seq, rel_pos, pwm_list, min.score = min.score)
  alt_results = alt_results[!(names(alt_results) %in% names(ref_results))] # only the binding events NOT already spotted in ref
  
  # accounts for one de novo hitting same TF binding motif twice
  n_ref_bindings = sapply(ref_results, length)
  n_alt_bindings = sapply(alt_results, length) # these are actually only the alts that are true GOB (i.e. don't show up in ref)
  
  # scan ref_pwms for change due to mutation (alt) - tag with LOB if score decreases and GOB if score increases
  ref_binding_change = lapply(ref_results, function(r) binding_change(r, rel_pos, ref, alt))
  ref_binding_change = do.call(rbind, ref_binding_change)
  
  # look at score change for alt_results - these must be higher in alt and lower score in ref (below min.score threshold)
  # note, the score output will be transposed! (alt_score, ref_score)
  alt_binding_change = lapply(alt_results, function(r) binding_change(r, rel_pos, alt, ref))
  alt_binding_change = do.call(rbind, alt_binding_change)
  
  all_names = c(rep(names(ref_results), n_ref_bindings), rep(names(alt_results), n_alt_bindings))
  
  if (length(all_names) == 0) { # no binding results for ref or alt
    return( NULL )
  }
  
  if (length(names(ref_results)) > 0 & length(names(alt_results)) > 0) { 
    if (length(all_names) != length(c(ref_binding_change[,1], alt_binding_change[,2]))){
      print(ref_seq)
    }
    binding_changes = data.frame("jaspar_internal" = all_names, 
                               "ref_score" = c(ref_binding_change[,1], alt_binding_change[,2]),
                               "alt_score" = c(ref_binding_change[,2], alt_binding_change[,1]))
  } else if (length(names(ref_results)) > 0) { # only binding for ref sequence
    binding_changes = data.frame("jaspar_internal" = rep(names(ref_results), n_ref_bindings), 
                                 "ref_score" = ref_binding_change[,1],
                                 "alt_score" = ref_binding_change[,2])
  } else { # only binding for alt sequence
    binding_changes = data.frame("jaspar_internal" = rep(names(alt_results), n_alt_bindings), 
                                 "ref_score" = alt_binding_change[,2],
                                 "alt_score" = alt_binding_change[,1])
  }

  binding_changes$diff = binding_changes$ref_score - binding_changes$alt_score
  
  binding_changes$result = ifelse(binding_changes$diff > 0, "LOSS", "GAIN")
  binding_changes$result[binding_changes$diff == 0] = "SILENT"
  
  return(binding_changes)
  
}

split_site_set <- function(ss){
  
  # input: SiteSet
  # output: two or more SiteSets with one row each from original site set
  
  l = length(ss)
  s = list(ss[1])
  if (l > 1){
    for (i in seq(2, l)){
      s = c(s, ss[i])
    }
  }
  return(s)
}

### calculate ref vs. alt change in binding

binding_change <- function(site_set, rel_pos, ref, alt, min.score = "95%"){
  
  # input: SiteSet object (TFBSTools), relative position of de novo, ref, alt
  # returns: ref_score, alt_score in 2x1 vector
  
  # TODO: reformulate to allow indels
    
  # get the position of the de novo within the motif
  motif_pos = rel_pos - start(site_set@views@ranges) + 1
  
  # get position weight matrix
  pwm = site_set@pattern@profileMatrix
  
  # get score with reference allele
  ref_score = site_set@score
  
  # compute change with alt allele and add to ref to calculate new alt allele binding score
  change = as.numeric(pwm[alt, motif_pos] - pwm[ref, motif_pos])
  alt_score = ref_score + change
  
  return(cbind(ref_score, alt_score))
}
