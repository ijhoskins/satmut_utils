# Functions useful for summarization of data generated by satmut_utils

library(data.table)
library(e1071)

VCF_SUMMARY_V1.0_COLCLASSES = c("#CHROM"="character", "POS"="integer", "ID"="character", "REF"="character", "ALT"="character", "QUAL"="numeric", "FILTER"="character", 
                                "POS_NT"="integer", "REF_NT"="character", "ALT_NT"="character", "DP"="integer", "CAO"="integer", "NORM_CAO"="numeric", "CAF"="numeric",
                                "R1_PLUS_AO"="integer", "R1_MINUS_AO"="integer", "R2_PLUS_AO"="integer", "R2_MINUS_AO"="integer",
                                "R1_PLUS_MED_RP"="numeric", "R1_MINUS_MED_RP"="numeric", "R2_PLUS_MED_RP"="numeric", "R2_MINUS_MED_RP"="numeric",
                                "R1_PLUS_MED_BQ"="numeric", "R1_MINUS_MED_BQ"="numeric", "R2_PLUS_MED_BQ"="numeric", "R2_MINUS_MED_BQ"="numeric",
                                "R1_PLUS_MED_NM"="numeric", "R1_MINUS_MED_NM"="numeric", "R2_PLUS_MED_NM"="numeric", "R2_MINUS_MED_NM"="numeric",
                                "LOCATION"="character", "REF_CODON"="character", "ALT_CODON"="character", "REF_AA"="character", "ALT_AA"="character", 
                                "AA_CHANGE"="character", "AA_POS"="character", "MATCHES_MUT_SIG"="character")

AA_CODONS_LIST<- list("A"=c("GCC", "GCT", "GCA", "GCG"), "C"=c("TGC", "TGT"), "D"=c("GAC", "GAT"), "E"=c("GAG", "GAA"), "F"=c("TTC", "TTT"), 
                 "G"=c("GGC", "GGG", "GGA", "GGT"), "H"=c("CAC", "CAT"), "I"=c("ATC", "ATT", "ATA"), "K"=c("AAG", "AAA"), 
                 "L"=c("CTG", "CTC", "TTG", "CTT", "CTA", "TTA"), "M"=c("ATG"), "N"=c("AAC", "AAT"), 
                 "P"=c("CCC", "CCT", "CCA", "CCG"), "Q"=c("CAG", "CAA"), "R"=c("CGC", "AGG", "CGG", "AGA", "CGA", "CGT"), 
                 "S"=c("AGC", "TCC", "TCT", "AGT", "TCA", "TCG"), "T"=c("ACC", "ACA", "ACT", "ACG"), "V"=c("GTG", "GTC", "GTT", "GTA"), 
                 "W"=c("TGG"), "Y"=c("TAC", "TAT"), "*"=c("TGA", "TAG", "TAA"))

STOP_CODONS<- c("TGA", "TAG", "TAA")

aa_codon_counts<- sapply(AA_CODONS_LIST, length)

aa_codon_proportions<- aa_codon_counts/sum(aa_codon_counts)

positive_charged_aas<- c("R", "H", "K")
negative_charged_aas<- c("D", "E")
polar_uncharged_aas<- c("S", "T", "N", "Q") # C
nonpolar_uncharged_aas<- c("A", "V", "I", "L", "M", "F", "Y", "W") # G, P
other_aas<- c("C", "G", "P")

aa_map<- list("A"="Ala", "C"="Cys", "D"="Asp", "E"="Glu", "F"="Phe", 
              "G"="Gly", "H"="His", "I"="Ile", "K"="Lys", "L"="Leu", 
              "M"="Met", "N"="Asn", "P"="Pro", "Q"="Gln", "R"="Arg", 
              "S"="Ser", "T"="Thr", "V"="Val", "W"="Trp", "Y"="Tyr",
              "*"="Ter")

aa_class_list<- list("Positive_charged"=positive_charged_aas, "Negative_charged"=negative_charged_aas, 
                     "Polar_uncharged"=polar_uncharged_aas, "Nonpolar_uncharged"=nonpolar_uncharged_aas, 
                     "Other"=other_aas, "Stop"=c("*"))

positive_charged_codons<- c("CGC", "AGG", "CGG", "AGA", "CGA", "CGT",   "CAC", "CAT",   "AAG", "AAA")
negative_charged_codons<- c("GAC", "GAT",   "GAG", "GAA")
polar_uncharged_codons<- c("AGC", "TCC", "TCT", "AGT", "TCA", "TCG",   "ACC", "ACA", "ACT", "ACG",   "AAC", "AAT",   "CAG", "CAA")
nonpolar_uncharged_codons<- c("GCC", "GCT", "GCA", "GCG",   "GTG", "GTC", "GTT", "GTA",   "ATC", "ATT", "ATA",   
                              "CTG", "CTC", "TTG", "CTT", "CTA", "TTA",   "ATG",   "TTC", "TTT",   "TAC", "TAT",   "TGG")


#' Reads in variant calls in vcf.summary.txt format
#'
#' @param indir character path of input directory
#' @param colclasses named character vector specifying class of each column 
#' @return data.table of all samples in the input dir
read_vcf_summary_files<- function(indir=".", colclasses=VCF_SUMMARY_V1.0_COLCLASSES){
  
  vcf_summary_file_pattern<- "var.cand.vcf.summary.txt"
  vcf_summary_files<- dir(path=indir, pattern=vcf_summary_file_pattern, full.names=TRUE)
  
  if(length(vcf_summary_files)==0){
    stop(paste("No var.cand.vcf.summary.txt files found in directory", indir))
  }
         
  vcf_summary_files_split<- sapply(strsplit(as.character(vcf_summary_files), "/", fixed=TRUE), tail, 1)
  sample_vector<- sapply(strsplit(as.character(vcf_summary_files_split), ".", fixed=TRUE), head, 1)
  
  vcf_summary_file_list<- lapply(vcf_summary_files, fread, sep="\t", header=TRUE, colClasses=colclasses)
  vcf_summary_file_list_nrows<- sapply(vcf_summary_file_list, nrow)
  
  # Run with fill=TRUE just in case any files differ in the number of columns
  vcf_summary_file_list_all<- rbindlist(vcf_summary_file_list, fill=TRUE)
  vcf_summary_file_list_all[,Sample:=factor(rep(sample_vector, vcf_summary_file_list_nrows), levels=sample_vector)]
  return(vcf_summary_file_list_all)
}


#' Gets the type of variant given a VCF REF and ALT field
#'
#' @param ref character REF
#' @param alt character ALT
#' @return character one of {"SNP", "di_nt_MNP", "tri_nt_MNP", "MNP", "[123456]_nt_HAPLO", "INS", "DEL"}
get_var_type<- function(ref, alt, split_mnps=TRUE){
  
  len_ref<- nchar(ref)
  len_alt<- nchar(alt)
  ref_c<- unlist(strsplit(ref, ""))
  alt_c<- unlist(strsplit(alt, ""))
  
  if(len_ref == len_alt){
    if(len_ref == 1){
      return("SNP")
      }
    else if(split_mnps) {
      if(len_ref == 2){
        return("di_nt_MNP")
      }
      hd<- hamming.distance(sub(",", "", ref_c, fixed=T), sub(",", "", alt_c, fixed=T))
      if(len_ref == 3){
        if(hd == 2){
          if(grepl(",", ref, fixed=T)){
            return("2_nt_HAPLO")
          } else {
            return("di_nt_MNP")
          }
        } else {
          return("tri_nt_MNP")
        }
      } else {
        if(hd == 2){
          return("2_nt_HAPLO")
        } else if (hd == 3){
          return("3_nt_HAPLO")
        } else if (hd == 4){
          return("4_nt_HAPLO")
        } else if (hd == 5){
          return("5_nt_HAPLO")
        } else if (hd == 6){
          return("6_nt_HAPLO")
        }
      }
    } else{
    return("MNP")
    }
  } else if(len_ref < len_alt){
    return("INS")
  } else if(len_ref > len_alt){
    return("DEL")
  }
}


#' Annotates data.table of vcf.summary.txt formatted data
#'
#' @param in_dt data.table in vcf.summary.txt format
#' @param mapping_dt data.table mapping sample name to desired fraction label
#' @param tile_positions data.table with POS and TILE columns, for annotating PCR tiles
#'   NULL if no tile annotation desired.
#' @return data.table of all samples in the input directory.
annotate_calls<- function(in_dt, mapping_dt=NULL, tile_positions_dt=NULL){
  
  copy_dt<- copy(in_dt)
  
  copy_dt[,VAR_ID:=as.factor(paste(`#CHROM`,POS,REF,ALT,sep=":"))]
  
  if(!is.null(tile_positions_dt)){
    tile_positions_match<- match(copy_dt$POS, tile_positions_dt$POS)
    copy_dt[,TILE:=as.factor(tile_positions_dt[tile_positions_match,TILE])]
  }
  
  copy_dt[,TYPE:=as.factor(get_var_type(REF, ALT)), by=seq_len(nrow(copy_dt))]
  copy_dt[,SUBST:=as.factor(paste(REF_NT, ALT_NT, sep=":"))]
  copy_dt[,REF_NT:=as.factor(REF_NT)]
  copy_dt[,ALT_NT:=as.factor(ALT_NT)]
  copy_dt[,`#CHROM`:=as.factor(`#CHROM`)]
  copy_dt[,FILTER:=as.factor(FILTER)]
  copy_dt[,UP_REF_NT:=as.factor(UP_REF_NT)]
  copy_dt[,DOWN_REF_NT:=as.factor(DOWN_REF_NT)]
  copy_dt[,AA_POS2:=as.numeric(sapply(
    as.character(AA_POS), summarize_str_list)),
          by=seq_len(nrow(copy_dt))]
  
  copy_dt[,VAR_TYPE:=ifelse(REF_AA==ALT_AA, "Silent", "Missense"), 
          by=seq_len(nrow(copy_dt))]
  
  copy_dt[ALT_AA=="*", VAR_TYPE:="Nonsense"]
  copy_dt[,VAR_TYPE:=factor(VAR_TYPE, levels=c("Silent", "Missense", "Nonsense"))]
  
  if(!is.null(mapping_dt)){
    copy_dt<- add_fraction_ids(copy_dt, mapping_dt)
  }
  
  copy_dt[,log10_CAF:=log10(CAF)]
  
  return(copy_dt)
}


#' Gets the median and returns as numeric
#'
#' @param x integer or numeric vector
#' @return numeric value of median
numeric_median<- function(x){
  res<- as.numeric(median(x))
  return(res)
}


#' Collapses per-bp records into per-variant records
#'
#' @param in_dt data.table in vcf.summary.txt format
#' @param key_features character vector of keys to merge on. (Include your other variables).
#' @param agg_fun function to use to aggregate for each level.
#' @return data.table with per-bp statistics collapsed/aggregated
#' @details needed to collapse multiple records for MNPs into a single record
collapse_mnps<- function(in_dt, key_features=c("Sample", "VAR_ID"), agg_fun=numeric_median){
  
  copy_dt<- copy(in_dt)
  setkeyv(copy_dt, key_features)
  
  # The following is just for reference
  per_bp_features<- c("POS_NT", "REF_NT", "ALT_NT", "UP_REF_NT", "DOWN_REF_NT")
  
  numeric_col_classes<- sapply(copy_dt, class)
  which_numeric<- numeric_col_classes%in%c("integer", "numeric", "single", "double")
  
  numeric_cols_dt<- copy_dt[,which_numeric,with=FALSE]
  nonnumeric_cols_dt<- copy_dt[,!which_numeric,with=FALSE]
  
  # Add back in the key columns
  for(i in 1:length(key_features)){
    numeric_cols_dt[,key_features[i]:=nonnumeric_cols_dt[,key_features[i],with=FALSE]]
  }
  #numeric_cols_dt[,Sample:=nonnumeric_cols_dt[,Sample]]
  #numeric_cols_dt[,VAR_ID:=nonnumeric_cols_dt[,VAR_ID]]
  
  numeric_features<- names(numeric_cols_dt)
  # numeric_features is evaluated lazily in .SDcols, so explicitly remove the keys
  numeric_features<- numeric_features[
    !numeric_features%in%key_features]
  
  nonnumeric_features<- names(nonnumeric_cols_dt)
  nonnumeric_features<- nonnumeric_features[
    !nonnumeric_features%in%c(key_features)]
  
  # Aggregate numeric features; most features have duplicated values, 
  # but BQs are per-bp features to be aggregated
  numeric_cols_agg_dt<- numeric_cols_dt[
    ,lapply(.SD, agg_fun), 
    by=key_features, 
    .SDcols=numeric_features]
  
  # Similarly, we must collapse the nonnumeric features
  nonnumeric_cols_collapse_dt<- nonnumeric_cols_dt[
    ,lapply(.SD, head, 1), 
    by=key_features, 
    .SDcols=nonnumeric_features]
  
  merged_dt<- merge(nonnumeric_cols_collapse_dt, 
                    numeric_cols_agg_dt, 
                    by=key_features)
  
  return(merged_dt)
}


#' Applies a simple AF subtraction approach for removal of error
#'
#' @param in_dt data.table in vcf.summary.txt format
#' @return data.table with AF subtract of the negative applied to other samples
subtract_af<- function(in_dt){
  
  copy_dt<- copy(in_dt)
  neg_dt<- copy_dt[Input=="Negative_control", .(VAR_ID, CAO, CAF, NORM_CAO)]
  
  # Merge the two tables so we can easily subtract the negative frequencies via a vectorized operation
  merged_dt<- merge(copy_dt[Input!="Negative_control"], neg_dt, by="VAR_ID", all.x=TRUE)
  
  # Set baseline to 0 for variants unique to a polysomal fraction
  merged_dt[is.na(CAF.y), CAF.y:=0.0]
  merged_dt[is.na(CAO.y), CAO.y:=0]
  merged_dt[is.na(NORM_CAO.y), NORM_CAO.y:=0.0]

  # Finally subtract the frequencies and remove those variants with negative resultant frequencies
  merged_dt[,CAF:=CAF.x-CAF.y]
  merged_dt[,CAO:=CAO.x-CAO.y]
  merged_dt[,NORM_CAO:=NORM_CAO.x-NORM_CAO.y]
  filtered_dt<- merged_dt[CAF>0,]
  filtered_dt[,log10_CAF:=log10(CAF)]
  
  return(filtered_dt)
}


#' Reverse complements a sequence
#'
#' @param nt_seq nucleotide sequence, should contain only {A,T,C,G,N}, case insensitive
#' @return character reverse complemented sequence
revcomp<- function(nt_seq){
  
  if(grepl("[^ATCGUNatcgun]", nt_seq)){
    stop("Sequence has invalid characters")
  }
  
  seq_rev<- paste(rev(unlist(strsplit(nt_seq, ""))), collapse="")
  seq_rev_sub<- chartr("ATCGUNatcgun", "TAGCANtagcan", seq_rev)
  return(seq_rev_sub)
}


#' Summarizes values in a delimited list
#'
#' @param str_list character of delimited integers to summarize
#' @param split_sep character delimiter of values
#' @param summarization_stat function for summarization
#' @return numeric summarized values
summarize_str_list<- function(str_list, split_sep=",", summarization_stat="median"){
  split_unlisted<- as.integer(unlist(strsplit(str_list, split_sep)))
  summarized<- get(summarization_stat)(split_unlisted)
  return(summarized)
}


#' Annotates and collapses per-bp records into per-variant records
#'
#' @param in_dt data.table in vcf.summary.txt format
#' @param mapping_dt data.table mapping sample name to desired fraction label
#' @param tile_positions data.table with POS and TILE columns, for annotating PCR tiles. 
#' NULL if no tile annotation desired.
#' @param key_features character vector of keys to merge on. (Include your other variables).
#' @return data.table
postprocess_vcf_summary_dt<- function(in_dt, mapping_dt=NULL, tile_positions_dt=NULL, 
                                      key_features=c("Sample", "VAR_ID")){
  
  copy_filt_dt<- in_dt[LOCATION=="CDS"]
  annot_dt<- annotate_calls(copy_filt_dt, mapping_dt, tile_positions_dt)
  collapse_dt<- collapse_mnps(in_dt=annot_dt, key_features=key_features)
  return(collapse_dt)
}


#' Computes cumulative coverage depth over positions
#'
#' @param in_dt data.table in vcf.summary.txt format
#' @param key_features character vector of keys to merge on. (Include your other variables).
#' @param integer stepsize: compute proportion between 0 and max coverage with this stepsize.
#' a smaller stepsize returns a more granular profile but takes longer time to compute.
#' @return data.table
waterfall_ecdf<- function(in_dt, key_features=c("Sample"), stepsize=10){
  
  max_cov<- in_dt[,max(V4)]
  
  # Modify max_cov slighlty so we can have an appropriate stepsize
  # to split data evenly
  while(max_cov%%stepsize == 0){
    max_cov<- max_cov + 1
  }
  
  # We compute the ECDF then plot values for 0 to max coverage
  ecdf_res_dt<- in_dt[,.(Proportion_greater_than=1-ecdf(V4)(seq(0,max_cov,stepsize)),
                         Coverage=seq(0,max_cov,stepsize)), by=key_features]
  return(ecdf_res_dt)
}


#' Re-computes CAF by averaging CAF over positions in MNPs
#'
#' @param in_dt data.table vcf.summary.txt per-bp calls prior to collapse
#' @return data.table with CAF updated for MNPs
recompute_caf<- function(in_dt){
  copy_dt<- copy(in_dt)
  recomp_dt<- copy_dt[,CAF:=max(CAO/DP), by=.(VAR_ID)]
  return(recomp_dt)
}


#' Converts long to short AA formats
#'
#' @param aa_change character HGVS amino acid change, e.g. p.Ala41Thr
#' @return character compact (short) amino acid change annotation, e.g. p.A41T
convert_aa_notation<- function(aa_change){
  
  aa_change_strip<- strsplit(as.character(aa_change), ".", fixed=T)
  
  if(length(aa_change_strip[[1]])==1){
    aa_change_strip<- aa_change
  } else{
    aa_change_strip<- sapply(aa_change_strip, "[", 2)
  }
  
  pos<- gsub(pattern="[[:alpha:]]|[[:punct:]]", replacement="", aa_change_strip)
  
  aa_change_no_int<- gsub(pattern="[[:digit:]]", replacement=":", aa_change_strip)
  
  aa_change_no_int_split<- unlist(strsplit(aa_change_no_int, ":+"))
  
  final_res<- c()
  for(e in 1:length(aa_change_no_int_split)){
    
    if(e==2){
      final_res<- append(final_res, pos)
      
      if(grepl("[[:punct:]]", aa_change_no_int_split[e])){
        match_res<- names(aa_map)[which(aa_change_no_int_split[1]==aa_map)]
        final_res<- append(final_res, match_res)
        next
      }
    }
    match_res<- names(aa_map)[which(aa_change_no_int_split[e]==aa_map)]
    final_res<- append(final_res, match_res)
  }
  
  final_res<- paste0(c("p.",final_res), collapse="")
  return(final_res)
}


#' Gets the type of amino acid
#'
#' @param x character single-letter amino acid code
#' @return character one of {"Positive", "Negative", "Polar", "Nonpolar", "Other"}
get_aa_change_type<- function(x){
  
  if(x%in%positive_charged_aas){
    return("Positive")
  } else if(x%in%negative_charged_aas){
    return("Negative")
  } else if(x%in%polar_uncharged_aas){
    return("Polar")
  }else if(x%in%nonpolar_uncharged_aas){
    return("Nonpolar")
  } else{
    return("Other")
  }
}


#' lapply-friendly version of %in%
#'
#' @param e_set character vector containing values to search in
#' @param x character to search
#' @return logical whether or not the x was found in e_set
is_in<- function(e_set, x){
  return(x%in%e_set)
}


#' Reads individual count tables generated by run_satmut_utils_to_dimsum.py and merges them 
#'
#' @param indir character path of input directory containing DiMSum_counts.txt tables
#' @return data.table merged tables
merge_dimsum_tables<- function(indir="."){
  
  file_pattern<- "DiMSum_counts.txt"
  dimsum_files<- dir(path=indir, pattern=file_pattern, full.names=TRUE)
  
  if(length(dimsum_files)==0){
    stop(paste("No DiMSum_counts.txt files found in directory", indir))
  }
  
  dimsum_files_split<- sapply(strsplit(as.character(dimsum_files), "/", fixed=TRUE), tail, 1)
  sample_vector<- sapply(strsplit(as.character(dimsum_files_split), ".", fixed=TRUE), head, 1)
  
  dimsum_files_list<- lapply(dimsum_files, fread, sep="\t", header=TRUE)
  dimsum_files_list_nrows<- sapply(dimsum_files_list, nrow)
  
  dimsum_files_list_all_dt<- rbindlist(dimsum_files_list)
  dimsum_files_list_all_dt[,Sample:=factor(rep(sample_vector, dimsum_files_list_nrows), levels=sample_vector)]
  
  # Get mapping between VAR_ID and nt_seq
  var_id_to_nt_string_dt<- unique(dimsum_files_list_all_dt[,.(VAR_ID, nt_seq)])
  
  # Cast counts using VAR_ID; then add nt string annotations back
  dimsum_files_list_cast_dt<- dcast(
    data=dimsum_files_list_all_dt[,-c("nt_seq")],
    formula=VAR_ID~Sample, value.var="counts", fill=0)
  
  nt_seqs<- var_id_to_nt_string_dt[
    match(dimsum_files_list_cast_dt[,VAR_ID], VAR_ID), nt_seq]
  
  res_dt<- cbind(data.table(nt_seq=nt_seqs), dimsum_files_list_cast_dt)
  
  return(res_dt)
}


numeric_to_integer<- function(x){
  res<- as.integer(round(x))
  return(res)
}

#' Correct test sample variant counts by subtraction of negative control background 
#'
#' @param dt data.table merged DiMSum_counts.txt tables generated from merge_dimsum_tables()
#' @param test_samples list test samples, paired one-to-one with NC sample(s) in nc_samples 
#' @param nc_samples list negative control samples, paired one-to-one with elements in test_samples
#' @return data.table of rounded, corrected CPMs for test samples only (no NC samples)
subtract_neg_control_cpms<- function(dt, test_samples, nc_samples){
  
  if(length(test_samples)!=length(nc_samples)){
    stop("Length mismatch between test_samples and nc_samples lists.")
  }
  
  copy_dt<- copy(dt)
  noncount_cols<- which(sapply(copy_dt, class)!="integer")
  count_cols<- which(sapply(copy_dt, class)=="integer")
  
  copy_noncounts_dt<- copy_dt[,noncount_cols,with=FALSE]
  copy_counts_dt<- copy_dt[,count_cols,with=FALSE]
  
  # Generate CPM dataset prior to background subtraction
  # Use total counts to estimate library size normalization factors
  colsums<- colSums(copy_counts_dt)
  colsums_cpm_factor<- 1e6/colsums
  copy_cpm_dt<- as.data.table(t(t(copy_counts_dt)*colsums_cpm_factor))
  
  # Determine which variant/sequence is the WT count
  wt_var_idx<- which.max(rowSums(copy_counts_dt))
  copy_cpm_wt_dt<- copy_cpm_dt[wt_var_idx,]
  copy_cpm_nonwt_dt<- copy_cpm_dt[!wt_var_idx,]
  
  # Subtract negative control CPM from samples, and add this CPM
  # to the wildtype sequence (assumed max count among variants/sequences)
  # there may be multiple pairings between test and NC samples
  copy_cpm_nonwt_subtract_dt<- NULL
  copy_cpm_wt_add_dt<- NULL
  for(i in 1:length(nc_samples)){
    
    # Subtract NC from test samples
    sample_cpm_dt<- copy_cpm_nonwt_dt[,test_samples[[i]], with=FALSE]
    nc_cpms<- unlist(copy_cpm_nonwt_dt[,nc_samples[[i]], with=FALSE])
    sample_cpm_subtract_dt<- sample_cpm_dt-nc_cpms
    copy_cpm_nonwt_subtract_dt<- cbind(copy_cpm_nonwt_subtract_dt, sample_cpm_subtract_dt)
    
    # Add back sum of subtracted CPMs to the WT
    wt_cpm_dt<- copy_cpm_wt_dt[,test_samples[[i]], with=FALSE]
    wt_cpm_add_dt<- wt_cpm_dt + sum(nc_cpms)
    copy_cpm_wt_add_dt<- cbind(copy_cpm_wt_add_dt, wt_cpm_add_dt)
  }
  
  # Round the numeric values and convert to integer
  copy_cpm_nonwt_subtract_int_dt<- as.data.table(apply(
    copy_cpm_nonwt_subtract_dt, 2, numeric_to_integer))
  
  copy_cpm_wt_add_int_dt<- as.data.table(lapply(
    copy_cpm_wt_add_dt, numeric_to_integer))
  
  # Now add back in the VAR_ID and nt_seq labels
  copy_cpm_nonwt_subtract_dt<- cbind(
    copy_noncounts_dt[!wt_var_idx,], copy_cpm_nonwt_subtract_int_dt)
  
  copy_cpm_wt_add_dt<- cbind(
    copy_noncounts_dt[wt_var_idx,], copy_cpm_wt_add_int_dt)
  
  copy_cpm_corrected_dt<- rbind(copy_cpm_wt_add_dt, copy_cpm_nonwt_subtract_dt)
  
  return(copy_cpm_corrected_dt)
}


