### From Princy: /work-zfs/abattle4/parsana/gtex_gender/variant_level_analysis/src/gene_level_v2.R

require(data.table)
require(tidyverse)
require(magrittr)
require(foreach)
require(rstan)
require(doMC)
require(eagle1S)
require(doParallel)
registerDoParallel()
getDoParWorkers()


#inputargs = commandArgs(TRUE)
inputargs = c("8", "LIVER.txt")
print(inputargs)
which_chr = paste("chr", inputargs[1], sep = "")
tissue = gsub("\\.txt", "", inputargs[2])
cat(tissue, ":", which_chr, "\n")
topresDir = "/work-zfs/abattle4/parsana/gtex_gender/variant_level_analysis/results/gene_level_with_env/"
resDir = paste(topresDir, tissue, "/", sep = "")
aseDir = "/work-zfs/abattle4/parsana/gtex_gender/variant_level_analysis/data/ase/"
ciseqtlDir = "/work-zfs/abattle4/parsana/gtex_gender/variant_level_analysis/data/ciseqtls/from_meri/"
checkpoint_dir = paste0(resDir, which_chr, "/")
#ase_fn = paste(aseDir, inputargs[2], sep = "")
ase_fn = "eagle_ase.txt"
meta_fn = "/work-zfs/abattle4/parsana/gtex_gender/variant_level_analysis/data/subject_sex.txt"
phased_fn = paste(
  "/work-zfs/abattle4/parsana/gtex_gender/variant_level_analysis/data/phased_genotypes/",
  which_chr,
  ".txt.gz",
  sep = ""
)
gtex_abbrv = "/work-zfs/abattle4/parsana/gtex_gender/variant_level_analysis/data/etc/gtex_colors.txt"
ciseqtl_tiss = read_tsv(gtex_abbrv) %>% filter(tissue_abbrv == tissue) %>% .$tissue_id
ciseqtl_fn = paste(ciseqtlDir, "independent_eqtls_txt.txt", sep = "")

bb_stan = eagle1S:::stanmodels$bb
min_reads = 1000
min_samples = 10

ase = read_tsv(ase_fn) %>% filter(CHR == which_chr)
ase = ase %>% rename(
  CHROM = CHR,
  individual = SUBJECT_ID,
  REF = REF_ALLELE,
  ALT = ALT_ALLELE,
  r = REF_COUNT,
  a = ALT_COUNT
)
ase = ase %>% filter(r + a >= 10)
ase = ase %>% filter(pmin(r,a) >= 3)

individuals_tissue_chr = unique(ase$individual)


meta = read_tsv(meta_fn) %>% rename(individual = SUBJID, x = SEX)
meta = meta[meta$individual %in% individuals_tissue_chr, ]
if (min(table(meta$x)) < 10 | length(table(meta$x)) == 1) {
  print("very few samples")
  quit()
} else{
  dir.create(resDir, showWarnings = F)
  dir.create(checkpoint_dir,
             recursive = T,
             showWarnings = F)
}

phased = read_tsv(phased_fn, comment = "##") %>% rename(CHROM =
                                                          `#CHROM`)

ciseqtl = read.delim(ciseqtl_fn, stringsAsFactors = F)
ciseqtl = ciseqtl %>%
  filter(tissue == ciseqtl_tiss) %>%
  filter(str_detect(variant_id, paste(which_chr, "_", sep = "")))
# filter(abs(tss_distance) <= 1e6)

ciseqtl$gene_id = sapply(strsplit(ciseqtl$gene_id, "\\."), function(x)
  x[1])

# ciseqtl = ciseqtl %>% semi_join(ase, by = c("gene_id" = "GENE_ID"))

### select genes with enough samples
ase_genes = ase %>% group_by(GENE_ID) %>% summarise(n = n()) %>%
  filter(n >= min_samples) %>% drop_na %>% .$GENE_ID
### select genes with enough reads
ase_genes =  ase %>% filter(GENE_ID %in% ase_genes) %>%
  group_by(GENE_ID) %>% summarise(total_count_snp = sum(TOTAL_COUNT)) %>%
  filter(total_count_snp >= 1000) %>% .$GENE_ID

ase_cieqtl_genes = intersect(ase_genes, ciseqtl$gene_id)
ciseqtl = ciseqtl %>% filter(gene_id %in% ase_cieqtl_genes)
ase = ase %>% filter(GENE_ID %in% ase_cieqtl_genes)

exonic_snps = unique(ase$VARIANT_ID)
select_snps_phased = data.frame(
  variant_id = c(ciseqtl$variant_id, exonic_snps),
  stringsAsFactors = F
)

ase = ase %>% select(CHROM, POS, REF, ALT, individual, r, a, GENE_ID)
phased = semi_join(phased, select_snps_phased, by = c("ID" = "variant_id"))

phased_types = c("0|0", "0|1", "1|0", "1|1")
phased_hets = c("0|1", "1|0")

phased = phased %>% unite(pos_alt, POS, ALT, sep = "_", remove =
                            F) %>%
distinct(pos_alt, .keep_all = TRUE)
class(phased) = "data.frame"
rownames(phased) = phased$pos_alt

ase %<>%  mutate(pos_alt = paste(POS, ALT, sep = "_"),
                 individual = as.character(individual)) %>%
  filter(pos_alt %in% rownames(phased)) %>%
  mutate(geno = phased[cbind(as.character(pos_alt), individual)])
unique_genes = unique(ase$GENE_ID)

logit = function(p) {
  log(p / (1 - p))
}
inv_logit = function(g) {
  1 / (1 + exp(-g))
}
min_samples = 10
concShape = 1.001
concRate = 0.001
min_reads = 1000
max_iter = 30
allres = foreach(
  exonic_snp_pos = unique_genes,
  .errorhandling = "remove",
  .combine = bind_rows
) %dopar% {
  # checkpointing per "gene"
  if (!is.null(checkpoint_dir)) {
    check_fn = paste0(checkpoint_dir, exonic_snp_pos, ".txt.gz")
    if (file.exists(check_fn) & !interactive()) {
      return(read_tsv(check_fn) %>%
               mutate(exonic_snp_pos = exonic_snp_pos))
    }
  }
  
  gene_ase = ase %>%
    filter(GENE_ID == exonic_snp_pos) %>%
    select(CHROM, POS, REF, ALT, r, a, individual, geno)
  
  if (nrow(gene_ase) < min_samples) {
    cat(exonic_snp_pos, "sample_size not enough\n")
    return(NULL)
  }
  
  allelic_count_total = sum(gene_ase$a + gene_ase$r)
  cat("Allelic total count ", allelic_count_total, "\n")
  if (allelic_count_total < min_reads) {
    cat(exonic_snp_pos, "readcounts not enough\n")
    return(NULL)
  }
  cis_snps = ciseqtl %>% filter(gene_id == exonic_snp_pos) %>% .$variant_id
  cat("Number of ciseQTLs:", exonic_snp_pos, length(cis_snps), "\n")
  if (length(cis_snps) == 0) {
    cat("no snps for gene", exonic_snp_pos, "\n")
    return(NULL)
  }
  
  # iterate over cis SNPs
  temp_results = foreach(
    snp_pos = cis_snps,
    .errorhandling = if (interactive())
      "stop"
    else
      "remove",
    .combine = bind_rows
  ) %do% {
    cat(exonic_snp_pos, snp_pos, "\n")
    
    reg_geno = (phased %>% filter(ID == snp_pos))[, 11:ncol(phased)] %>% as.matrix()
    
    if (nrow(reg_geno) != 1) {
      print("Skipping >biallelic site")
      return(NULL)
    }
    
    reg_geno = data_frame(individual = colnames(reg_geno),
                          reg_geno = as.character(reg_geno))
    
    # join the ASE and cisSNP phased genotypes
    ase_temp = gene_ase %>% inner_join(reg_geno, by = "individual") %>%
      filter(geno %in% phased_hets, # signal only comes from het exonic SNPs
             reg_geno %in% phased_types, # require phased cisSNP
             (r + a) > 0) %>% # must have some coverage
      mutate(het_x = ifelse(reg_geno %in% phased_hets, # if cisSNP is het
                            ifelse(geno == reg_geno, 1, -1), # is it in phase with the exonicSNP?
                            0))
    
     if (min(table(ase_temp$x)) < min_samples) {
      cat(exonic_snp_pos, snp_pos, "sample_size not enough\n")
      return(NULL)
    }

    if (nrow(ase_temp) < min_samples) {
      cat(exonic_snp_pos, snp_pos, "sample_size not enough\n")
      return(NULL)
    }
    num_het_snps = sum(ase_temp$het_x != 0)
    if (num_het_snps < min_samples) {
      cat(exonic_snp_pos, snp_pos, "sample_size not enough\n")
      return(NULL) # no heterozygous regulatory SNPs
    }
    
    ase_temp %<>% left_join(meta, by = "individual")
    ase_temp$x <- ifelse(ase_temp$x==2, 0, 1)
    # LM for initialization of GLM
    coverage =  with(ase_temp, a + r)
    y = logit(ase_temp$a / coverage)
    
    get_init = function(x_mat) {
      beta_init = solve(t(x_mat) %*% x_mat, t(x_mat) %*% y) # LM
      # now try to get a sensible initialization for the concentration param
      # BB variance is n*p*(1-p)*(conc+n)/(conc+1).
      # np(1-p) is binomial variance.
      # Second term: (conc+1+(n-1))/(conc+1)=1+(n-1)/(conc+1).
      prob = inv_logit(x_mat %*% beta_init)
      fitted_a = coverage * prob
      var_a = (ase_temp$a - fitted_a) ^ 2
      bin_var = coverage * prob * (1 - prob)
      #(coverage - 1) / (var_a / bin_var - 1) - 1 # method of moments
      conc_init = mean((coverage + 1) / (var_a / bin_var + 1) - 1) # like adding Gamma(2,2) pseudocounts
      if (conc_init < 1.5 | conc_init > 100.)
        conc_init = 10.
      list(beta = as.numeric(beta_init) %>% as.array(),
           conc = conc_init)
    }
    
    # are we testing the exonic SNP itself (or a perfectly linked SNP)
    testing_self_flag = length(unique(ase_temp$het_x)) == 1
    
    # full model with eqtl and gxe
    x_full = if (!testing_self_flag)
      model.matrix(~ x + het_x + x:het_x, data = ase_temp) else
      model.matrix(~ x, data = ase_temp)
    
    stan_dat = list(
      N = nrow(x_full),
      P = ncol(x_full),
      x = x_full,
      ys = ase_temp$a,
      ns = ase_temp$a + ase_temp$r,
      concShape = concShape,
      concRate = concRate
    )
    
    fit_full = if (det(t(x_full) %*% x_full) > 0) {
      optimizing(bb_stan,
                 data = stan_dat)
    } else {
      list(value = NA,
           par = rep(NA, 1 + ncol(x_full)))
    }
    # TODO: does this work when testing the exonic SNP itself?
    names(fit_full$par) = c("conc", "intercept", if (!testing_self_flag)
      "b_eqtl", "b_gxe", "env")
    
    # eqtl but no gxe model
    x_eqtl = if (!testing_self_flag)
      model.matrix(~ het_x , data = ase_temp)
    else
      model.matrix(~ 1, data = ase_temp)
    stan_dat$x = x_eqtl
    stan_dat$N = nrow(x_eqtl)
    stan_dat$P = ncol(x_eqtl)
    fit_eqtl = optimizing(bb_stan,
                          data = stan_dat)$value
    
    # null model, no eQTL or GxE
    x_0 = model.matrix(~ 1, data = ase_temp)
    stan_dat$x = x_0
    stan_dat$N = nrow(x_0)
    stan_dat$P = ncol(x_0)
    fit_0 = optimizing(bb_stan,
                       data = stan_dat)$value
    
    df = ncol(x_full) - ncol(x_eqtl)
    
    # if (F) { # debugging code
    #   ase_temp %>% mutate(het=reg_geno,
    #                       coverage=r+a,
    #                       ar = a/coverage,
    #                       in_phase=geno == reg_geno) %>%
    #     ggplot(aes( x,ar,size=coverage,col=factor(het_x))) + geom_point() + ylim(0,1) +
    #     xlab("Environmental factor") + ylab("Phased allelic ratio")
    #   pchisq( 2.0*(fit_full$value - fit_eqtl), df = df, lower.tail = F)
    #   pchisq( 2.0*(fit_eqtl - fit_0), df = 1, lower.tail = F)
    # }
    
    num_het_ind = length(unique(ase_temp %>% filter(het_x != 0) %>% .$individual)) # for performance analysis later
    
    data_frame(
      total_count = sum(coverage),
      num_het_snps = num_het_snps,
      num_het_ind = num_het_ind,
      reg_snp_pos = snp_pos,
      df = df,
      l0 = fit_0,
      l_geno = fit_eqtl,
      l_interact = fit_full$value) %>% cbind(as_data_frame(as.list(fit_full$par)))
  }
  
  if (!is.null(checkpoint_dir)) {
    print("Saving results")
    checkpoint_file = gzfile(check_fn, "w")
    temp_results %>% write_tsv(checkpoint_file) # write_tsv
    close(checkpoint_file)
  }
  
  if (is.null(temp_results)) {
    return(NULL)
  } else{
    temp_results %>%
      mutate(exonic_snp_pos = exonic_snp_pos)
  }
}


#res_file= gzfile( paste0(outputdir,which_chr,".txt.gz"), "w" )
#allres %>% format(digits=5) %>% write.table(res_file, quote = F, row.names = F, col.names = T, sep="\t")
#close(res_file)
#        allres = ase %>%
#          eagle1S(phased,
#                  meta,
#                  cisdist = cisdist,
#                  checkpoint_dir=checkpoint_dir)

save(allres, file = paste(resDir, which_chr, ".Rdata", sep = ""))
