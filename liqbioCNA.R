#!/usr/bin/env Rscript

###########################################

library(getopt)
library(RJSONIO)
library(VariantAnnotation)
#long, short(NA), argmask, datatype, desc
#argmask 0=no arg, 1=req, 2=optional
args <- rbind(
  c("tumor_cnr", NA, 1, "character", "tumor bin file from CNVkit"),
  c("tumor_cns", NA, 1, "character", "tumor segment file from CNVkit"),
  c("normal_cnr", NA, 1, "character", "normal bin file from CNVkit"),
  c("normal_cns", NA, 1, "character", "normal segment file from CNVkit"),
  c("het_snps_vcf", NA, 1, "character", "heterozygous SNPs .vcf file"),
  c("purecn_csv", NA, 1, "character", "PureCN result .csv file"),
  c("purecn_genes_csv", NA, 1, "character", "PureCN result _genes.csv"),
  c("purecn_loh_csv", NA, 1, "character", "PureCN result _loh.csv"),
  c("purecn_variants_csv", NA, 1, "character", "PureCN result _variants.csv"),
  c("svcaller_T_DEL", NA, 1, "character", "Tumor SV caller DEL-events.gtf"),
  c("svcaller_T_DUP", NA, 1, "character", "Tumor SV caller DUP-events.gtf"),
  c("svcaller_T_INV", NA, 1, "character", "Tumor SV caller INV-events.gtf"),
  c("svcaller_T_TRA", NA, 1, "character", "Tumor SV caller TRA-events.gtf"),
  c("svcaller_N_DEL", NA, 1, "character", "Normal SV caller DEL-events.gtf"),
  c("svcaller_N_DUP", NA, 1, "character", "Normal SV caller DUP-events.gtf"),
  c("svcaller_N_INV", NA, 1, "character", "Normal SV caller INV-events.gtf"),
  c("svcaller_N_TRA", NA, 1, "character", "Normal SV caller TRA-events.gtf"),
  c("germline_mut_vcf", NA, 1, "character", "germline mutation vcf file"),
  c("somatic_mut_vcf", NA, 1, "character", "somatic mutation vcf file"),
  c("plot_png", NA, 1, "character", "plot .png file name"),
  c("cna_json", NA, 1, "character", "CNA output json file name"),
  c("purity_json", NA, 1, "character", "purity output json file name")
)


opts <- getopt(args)


chrsizes=structure(list(
  chr = c("1", "2", "3", "4", "5", "6", "7", "8", 
          "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", 
          "20", "21", "22", "X", "Y", "MT"), 
  size = c(249250621L, 243199373L, 
           198022430L, 191154276L, 180915260L, 171115067L, 159138663L, 146364022L, 
           141213431L, 135534747L, 135006516L, 133851895L, 115169878L, 107349540L, 
           102531392L, 90354753L, 81195210L, 78077248L, 59128983L, 63025520L, 
           48129895L, 51304566L, 155270560L, 59373566L, 16569L), 
  cumstart = c(0, 
               249250621, 492449994, 690472424, 881626700, 1062541960, 1233657027, 
               1392795690, 1539159712, 1680373143, 1815907890, 1950914406, 2084766301, 
               2199936179, 2307285719, 2409817111, 2500171864, 2581367074, 2659444322, 
               2718573305, 2781598825, 2829728720, 2881033286, 3036303846, 3095677412
  )), 
  .Names = c("chr", "size", "cumstart"), 
  row.names = c(NA,-25L), 
  class = "data.frame")

genes=data.frame(label = c("APC", "AR", "ATM", "BRCA1", "BRCA2", 
                           "CCND1", "CDK12", "CDKN1B", "CDKN2A", "CDKN2B", "CHD1", "CHEK2", 
                           "FANCA", "HDAC2", "MLH1", "MLH3", "MSH2", "MSH3", "MSH6", "MYC", 
                           "NKX3-1", "PIK3CA", "PIK3R1", "PMS2", "PPP2R2A", "PTEN", "RB1", 
                           "TMPRSS2", "TP53", "ZBTB16"), 
                 chromosome = c("5", "X", "11", 
                                "17", "13", "11", "17", "12", "9", "9", "5", "22", "16", "6", 
                                "3", "14", "2", "5", "2", "8", "8", "3", "5", "7", "8", "10", 
                                "13", "21", "17", "11"), 
                 start = c(112043201, 66763873, 108093558, 
                           41196311, 32889616, 69455872, 37617738, 12870301, 21967750, 22002901, 
                           98190907, 29083730, 89864766, 114257319, 37034840, 75480466, 
                           47630205, 79950466, 48010220, 128748314, 23536205, 178866310, 
                           67584251, 6012869, 26149006, 89623194, 48877882, 4.1e+07, 7571719, 
                           113931287), 
                 end = c(112181936L, 66950461L, 108239826L, 41277468L, 
                         32973809L, 69469242L, 37690800L, 12875305L, 21975132L, 22009312L, 
                         98262238L, 29137822L, 89883065L, 114291888L, 37092337L, 75518235L, 
                         47710367L, 80172634L, 48034092L, 128753680L, 23540402L, 178952497L, 
                         67597649L, 6048737L, 26230195L, 89728532L, 49056026L, 42880085L, 
                         7590868L, 114121397L),
                 stringsAsFactors = F)
genes$cumstart <- genes$start + chrsizes$cumstart[match(genes$chromosome,chrsizes$chr)]
genes$cumend <- genes$end + chrsizes$cumstart[match(genes$chromosome,chrsizes$chr)]


### Read SNP allele ratio
{
  vcf <- readVcf(opts$het_snps_vcf,genome = "GRCh37")
  g <- geno(vcf)
  chr <- pos <- rownames(g$DP)
  alf=NULL
  smoothedAi=NA
  if (length(pos)>0) { 
    for (i in 1:length(pos)) {
      temp <- strsplit(chr[i],':')[[1]]
      chr[i] <- as.character(temp[1])
      pos[i] <- as.numeric(strsplit(temp[2],'_')[[1]][1])
    }
    pos=as.numeric(pos)
    alf <- data.frame(chromosome=chr, start=pos, end=pos, stringsAsFactors = F, cumstart=NA, cumend=NA)
  }
  
  if (!is.null(alf)) if( ! all(is.na(alf$chromosome)) ) {
    
    for(chr in chrsizes$chr){
      ix <- which(alf$chromosome == chr)
      alf$cumstart[ix] <- alf$start[ix] + chrsizes$cumstart[chrsizes$chr==chr]
      alf$cumend[ix] <- alf$end[ix] + chrsizes$cumstart[chrsizes$chr==chr]
    }
    
    alf$t <- as.numeric(g$AO[,2])/as.numeric(g$DP[,2])
    alf$t[is.nan(alf$t)]=NA # allele freq becomes NaN if cov=0. Then set to NA
    alf$n <- as.numeric(g$AO[,1])/as.numeric(g$DP[,1])
    alf$n[is.nan(alf$n)]=NA # allele freq becomes NaN if cov=0. Then set to NA
    alf$td <- as.numeric(g$DP[,2])
    alf$nd <- as.numeric(g$DP[,1])
    
    alf$ai=2*abs(alf$t-0.5)
  }
}


### Read somatic point mutations:
{
  vcf <- readVcf(opts$somatic_mut_vcf,genome = "GRCh37")
  g <- geno(vcf)
  r=rowRanges(vcf)
  if (length(g)>0) { # if there are any somatic mutations...
    chr=as.character(seqnames(r))
    pos=data.frame(ranges(r))$start
    salf <- data.frame(N=1:length(chr),chromosome=chr,pos=pos,stringsAsFactors = F)
    rownames(salf)=names(r)
    salf$REF=as.data.frame(r$REF)[,1]
    salf$ALT=as.data.frame(r$ALT)[,3]
    salf$cumpos <- NA
    for(chr in chrsizes$chr){
      ix <- which(salf$chromosome == chr)
      salf$cumpos[ix] <- salf$pos[ix] + chrsizes$cumstart[chrsizes$chr==chr]
    }
    
    salf$AF.T <- as.numeric(g$VD[,1])/as.numeric(g$DP[,1])
    salf$AO.T <- as.numeric(g$VD[,1])
    salf$DP.T <- as.numeric(g$DP[,1])
    salf$AO.N <- as.numeric(g$VD[,2])/as.numeric(g$DP[,2])
    salf$DP.N <- as.numeric(g$DP[,2])
    
    salf$type='other'
    salf$type[isSNV(vcf)]='snv'
    salf$type[isDeletion(vcf)]='del'
    salf$type[isInsertion(vcf)]='ins'
    
    
    header=info(header(vcf))$Description
    ix=grep('Consequence annotations from Ensembl',header)
    header=strsplit(header[ix],'\\|')[[1]]
    header[1]='Allele'

    vep=info(vcf)$CSQ
    
    ## This loop creates a new "table" with all mutation effects.
    rowspermut=unlist(lapply(vep,length))
    table=#data.frame(
      matrix('',nrow = sum(rowspermut),ncol = length(header)+1)#,stringsAsFactors = F)
    colnames(table)=c('N',header)  # blir detta fel ibland?????
    for (i in 1:length(vep)) {
      #if (i %% 100 ==0) cat(i,'..')
      for (j in 1:length(vep[[i]])) { # for each effect
        t2=vep[[i]][j]
        t2=strsplit(t2,'[|]')[[1]] # pipe separated line with one effect of the mutation
        t2=c((i),t2)
        thisrow=sum(rowspermut[1:i])-rowspermut[i]+j
        table[thisrow,1:length(t2)]=t2
      }
    }
    
    salf=merge(salf,table,by='N',all=T)
  } #end somatic mutations
  salf=(salf[,-1])
  
  # mark the type
  salf$pch=rep(0,nrow(salf))
  salf$pch[salf$type=='snv']=21
  salf$pch[salf$type=='del']=24
  salf$pch[salf$type=='ins']=25
  
  # Icke-NA/intron på konsekvens
  salf$hasConsequence=!salf$Consequence %in% c("intron_variant", "synonymous_variant", 
                                             "splice_region_variant&intron_variant", 
                                             "3_prime_UTR_variant", "intergenic_variant", 
                                             "regulatory_region_variant", "upstream_gene_variant", 
                                             "downstream_gene_variant", 
                                             "intron_variant&non_coding_transcript_variant", 
                                             "5_prime_UTR_variant", "splice_region_variant&synonymous_variant", 
                                             "non_coding_transcript_exon_variant&non_coding_transcript_variant", 
                                             "intron_variant&NMD_transcript_variant", "TF_binding_site_variant", 
                                             "splice_region_variant&non_coding_transcript_exon_variant&non_coding_transcript_variant", 
                                             "mature_miRNA_variant", "coding_sequence_variant&5_prime_UTR_variant")
  salf$hasConsequence[is.na(salf$hasConsequence)]=F

    
  ## Vilka är tydligt deleterious? (inkl inframe)
  salf$deleterious=rep(F,nrow(salf))
  salf$deleterious[unique(c(
    grep('start',salf$Consequence),
    grep('stop',salf$Consequence),
    grep('frame',salf$Consequence),
    grep('acceptor',salf$Consequence),
    grep('donor',salf$Consequence)
  ))]=T
  salf$deleterious[!salf$hasConsequence]=F
  
}

### Read germline point mutations
{
  vcf <- readVcf(opts$germline_mut_vcf,genome = "GRCh37")
  g <- geno(vcf)
  r=rowRanges(vcf)
  if (length(g)>0) { # if there are any somatic mutations...
    chr=as.character(seqnames(r))
    pos=data.frame(ranges(r))$start
    galf <- data.frame(N=1:length(chr),chromosome=chr,pos=pos,stringsAsFactors = F)
    rownames(galf)=names(r)
    galf$cumpos <- NA
    galf$REF=as.data.frame(r$REF)[,1]
    galf$ALT=as.data.frame(r$ALT)[,3]
    for(chr in chrsizes$chr){
      ix <- which(galf$chromosome == chr)
      galf$cumpos[ix] <- galf$pos[ix] + chrsizes$cumstart[chrsizes$chr==chr]
    }
    
    galf$AF <- as.numeric(g$AO[,1])/as.numeric(g$DP[,1])
    galf$AO <- as.numeric(g$AO[,1])
    galf$DP <- as.numeric(g$DP[,1])
    
    galf$type='other'
    galf$type[isSNV(vcf)]='snv'
    galf$type[isDeletion(vcf)]='del'
    galf$type[isInsertion(vcf)]='ins'
    
    
    header=info(header(vcf))$Description
    ix=grep('Consequence annotations from Ensembl',header)
    header=strsplit(header[ix],'\\|')[[1]]
    header[1]='Allele'
    
    vep=info(vcf)$CSQ
    
    ## This loop creates a new "table" with all mutation effects.
    rowspermut=unlist(lapply(vep,length))
    table=#data.frame(
      matrix('',nrow = sum(rowspermut),ncol = length(header)+1)#,stringsAsFactors = F)
    colnames(table)=c('N',header)  # blir detta fel ibland?????
    for (i in 1:length(vep)) {
      for (j in 1:length(vep[[i]])) { # for each effect
        t2=vep[[i]][j]
        t2=strsplit(t2,'[|]')[[1]] # pipe separated line with one effect of the mutation
        t2=c((i),t2)
        thisrow=sum(rowspermut[1:i])-rowspermut[i]+j
        table[thisrow,1:length(t2)]=t2
      }
    }
    
    galf=merge(galf,table,by='N',all=T)
  } #end germline mutations
  galf=(galf[galf$ALT>=12 & galf$AF>=0.2,-1])
  # mark the type
  galf$pch=rep(0,nrow(galf))
  galf$pch[galf$type=='snv']=21
  galf$pch[galf$type=='del']=24
  galf$pch[galf$type=='ins']=25
  
  # Icke-NA/intron på konsekvens
  galf$hasConsequence=!galf$Consequence %in% c("intron_variant", "synonymous_variant", 
                                               "splice_region_variant&intron_variant", 
                                               "3_prime_UTR_variant", "intergenic_variant", 
                                               "regulatory_region_variant", "upstream_gene_variant", 
                                               "downstream_gene_variant", 
                                               "intron_variant&non_coding_transcript_variant", 
                                               "5_prime_UTR_variant", "splice_region_variant&synonymous_variant", 
                                               "non_coding_transcript_exon_variant&non_coding_transcript_variant", 
                                               "intron_variant&NMD_transcript_variant", "TF_binding_site_variant", 
                                               "splice_region_variant&non_coding_transcript_exon_variant&non_coding_transcript_variant", 
                                               "mature_miRNA_variant", "coding_sequence_variant&5_prime_UTR_variant")
  galf$hasConsequence[is.na(galf$hasConsequence)]=F
  
  
  ## Vilka är tydligt deleterious? (inkl inframe)
  galf$deleterious=rep(F,nrow(galf))
  galf$deleterious[unique(c(
    grep('start',galf$Consequence),
    grep('stop',galf$Consequence),
    grep('frame',galf$Consequence),
    grep('acceptor',galf$Consequence),
    grep('donor',galf$Consequence)
  ))]=T
  galf$deleterious[!galf$hasConsequence]=F
  
  
  
}

### read CNVkit copy number data
{
  segments <- read.delim(opts$tumor_cns,stringsAsFactors = F)
  bins <- read.delim(opts$tumor_cnr,stringsAsFactors = F)

  ## Segment start and end pos
  segments$cumstart <- NA
  segments$cumend   <- NA
  for(chr in chrsizes$chr){
    ix <- which(segments$chromosome == chr)
    segments$cumstart[ix] <- segments$start[ix] + chrsizes$cumstart[chrsizes$chr==chr]
    segments$cumend[ix] <- segments$end[ix] + chrsizes$cumstart[chrsizes$chr==chr]
  }
  segments$centerpos <- segments$cumstart+(segments$cumend-segments$cumstart)/2
  
  ## Bin start and end pos
  bins$cumstart <- NA
  bins$cumend   <- NA
  for(chr in chrsizes$chr){
    ix <- which(bins$chromosome == chr)
    bins$cumstart[ix] <- bins$start[ix] + chrsizes$cumstart[chrsizes$chr==chr]
    bins$cumend[ix] <- bins$end[ix] + chrsizes$cumstart[chrsizes$chr==chr]
  }
  bins$centerpos <- bins$cumstart+(bins$cumend-bins$cumstart)/2
  bins=bins[order(bins$cumstart),]
}

### Read structural variant files
{ # tumor
  t_strvs=NULL
  try( { 
    sv <- read.delim(opts$svcaller_T_DEL,header=F,stringsAsFactors = F)
    colnames(sv)[c(1,4,5)]=c('chr','start','end')
    sv$type='DEL'
    t_strvs=rbind(t_strvs,sv)
  }, silent=T)
  try( { 
    sv <- read.delim(opts$svcaller_T_DUP,header=F,stringsAsFactors = F)
    colnames(sv)[c(1,4,5)]=c('chr','start','end')
    sv$type='DUP'
    t_strvs=rbind(t_strvs,sv)
  }, silent=T)
  try( { 
    sv <- read.delim(opts$svcaller_T_INS,header=F,stringsAsFactors = F)
    colnames(sv)[c(1,4,5)]=c('chr','start','end')
    sv$type='INS'
    t_strvs=rbind(t_strvs,sv)
  }, silent=T)
  try( { 
    sv <- read.delim(opts$svcaller_T_TRA,header=F,stringsAsFactors = F)
    colnames(sv)[c(1,4,5)]=c('chr','start','end')
    sv$type='TRA'
    t_strvs=rbind(t_strvs,sv)
  }, silent=T)
  t_strvs$cumstart=NA
  t_strvs$cumend=NA
  for(i in 1:nrow(chrsizes)){
    ix <- which(t_strvs$chr == chrsizes$chr[i])
    t_strvs$cumstart[ix] <- t_strvs$start[ix] + chrsizes$cumstart[i]
    t_strvs$cumend[ix] <- t_strvs$end[ix] + chrsizes$cumstart[i]
  }
}
{ # normal
  n_strvs=NULL
  try( { 
    sv <- read.delim(opts$svcaller_N_DEL,header=F,stringsAsFactors = F)
    colnames(sv)[c(1,4,5)]=c('chr','start','end')
    sv$type='DEL'
    n_strvs=rbind(n_strvs,sv)
  }, silent=T)
  try( { 
    sv <- read.delim(opts$svcaller_N_DUP,header=F,stringsAsFactors = F)
    colnames(sv)[c(1,4,5)]=c('chr','start','end')
    sv$type='DUP'
    n_strvs=rbind(n_strvs,sv)
  }, silent=T)
  try( { 
    sv <- read.delim(opts$svcaller_N_INS,header=F,stringsAsFactors = F)
    colnames(sv)[c(1,4,5)]=c('chr','start','end')
    sv$type='INS'
    n_strvs=rbind(n_strvs,sv)
  }, silent=T)
  try( { 
    sv <- read.delim(opts$svcaller_N_TRA,header=F,stringsAsFactors = F)
    colnames(sv)[c(1,4,5)]=c('chr','start','end')
    sv$type='TRA'
    n_strvs=rbind(n_strvs,sv)
  }, silent=T)
  n_strvs$cumstart=NA
  n_strvs$cumend=NA
  for(i in 1:nrow(chrsizes)){
    ix <- which(n_strvs$chr == chrsizes$chr[i])
    n_strvs$cumstart[ix] <- n_strvs$start[ix] + chrsizes$cumstart[i]
    n_strvs$cumend[ix] <- n_strvs$end[ix] + chrsizes$cumstart[i]
  }
}


### Read PureCN files:
{
  purecn_stat=read.delim(opts$purecn_csv,sep = ',',stringsAsFactors = F) # purity/ploidy
  purecn_vars=read.delim(opts$purecn_variants_csv,sep = ',',stringsAsFactors = F) # mutations and snps
  purecn_genes=read.delim(opts$purecn_genes_csv,sep = ',',stringsAsFactors = F) # gene copy number and LOH
  purecn_loh=read.delim(opts$purecn_loh_csv,sep = ',',stringsAsFactors = F) # segmented copy number and LOH
  purecn_loh$cumstart=NA
  purecn_loh$cumend=NA
  for(i in 1:nrow(chrsizes)){
    ix <- which(purecn_loh$chr == chrsizes$chr[i])
    purecn_loh$cumstart[ix] <- purecn_loh$start[ix] + chrsizes$cumstart[i]
    purecn_loh$cumend[ix] <- purecn_loh$end[ix] + chrsizes$cumstart[i]
  }
}

### Calculate/extract purity
{
  # select high-confidence mutations for purity estimate (not AR)
  ix=salf$AO.T>=12 & salf$hasConsequence & salf$SYMBOL %in% c("ATM", "BRCA1", "BRCA2", "CCND1", "CDK12", "CDKN2A", 
                                             "CDKN2B", "CHD1", "CHEK2", "FANCA", "HDAC2", "MYC", "NKX3-1", 
                                             "PIK3CA", "PIK3R1", "PPP2R2A", "PTEN", "RB1", "TMPRSS2", "TP53", 
                                             "ZBTB16","SPOP","MED12","PIK3CA","FOXA1","COL5A1")
  median_af=median(salf$AF.T[ix],na.rm = T)
  # consider the mutation(s) being present @ 1 copy and there being 2 normal copies, what's the purity?
  t=median_af
  n=(1-median_af)/2
  p=t/(t+n)
  
  purity=purecn_stat
  purity$mutation_Purity=round(p,2)
  
  exportJson <- toJSON(purity)
  write(exportJson, opts$purity_json)
  
}

## Collect copy number for selected genes
cn_calls=rep('NO_CALL',nrow(genes))
names(cn_calls) = genes$label
## Get the AR copy number
g=genes[1,]
data_line=data.frame(g)
cbins=bins[bins$chromosome==g$chromosome,] # this chromosome
left=cbins$end<g$start-3e6 # left control
left.median=median(cbins$log2[left],na.rm=T)
target=cbins$start>g$start-.1e6 & cbins$end<g$end+.1e6
target.median=median(cbins$log2[target],na.rm=T)
right=cbins$start>g$end+3e6 # right control
right.median=median(cbins$log2[right],na.rm=T)
target.diff=target.median-c(left.median,right.median) # difference compared to both controls
control.toUse=which.min(abs(target.diff)) # which control is most similar? left:1 and right:2
controlpoints=list(cbins$log2[left],cbins$log2[right]) # list of left and right control points
pval=wilcox.test(x = cbins$log2[target],y = controlpoints[[control.toUse]])$p.value # p value with correct control points
data_line$cnvkit.logr=target.median
data_line$cnvkit.to.control=round(target.diff[control.toUse],3)
data_line$cnvkit.pval=round(pval,5)
# For AR amp, require the AR signal to exceed control by 0.5 (there is no PureCN data for AR)
if (data_line$cnvkit.to.control > 0.5) cn_calls['AR']='AMPLIFIED'

## Check for amplifications
for (g in c('CCND1','MYC', 'PIK3CA')) {
  t=purecn_genes[purecn_genes$gene.symbol==g,]
  if (t$focal & t$C>3) cn_calls[g]='AMPLIFIED' #   calls amplified if focal and >3 copies
  if (t$C >= 7) cn_calls[g]='AMPLIFIED'            #  also if ≥7 copies
  if (!is.na(t$type)) if (t$type=='AMPLIFIED') cn_calls[g]='AMPLIFIED'   #           or if >5 copies regardless of focal status
}
## Check for focal deletion and LOH of remaining
for (g in c("APC", "ATM", "BRCA1", "BRCA2", "CCND1", "CDK12", "CDKN1B", 
            "CDKN2A", "CDKN2B", "CDH1", "CHEK2", "FANCA", "HDAC2", "MLH1", 
            "MLH3", "MSH2", "MSH3", "MSH6", "NKX3-1", "PIK3R1", 
            "PMS2", "PPP2R2A", "PTEN", "RB1", "TMPRSS2", "TP53", "ZBTB16")) {
  t=purecn_genes[purecn_genes$gene.symbol==g,] 
  if (nrow(t)>0){
    if (isTRUE(t$loh)) cn_calls[g]='LOSS_OF_HETEROZYGOSITY'
    if (t$focal & t$C==1) cn_calls[g]='FOCAL_DELETION'
    if (t$C==0) cn_calls[g]='HOMOZYGOUS_DELETION'
  }
}
## TMPRSS2 fusion is special case
g=genes[genes$label=='TMPRSS2',]
data_line=data.frame(g)
cbins=bins[bins$chromosome==g$chromosome,] # this chromosome
left=cbins$end<g$start # left control
left.median=median(cbins$log2[left],na.rm=T)
target=cbins$start>g$start-.1e6 & cbins$end<g$end+.1e6
target.median=median(cbins$log2[target],na.rm=T)
right=cbins$start>g$end # right control
right.median=median(cbins$log2[right],na.rm=T)
target.diff=target.median-c(left.median,right.median) # difference compared to both controls
control.toUse=which.min(abs(target.diff)) # which control is most similar? left:1 and right:2
controlpoints=list(cbins$log2[left],cbins$log2[right]) # list of left and right control points
pval=wilcox.test(x = cbins$log2[target],y = controlpoints[[control.toUse]])$p.value # p value with correct control points
data_line$cnvkit.logr=target.median
data_line$cnvkit.to.control=round(target.diff[control.toUse],3)
data_line$cnvkit.pval=round(pval,5)
# For TMPRSS-ERG del, require the signal to be 0.1 below nearest control
if (data_line$cnvkit.to.control < -.01) cn_calls['TMPRSS2']='IMPLIED_FUSION'

# Write copy numbers to JSON file
exportJson <- toJSON(cn_calls)
write(exportJson, opts$cna_json)


### Frankenplot code
{
  try( {
    # Prepare the per-SNP logratio
    alf$log2 <- NA
    delta <- 3e6; for (i in 1:nrow(alf)) {
      ix <- bins$chromosome==alf$chromosome[i] & bins$start>alf$start[i]-delta & bins$end<alf$start[i]+delta
      alf$log2[i]=median(bins$log2[ix],na.rm=T)
    }
    # Prepare the smoothed AI
    ai <- smoothedAi <- alf$ai
    for (i in 1:length(ai)) {
      ss=max(1,i-5); e=min(i+5,length(ai))
      smoothedAi[i]=median(ai[ss:e])
    }
  }, silent=T)
  
  ## File name here:
  png(filename = opts$plot_png,width=11.7,height=8.3,units="in",res=600)
  
  ## set screens
  split.screen(figs=c(2,1))
  split.screen(as.matrix(data.frame(
    left=c(0,0.25),
    right=c(0.2,.9635),
    bottom=c(0,0.097),
    top=c(1,0.903))),1)
  split.screen(figs=c(2,1),3)
  split.screen(as.matrix(data.frame(left=c(rep(0.015,3)),
                                    right=c(rep(1,3)),
                                    bottom=c(0.10,0.5,0.55),
                                    top=c(0.5,0.55,1))),2)
  split.screen(figs=c(4,6),4)
  
  
  screen(1)
  # Top margin with text
  try({
    plot(1,type='n',axes=F,xlab='',ylab='')
    snpcov='NA'; if (!is.null(alf)) if (nrow(alf)>0) 
      snpcov=paste0('T',round(median(alf$td)),'/N',round(median(alf$nd)))
    mtext(
      text = paste0(purity$Sampleid,' ----- ', format(Sys.time(), "%F %H:%M:%S"),'',
                    '  SNPcov:',snpcov,
                    '  purCn.Ploidy:',round(purity$Ploidy,2),
                    '  purCn.Purity:',purity$Purity,
                    '  smut.Purity:',purity$mutation_Purity),
      side = 3,padj=-5.5,adj=1,cex=0.7)
  },silent=T)
  
  cex.axis <- .6
  cex.mtext <- 1.5
  cex.main <- 2
  cex.text <- .7
  axis1padj <- -2.35
  axis2hadj <- 0.3
  text1padj <- 1.7
  text2padj<- -3
  padj <- -3
  lwd=1
  color <- '#00000010'
  
  
  
  screen(5)
  #plot snp cov vs alf
  xlim=c(1,3000)
  ylim=(0:1)
  try( {
    par(mar=c(2,3,2,.5),xaxs='i',yaxs='i',las=1)
    plot(1,type='n',xlim=xlim,ylim=ylim,axes=F,main='Variants',cex.main=0.5,log='x')
    text(x=c(2,4,6),y=0.48,labels = c('1 alt','2','3'),cex=0.6,col='grey',srt=-55)
    points(x = 7:1000,y=3/(7:1000),type='l',col='grey')
    axis(1,c(1,10,100,1000),lwd=lwd,lend=1,cex.axis=cex.axis,padj=axis1padj,tck=-0.03)
    axis(2,seq(0,1,.25),lwd=lwd,lend=1,cex.axis=cex.axis,hadj=axis2hadj,tck=-0.03)
    mtext('Coverage',1,cex=cex.text,padj=text1padj)
    mtext('Alt allele ratio',2,cex=cex.text,las=0,padj=text2padj)
    points(alf$td,alf$t,cex=0.3,col='#00000080',xlim=xlim,ylim=ylim,pch=16,lwd=lwd)
    points(alf$nd,alf$n,cex=0.1,col='#60606080',xlim=xlim,ylim=ylim,pch=3,lwd=lwd)
    scol=rep('#C00000CC',nrow(salf))
    scol[salf$AO.T<6]='#500000CC'
    points(salf$DP.T,salf$AF.T,cex=0.4,col=scol,xlim=xlim,ylim=ylim,pch=salf$pch,lwd=lwd)
    segments(x0=median(alf$td,na.rm = T),y0=0,x1=median(alf$td,na.rm = T),y1=1,col='#00000080',lwd=lwd,lty=3)
    segments(x0=median(alf$nd,na.rm = T),y0=0,x1=median(alf$nd,na.rm = T),y1=1,col='#7070FF80',lwd=lwd,lty=3)
  }, silent=T)
  
  
  
  
  screen(6)
  try( {
    #plot AR only.
    #xlim=c(55e6,80e6)
    ylim=c(-1.2,2)
    
    ar=genes[1,]
    par(mar=c(2,3,2,.5),xaxs='i',yaxs='i',las=1)
    ix=bins$chromosome=='X' & bins$gene!='Background'
    t=bins[ix,]
    plot(t$log2,ylim=ylim,pch=16,cex=0.3,main='AR',cex.main=0.6,axes=F)
    axis(2,c(-1,0,1,2),lwd=lwd,lend=1,cex.axis=cex.axis,hadj=axis2hadj,tck=-0.03)
    points(t$log2,type='l',col='#00000050')
    # rect(xleft = c(66763873, 66863097, 66905851, 66931243, 66937319, 66941674, 
    #                66942668, 66943527),
    #      xright = c(66766604, 66863249, 66905968, 66931531, 66937464, 66941805, 
    #                 66942826, 66950461),
    #      ybottom = -0.1,ytop = 0.1)
    mtext('Log ratio',2,cex=cex.text,las=0,padj=text2padj)
  }, silent=T)
  
  
  
  #allelefreq vs logR
  
  screen(4)
  try( {
    par(mar=c(0,0,0,0))
    plot(1,type='n',axes=F,xlab='',ylab='')
    mtext('DNA ratio',1,padj=text1padj,cex=cex.text)
    mtext('Allelic imablance',2,padj=-3,cex=cex.text)
    
    for (c in 1:24)
    {
      cex=0.4
      screen(c+9) # for chrY, screen is 33.
      par(mar = c(0, 0, 0, 0))
      par(oma = c(0,0,0,0))
      par(mgp =c(1,0.5,0))
      #xlim=c(-1.1,1.1)
      xlim=c(0.3,2.1)
      ylim=c(0,1)
      ix <- alf$chromosome %in% c(chrsizes$chr[c],'X','Y')
      ixCurChr <-  alf$chromosome %in% chrsizes$chr[c]
      plot(2^alf$log2[!ix],smoothedAi[!ix],xlim=xlim,ylim=ylim,lwd=lwd,axes=F,ylab='',xlab='',pch=16,
           col='#B0B0B030',cex=cex)
      points(2^alf$log2[ixCurChr],smoothedAi[ixCurChr],col='#00800070',pch=16,cex=cex)
      
      
      if (c==24 & !is.null(alf)) {
        d=density(2^alf$log2[alf$chromosome %in% as.numeric(1:22)])
        points(d$x,d$y/max(d$y),type='l')
      }
      
      whole=(c(0.5,1,1.5,2))
      
      segments(
        x0=whole,x1=whole,
        y0=-0.05,y1=1.035,
        col='#D3D3D360',
        lwd=1)
      
      segments(
        x0=c(-2),x1=c(5),
        y0=c(1/3,1/2),y1=c(1/3,1/2),
        col='#D3D3D360',
        lwd=1)
      
      if(c == 23)
      {
        mtext("X", side = 3, line = -1, adj = 0.92, cex = 1)
      } else if(c == 24)
      {
        mtext("", side = 3, line = -1, adj = 0.92, cex = 1)
      } else
      {
        mtext(c, side = 3, line = -1, adj = 0.92, cex = 1)
      }
      
      if(c %in% 19:24) axis(side=1,cex.axis=0.5,at=c(0.5,1,1.5,2),tck=-0.05,padj=-1.6,#las=3,
                            labels=c('-50%','±0','+50%','+100%'),
                            col='white',col.ticks='black',lend=1)
      if(c %in% c(1,7,13,19)) axis(side=2,cex.axis=0.6,tck=-0.05,
                                   at=c(1/3,1/2),
                                   labels=c('2:1','3:1'),
                                   las=1,col='white',col.ticks='black',lend=1)
      if(c %in% c(6,12,18,24)) axis(side=4,cex.axis=0.6,tck=-0.05,
                                    at=c(1/3,1/2),
                                    labels=c('2:1','3:1'),
                                    las=1,col='white',col.ticks='black',lend=1)
      box(lwd=1)
    }
  }, silent=T)
  
  
  
  
  #LogR plot
  screen(9)
  try( {
    #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
    par(mar = c(0, 0, 0, 0))
    par(oma = c(0,0,0,0))
    par(mgp =c(1,0.5,0))
    par(lend=1)
    
    ymin = -2
    ymax = 2
    #Plot signal over complete genome
    plot(NA,NA,#mpos[ix],mval[ix],
         pch=16,
         cex=0.3,
         main='',
         xlab = "",
         ylab = "",
         col = '#00000003',
         xaxt="n",
         axes=F,
         ylim = c(ymin,ymax),
         xlim = c(0,3095370729)
    )
    seqminmax <- seq(ymin,ymax,by=0.5)
    #Add axis to the left & right of signal
    axis(side=2,tck=-0.025,at=seqminmax,cex.axis=0.6,pos=0,las=1)
    axis(side=4,tck=-0.025,at=seqminmax,cex.axis=0.6,pos=3095370729,las=1)
    mtext("Log ratio",side=2,line=0,cex=0.7,padj=1.15)
    
    whole=(c(.5,1,1.5,2))
    #Add grey segments
    segments(
      y0=log2(whole),y1=log2(whole),
      x0=0,x1=3e9,
      col='#00000020',
      lwd=1)
    
    
    #Add genes as lines
    segments(
      x0=(genes$cumstart+genes$cumend)/2,
      y0=-100,y1=100,
      col='#0000C020',
      lwd=2)
    
    if (!is.null(bins)) {
      ix=bins$gene=='Background'
      points((bins$cumstart[!ix]+bins$cumend[!ix])/2,bins$log2[!ix], #ontargets
             pch=1,
             cex=0.7,
             type='l',
             col = '#00000020'
      )
      points((bins$cumstart[!ix]+bins$cumend[!ix])/2,bins$log2[!ix], #ontargets
             pch=16,
             cex=0.5,
             #type='l',
             col = '#00000080'
      )
      
      #Add segments
      col='#00C000CC'
      segments(x0=segments$cumstart,x1=segments$cumend,
               y0=segments$log2,y1=segments$log2,
               col=col,
               lwd=3)
      #Add segments from PureCN 
      segments(x0=0,x1=3e9,
               y0=-1.8,y1=-1.8,
               col='#00000010',
               lwd=5)
      ix=purecn_loh$M==0 | purecn_loh$C==1
      if (sum(ix)>0) segments(x0=purecn_loh$cumstart[ix],x1=purecn_loh$cumend[ix],
                              y0=-1.8,y1=-1.8,
                              col='cyan',
                              lwd=5)
      col=rep(NA,nrow(purecn_loh))
      col[purecn_loh$C==1]='blue'
      col[purecn_loh$C==0]='violet'
      col[purecn_loh$C==3]='red'
      col[purecn_loh$C>=4]='orange'
      segments(x0=purecn_loh$cumstart,x1=purecn_loh$cumend,
               y0=-1.8,y1=-1.8,
               col=col,
               lwd=3)
      ix=purecn_loh$C==0
      if (sum(ix)>0) points(x = (purecn_loh$cumstart[ix]+purecn_loh$cumend[ix])/2,y = -1.8,pch=24,bg='lightblue')
    }
    
    #Add a bar between chromosomes to distinguish them
    segments(
      x0=chrsizes$cumstart+chrsizes$size,
      y0=-100,y1=100,
      col='#00000099',
      lwd=1)
  }, silent=T)
  
  ## add any structural variants
  try ( {
    text(x = (t_strvs$cumstart+t_strvs$cumend)/2,y = 1.8,labels = t_strvs$type,col='red',cex=0.4,srt=90)
    text(x = (n_strvs$cumstart+n_strvs$cumend)/2,y = 1.8,labels = n_strvs$type,col='blue',cex=0.4,srt=90)
    }, silent=T)
  

  
  #Genes
  screen(8)
  try( {
    #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
    par(mar = c(0, 0, 0, 0))
    par(oma = c(0,0,0,0))
    par(mgp =c(1,0.5,0))
    
    #Create a empty plot with the same ylim and xlim as the plot directly above it (signal) and below (AI)
    plot(0,0,xlab="",ylab="",main="",type="n",axes=F,xaxt="n",ylim=c(0,1),xlim=c(0,3095370729))
    
    chrsizes$labelpos[chrsizes$chr=='MT']=NA
    text(x = chrsizes$cumstart+0.5*chrsizes$size,y = 0.5,labels = chrsizes$chr,cex=.5)
    
    mtext(text='',side=2,las=1,line=-1.2)
    
    # Gene names
    text(x = (genes$cumstart+genes$cumend)/2,y = 0.5,labels = genes$label,srt=45,cex=0.4,col='#0000C0CC')
    
  }, silent=T)
  
  
  
  
  # BAF plot
  screen(7)
  
  try( {
    #Set marginals, outer marginals and mgp(which is for xlab,ylab,ticks and axis)
    par(mar = c(0, 0, 0, 0))
    par(oma = c(0,0,0,0))
    par(mgp =c(1,0.5,0))
    
    
    #Plot normal BAF over whole genome
    plot(NA,#alf$cumstart,alf$n,
         pch=16,
         cex=0.6,
         cex.axis=1,
         main='',
         xlab = "",
         ylab = "",
         col = "#80808080",
         xaxt="n",
         axes=F,
         ylim = c(-0.1,1.1),
         xlim = c(0,3095370729)
    )
    #Add axis to the left,right and below of AI. The below axis is the chromosome numbers 1-24.
    axis(side=2,tck=-0.04,at=seq(0,1,.2),cex.axis=0.6,pos=0,las=1) #at=c(0,0.25,0.33,0.5,0.67,0.75,1),labels=c('0','1/4','1/3','1/2','2/3','3/4','1')
    #axis(side=1,at=pre,pos=0,labels=c(seq(from=1,to=22),"X",'Y'),cex.axis=0.50,lty=0)#,tck=0,col.ticks='#00000000')
    axis(side=4,tck=-0.04,at=seq(0,1,.2),cex.axis=0.6,pos=3095370729,las=1) #at=c(0,0.25,0.33,0.5,0.67,0.75,1),labels=c('0','1/4','1/3','1/2','2/3','3/4','1')
    mtext("SNP allele ratio",side=2,line=0,cex=0.7,padj = 1.15)
    
    #Add a bar between chromosomes to distinguish them
    segments(
      x0=chrsizes$cumstart+chrsizes$size,
      y0=-100,y1=100,
      col='#00000099',
      lwd=1)
    
    #Add genes as lines
    segments(
      x0=(genes$cumstart+genes$cumend)/2,
      y0=-100,y1=100,
      col='#0000C020',
      lwd=2)
    
    segments(
      y0=c(0,1/4,1/3,2/3,3/4,1),y1=c(0,1/4,1/3,2/3,3/4,1),
      x0=0,x1=3095370729,
      # col='#D3D3D360',
      col='#00000020',
      lwd=1)
    
    points(alf$cumstart,alf$t,
           pch=16,
           cex=0.6,
           col = "#00000040")
    
    
    
  }, silent=T)
  try( {
    
    
    ## Add somatic mutations
    if (nrow(salf)>0) {
      scol=rep('#C00000CC',nrow(salf))
      ix=salf$AF.T<0.02 | salf$AO.T<6
      scol[ix]='#500000CC'
      points(salf$cumpos[ix],salf$AF.T[ix],
             pch=salf$pch[ix],
             cex=0.6,
             bg=scol[ix]
      )
      points(salf$cumpos[!ix],salf$AF.T[!ix],
             pch=salf$pch[!ix],
             cex=0.6,
             bg=scol[!ix]
      )
      
      g <- pos <- aa <- rep('',nrow(salf))
      for (j in 1:nrow(salf)) {
        g[j]=as.character(salf$SYMBOL[j])
        g[is.na(g)]=''
        pos[j]=strsplit(as.character(salf$Protein_position[j]),'/')[[1]][1]
        aa[j]=strsplit(as.character(salf$Amino_acids[j]),'/')[[1]][2]
        pos[is.na(pos)]=''; pos[pos=='-']=''; aa[is.na(aa)]=''
      }
      ix=pos!='' & g %in% genes$label  & aa!=''
        if (sum(ix)>0) text(x = salf$cumpos[ix],y=salf$AF.T[ix]-0.07,labels = paste0(g[ix],' ',pos[ix],aa[ix]),cex=0.4,srt=30,col=scol[ix])
    }
    
  }, silent=T)
  try( {
    
    ## Add germline mutations
    scol=rep('#0000C0CC',nrow(galf))
    
    ix=(galf$hasConsequence) & galf$ALT>=12 & galf$AF>=0.2
    
    points(galf$cumpos[ix],galf$AF[ix],
           pch=galf$pch[ix],
           cex=0.6,
           bg=scol[ix]
    )
    if (nrow(galf)>0) {
      g <- pos <- aa <- rep('',nrow(galf))
      for (j in 1:nrow(galf)) {
        g[j]=as.character(galf$SYMBOL[j])
        g[is.na(g)]=''
        pos[j]=strsplit(as.character(galf$Protein_position[j]),'/')[[1]][1]
        aa[j]=strsplit(as.character(galf$Amino_acids[j]),'/')[[1]][2]
        pos[is.na(pos)]=''; pos[pos=='-']=''; aa[is.na(aa)]=''
      }
      ix=ix & pos!='' & g %in% genes$label & aa!=''
      text(x = galf$cumpos[ix],y=galf$AF[ix]+0.07,labels = paste0(g[ix],' ',pos[ix],aa[ix]),cex=0.4,srt=30,col=scol[ix])
    }
    
  }, silent=T)
  
  #Close all the opened split.screens and release the figure
  
  close.screen(all.screens=T)
  dev.off()
}






