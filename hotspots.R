#### Possibly you need to install UniProt.ws and stringr first ####
# source("https://bioconductor.org/biocLite.R")
# biocLite("UniProt.ws")
# install.packages("stringr")

library(UniProt.ws)
library(stringr)
library(ggplot2)

options(stringsAsFactors = F)

### HERE: Specify parameters for search
hotspot_min_length = 5 ### minimum length of a hotspot in AS
n_hit_wonder = 1 ### number of hits two be classified as n_hit_wonder
hotspot_min_number_of_patients = 5 ### minimum number of patients mapping peptides onto the hotspot
hotspot_min_number_of_as = 1 ### minimum number of aminoacids to map peptides onto the hotspot ( maximum = hotspot_min_length - 1)


### input CSV with columns "Sequence", "Dignity" and "Single.Proteins" (other columns will be ignored)
df<-read.csv("Leon_Projects/Annika/CML_haema_classII_single_proteins_forCluster.csv")
outdir <- paste0("hotspot_analysis_") #,date())
dir.create(outdir)

### load uniprot human taxonomy
up <- UniProt.ws(taxId=9606)

### loop over all proteins in input df
for (protid in unique(df$Single.Proteins)){
  print(protid)
  seq<-tryCatch({select(up, keys=c(protid), columns=c("SEQUENCE"), keytype="UNIPROTKB")$SEQUENCE}, error=function(err){seq<-"XXX"})
  df_sub<-df[which(df$Single.Proteins %in% c(protid)),][c("Sequence","Dignity","Single.Proteins","Patient..Donor")]
  ### retrieve sequence of proteinid and match peptide sequence start and end point
  starts<-rep(NA,length(df_sub[,1]))
  ends<-rep(NA,length(df_sub[,1]))
  sums_benign<-rep(0,str_length(seq))
  sums_malign<-rep(0,str_length(seq))
  for (row in seq(1,length(df_sub[,1]))){
    peptide<-df_sub$Sequence[row]
    match<-str_locate_all(pattern=peptide, seq)[[1]]
    if (length(match)>0){
      start<-as.numeric(match[,1])
      end<-as.numeric(match[,2])
      ### collect coverage of protein sequence by summing over peptides
      if(df_sub$Dignity[row]=="benign"){                                                                       ######## SPELLING SENSITIVE ##########
        sums_benign[start:end]<-sums_benign[start:end]+1    
      } else {
        sums_malign[start:end]<-sums_malign[start:end]+1     
      }
      starts[row]<-start
      ends[row]<-end
    } else {
      starts[row]<- (-1)
      ends[row]<- (-1)
    }
  }

  ### test if there are malignant exclusive hotspots greater length 7
  exclusive_switch <- FALSE
  one_hit_switch <- FALSE
  sums_malign_norm<-sums_malign/sums_malign
  sums_benign_norm<-sums_benign/sums_benign
  sums_malign_norm[which(is.na(sums_malign_norm))]<-0
  sums_benign_norm[which(is.na(sums_benign_norm))]<-0
  sums_difference <- sums_malign_norm-sums_benign_norm
  sequential_intervals <-rle(diff(sort(which(sums_difference>0))))
  hotspot <- list()
  hotspot_pep <- list()
  if (any(sequential_intervals$lengths>=hotspot_min_length & sequential_intervals$values==1)) {   ### HERE change length 7
    ### map hotspot back to sequence ###
    sums_diff_not<-which(sums_difference < 1)
    seq_splt=strsplit(seq,"")[[1]]
    seq_splt[sums_diff_not]<-rep("X",length(sums_diff_not))
    seq_splt<-paste(seq_splt, collapse='')
    seq_split<-strsplit(seq_splt, 'X')[[1]]
    seq_split<-seq_split[which(seq_split!="")]
    seq_split<-seq_split[which(str_length(seq_split)>=hotspot_min_length)]                        ### here change length 7 as well
    print(seq_split)
    ### count number of patients with epitopes matching the hotspot
    for (seq_h in seq_split){
        match<-str_locate_all(pattern=seq_h, seq)[[1]]
        start<-as.numeric(match[,1]) + hotspot_min_number_of_as - 1
        end<-as.numeric(match[,2]) - hotspot_min_number_of_as + 1
        for (row in seq(1,length(starts))){
          if ((starts[row] < end) & (ends[row] > start)) {
            hotspot[[seq_h]]<-append(hotspot[[seq_h]], df_sub$Patient..Donor[row])
            hotspot[[seq_h]]<-unique(hotspot[[seq_h]])
            hotspot_pep[[seq_h]]<-append(hotspot_pep[[seq_h]], df_sub$Sequence[row])
            hotspot_pep[[seq_h]]<-unique(hotspot_pep[[seq_h]])
          }          
        }
    }
    print(sapply(hotspot, length))
    ### check if at least n patients match peptides in any hotspot and at least 2 unique peptide ids for this protein
    if (any(sapply(hotspot, length)>=hotspot_min_number_of_patients)) {      #### HERE: adjust threshold for number of patients
      exclusive_switch <- TRUE
    }
  }
  if (length(unique(df_sub$Sequence))== n_hit_wonder){     #### HERE: adjust threshold for n hit wonders
    one_hit_switch <- TRUE
  } 
    
  ### count benign and malignant samples
  peptide_counts<-table(df_sub$Dignity)
  if ("malignant" %in% names(peptide_counts)){
    malign_count<-as.numeric(peptide_counts[which(names(peptide_counts)=="malignant")])                     ######## SPELLING SENSITIVE ##########
    malign_c<-malign_count
  } else {
    ### set to 1 if none to avoid division by 0
    malign_count<-1
    malign_c<-malign_count-1
  }
  if ("benign" %in% names(peptide_counts)){
    benign_count<-as.numeric(peptide_counts[which(names(peptide_counts)=="benign")])                        ######## SPELLING SENSITIVE ##########
    benign_c<-benign_count
  } else {
    benign_count<-1
    benign_c<-benign_count-1
  }

  ### summarize sequence coverage counts, normalize over counts and multiply benign by -1 to visualize horizontal mirror depiction
  sm<-cbind(as.data.frame(sums_malign/malign_count),rep(paste0("malignant n=",malign_c),str_length(seq)),seq(1,str_length(seq)))
  sb<-cbind(as.data.frame(sums_benign*(-1)/benign_count),rep(paste0("benign n=",benign_c),str_length(seq)),seq(1,str_length(seq)))
  colnames(sm)<-c("Count",protid,"Sequence")
  colnames(sb)<-c("Count",protid,"Sequence")
  sums<-rbind(sm,sb)
  sums[[2]]<-factor(sums[[2]], levels = rev(levels(factor(sums[[2]]))))

  ### plot
  p<-ggplot(sums, aes_string(x="Sequence", y="Count", fill=protid)) + 
    geom_bar(stat="identity", position="identity", width = 1) +
    geom_hline(yintercept = 0) +
    theme_classic()
  
  ### output start and end mapping per peptide
  dir.create(paste0(outdir,"/only_benign"))
  dir.create(paste0(outdir,"/only_malign"))
  dir.create(paste0(outdir,"/tumor_associated"))
  dir.create(paste0(outdir,"/all"))
  dir.create(paste0(outdir,"/n_hit_wonders"))
  df_out<-cbind(df_sub,starts,ends)
  df_out_hot<-df_out[which(df_out$Sequence %in% unlist(hotspot_pep)),]
  if (exclusive_switch){
    capture.output(print(hotspot_pep), file = paste0(outdir,"/tumor_associated/",protid,"_hotspots_matches.csv"))
    write.csv(df_out_hot, file = paste0(outdir,"/tumor_associated/",protid,"_hotspot_peptides_only.csv"))
    write.csv(df_out, file = paste0(outdir,"/tumor_associated/",protid,"_all_peptides.csv"))
    ggsave(plot = p, filename = paste0(outdir,"/tumor_associated/",protid,"_hotspots.png"))
  } else if (one_hit_switch){
    write.csv(df_out, file = paste0(outdir,"/n_hit_wonders/",protid,"_all_peptides.csv"))
    ggsave(plot = p, filename = paste0(outdir,"/n_hit_wonders/",protid,"_hotspots.png"))    
  } else if (benign_c==0){
    write.csv(df_out, file = paste0(outdir,"/only_malign/",protid,"_all_peptides.csv"))
    ggsave(plot = p, filename = paste0(outdir,"/only_malign/",protid,"_hotspots.png")) 
  } else if (malign_c==0){
    write.csv(df_out, file = paste0(outdir,"/only_benign/",protid,"_all_peptides.csv"))
    ggsave(plot = p, filename = paste0(outdir,"/only_benign/",protid,"_hotspots.png")) 
  } else {
    write.csv(df_out, file = paste0(outdir,"/all/",protid,"_peptides.csv") )
    ggsave(plot = p, filename = paste0(outdir,"/all/",protid,"_hotspots.png"))
  }
}

