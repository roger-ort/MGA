
#conda update conda
#conda activate qiime2-amplicon-2024.10

library(rmarkdown)
library(tidyverse)
library(readxl)
library(dada2)
library(plotly)
library(ggtree)
library(qiime2R)
library(ShortRead)
sessionInfo()

#Raws should be in number folder
setwd("//wsl.localhost/Ubuntu/root/qiime2-MacaladyControls")
getwd()


for(l in readLines("settings.txt")){
  l<-str_split(gsub("export ", "", l), "=")[[1]]
  assign( l[1], l[2])
}
rm(l)

dir.create("Intermediates", showWarnings=FALSE)
dir.create("Output", showWarnings=FALSE)
dir.create("Figures", showWarnings=FALSE)
dir.create("Logs", showWarnings=FALSE)

SampleSheet <- read_excel(TemplateXLSX,sheet="SampleSheet.csv",skip=19) %>% mutate(Sample_ID=gsub("-","_",Sample_ID))
fastqs <- list.files(ReadDir,recursive=TRUE,full.names=TRUE,pattern="\\.fastq\\.gz$") %>% tibble(File=.) %>% filter(!grepl("Undetermined_",File)) %>% mutate(File=gsub("^\\\\wsl\\.localhost\\\\Ubuntu","",File),File=gsub("\\\\","/",File),Sample_ID=basename(File) %>% gsub("\\.fastq\\.gz$","",.) %>% gsub("_R[12]_001$","",.),Direction=case_when(grepl("_R1_001\\.fastq\\.gz$",File)~"forward",grepl("_R2_001\\.fastq\\.gz$",File)~"reverse")) %>% tidyr::pivot_wider(names_from=Direction,values_from=File)
fastqs <- fastqs %>% mutate(across(c(forward,reverse),~sub("^//wsl\\.localhost/Ubuntu","",.x)))
SampleSheet <- SampleSheet %>% dplyr::select(-dplyr::any_of(c("Forward_Read","Reverse_Read"))) %>% dplyr::left_join(fastqs %>% dplyr::select(Sample_ID,Forward_Read=forward,Reverse_Read=reverse),by="Sample_ID")
write_csv(SampleSheet,"Logs/Nreads_raw.txt")
fastqs %>% dplyr::select(`sample-id`=Sample_ID,forward,reverse) %>% tidyr::pivot_longer(c(forward,reverse),names_to="direction",values_to="filepath") %>% dplyr::rename(`absolute-filepath`=filepath) %>% dplyr::select(`sample-id`,`absolute-filepath`,direction) %>% dplyr::arrange(`sample-id`) %>% readr::write_csv("Intermediates/manifest.csv")


SampleSheet %>%
  sample_n(3) %>%
  dplyr::select(Forward_Read, Reverse_Read) %>%
  gather() %>%
  arrange(value) %>%
  pull(value) %>%
  plotQualityProfile() +
  facet_wrap(~file, ncol = 2)
ggsave("Figures/ReadQuality.pdf", device="pdf", height=8.5, width=11)






#Post qiimerun


read_qza("Output/ASV_sequences.qza")$data %>%
  names() %>%
  tibble(`feature-id`=.) %>%
  write_tsv("Intermediates/SVs_passing_size.txt")

read_qza("Output/ASV_sequences.qza")$data %>% sapply(., length) %>% tibble(ASV_Length=.) %>% mutate(Step="Post_filter") %>%
  bind_rows(
    read_qza("Intermediates/SVsequences.qza")$data %>% sapply(., length) %>% tibble(ASV_Length=.) %>% mutate(Step="Pre_filter")
  ) %>%
  ggplot(aes(x=ASV_Length, color=Step)) + geom_freqpoly() + theme_q2r() + ylab("N features")

ggsave("Figures/ASV_lengths.pdf", height=3, width=4, useDingbats=F)

read_qza("Output/ASV_sequences.qza")$data %>% sapply(., length) %>% tibble(ASV_Length=.) %>% group_by(ASV_Length) %>% summarize(N=n()) %>% knitr::kable()

print(paste0(
  "Size selection removed ",
  (sum(read_qza("Intermediates/SVtable.qza")$data)-sum(read_qza("Output/ASV_table.qza")$data))/sum(read_qza("Intermediates/SVtable.qza")$data)*100,
  "% of reads"
))

ReadTracking<-
  SampleSheet %>%
  dplyr::select(Sample_ID, Sample_Plate, Sample_Well)

ReadTracking<-
  ReadTracking %>%
  left_join(
    read_qza("Logs/Dada_stats.qza")$data %>%
      rownames_to_column("Sample_ID") %>%
      dplyr::select(Sample_ID, Primer_Trimming=input, Dada2_Filtering=filtered, Dada2_Denoising=denoised, Dada2_Overlap=merged, Dada2_ChimeraRemoval=non.chimeric)
  )

ReadTracking<-
  ReadTracking %>%
  left_join(
    read_qza("Output/ASV_table.qza")$data %>%
      colSums() %>% data.frame(Size_Selection=.) %>% rownames_to_column("Sample_ID")
  )

interactive_table(ReadTracking)




rplot<-
  ReadTracking %>%
  pivot_longer(!c(Sample_ID, Sample_Well, Sample_Plate),  names_to="Step", values_to ="Nreads") %>%
  mutate(Step=factor(Step, levels=colnames(ReadTracking)[-(1:3)])) %>%
  ggplot(aes(x=Step, y=Nreads, group=Sample_ID, label=Sample_ID)) +
  geom_line() +
  theme_q2r() +
  ylab("Number of Reads") +
  xlab("Pipeline Step") +
  theme(axis.text.x = element_text(angle=45, hjust=1))

ggplotly(rplot)

ggsave("Figures/Read_Tracking.pdf", rplot, device="pdf", height=6, width=4, useDingbats=F)
rm(rplot)






if(!file.exists("qiime merge/ASV_d2taxonomy.txt")){
  message(date(), " Assigning taxonomy with dada2")
  
  #seqs<-read_qza("qiime merge/MergeSeqs.qza")$data
  taxonomy <- assignTaxonomy(as.character(ASV_seq), "silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=as.integer(8))
  taxonomy <- addSpecies(taxonomy, "silva_species_assignment_v138.1.fa.gz", allowMultiple=TRUE)
  colnames(taxonomy)[8]<-"Species_ambiguous"
  
  as.data.frame(taxonomy) %>%
    rownames_to_column("Sequence") %>%
    left_join(tibble(ASV=names(ASV_seq), Sequence=as.character(ASV_seq))) %>%
    dplyr::select(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species, Species_ambiguous, Sequence) %>%
    write_tsv("qiime merge/ASV_d2taxonomy.txt")
  
} else {
  message(date(), " Skipping Dada2 taxonomy as already complete")
}

taxonomy<-read_tsv("qiime merge/ASV_d2taxonomy.txt")






