
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
setwd("//wsl.localhost/Ubuntu/root/qiime2-Macalady9")
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

SampleSheet<- read_excel(TemplateXLSX, sheet="SampleSheet.csv", skip=19) %>% mutate(Sample_ID=gsub("-","_", Sample_ID))

fastqs<-
  list.files(ReadDir, recursive = TRUE, full.names = TRUE, pattern="\\.fastq\\.gz") %>%
  tibble(File=.) %>%
  filter(!grepl("Undetermined_", File)) %>% # these are reads that could not be demultiplexed
  mutate(Sample_ID=basename(File) %>% gsub(Suffix, "", .) %>% gsub("\\.fastq\\.gz","", .) %>% gsub("-","_", .)) %>% #note: need to switch _ to - in file names
  mutate(Direction=case_when(
    grepl("_R1_", .$File) ~ "Forward_Read",
    grepl("_R2_", .$File) ~ "Reverse_Read"
  )) %>%
  spread(key=Direction, value=File)

SampleSheet<-SampleSheet %>% left_join(fastqs) %>% filter(grepl(Project, Sample_ID))

interactive_table(SampleSheet)

SampleSheet %>% 
  filter(is.na(Forward_Read)) %>%
  interactive_table()

SampleSheet<-SampleSheet %>% filter(!is.na(Forward_Read)) %>% mutate(Nreads_illumina=countFastq(Forward_Read)$records)

write_csv(SampleSheet, "Logs/Nreads_raw.txt")

fastqs %>%
  select(`sample-id`=1, forward=2, reverse=3) %>%
  gather(-`sample-id`, value=`absolute-filepath`, key=`direction`) %>%
  select(1,3,2) %>%
  arrange(`sample-id`) %>%
  filter(`sample-id`!="Undetermined") %>%
  filter(`sample-id` %in% SampleSheet$Sample_ID) %>%
  write_csv("Intermediates/manifest.csv")
    #can altered this last line so it doesnt destroy the manifest which i have to edit by hand to add absolute paths to

SampleSheet %>%
  sample_n(3) %>%
  select(Forward_Read, Reverse_Read) %>%
  gather() %>%
  arrange(value) %>%
  pull(value) %>%
  plotQualityProfile() +
  facet_wrap(~file, ncol = 2)
ggsave("Figures/ReadQuality.pdf", device="pdf", height=8.5, width=11)









#in qiime:

cd qiime2-Macalady9

if [ ! -f $PWD/Intermediates/Reads.qza ]; then
echo $(date) Generating Read Artifact

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path $PWD/Intermediates/manifest.csv \
--output-path $PWD/Intermediates/Reads.qza \
--input-format PairedEndFastqManifestPhred33

else
  echo $(date) Skipping making Read Artifact as already complete
fi







