# Create chromosome name mappings

The goal is to create a frame which maps the various genome assembly 
chromosome names to one another, for instance, from UCSC to ensembl. I prefer 
to do this in R.

## load dependencies

```{r setup}
library(rtracklayer)
library(tidyverse)
```

## save paths and the regex used to clean the seqnames from the fasta

These were downloaded from their respective sites on 20220812

```{r}
genome_paths = list(
  ensembl = "~/Downloads/homo_sapien/ensembl_genome.fa",
  ucsc = "~/Downloads/homo_sapien/ucsc_genome.fa",
  refseq = "~/Downloads/homo_sapien/refseq_genome.fa",
  genbank = "~/Downloads/homo_sapien/genbank_genome.fa"
)

These regex patterns were determined by looking at the names in the fasta file

extract_seqnames_regex = list(
  refseq = "GRCh38.p14 Primary Assembly$|mitochondrion, complete genome$",
  genbank = "[\\d,X,Y], GRCh38 reference primary assembly$|mitochondrion, complete genome$"
)
```

## create a list of the primary 22 autosomes, XY and mitochondria for each source
```{r}
genomes = list(
  ensembl = unlist(map(str_split(
    names(rtracklayer::import(genome_paths$ensembl, 
                              format = "fasta"))[1:25], " "), ~.[[1]][[1]])),
  ucsc = names(rtracklayer::import(genome_paths$ucsc, format = "fasta"))[1:25],
  refseq = 
    names(rtracklayer::import(genome_paths$refseq, format = "fasta")) %>% 
      .[str_detect(., extract_seqnames_regex$refseq)],
  genbank = 
    names(rtracklayer::import(genome_paths$genbank, format = "fasta")) %>% 
      .[str_detect(., extract_seqnames_regex$genbank)]
)
```

## Create a frame which allows mapping between sources

Note that gencode is the same as UCSC

```{r}
clean_seqnames_regex = list(
  refseq = "NC_\\d+.\\d+ Homo sapiens chromosome |NC_\\d+.\\d+ Homo sapiens|, GRCh38.p14 Primary Assembly|, complete genome",
  genbank = "CM\\d+.\\d+ Homo sapiens chromosome |J\\d+.\\d+ Homo sapiens|, GRCh38 reference primary assembly|, complete genome"
)
genomes_df_list = list(
  ensembl = tibble(
    ensembl = genomes$ensembl, 
    seq = genomes$ensembl),
  ucsc    = tibble(
    ucsc = genomes$ucsc, 
    seq = str_remove(genomes$ucsc, "chr")),
  # gencode is the same as ucsc
  gencode = tibble(
    gencode = genomes$ucsc, 
    seq = str_remove(genomes$ucsc, "chr")),
  refseq  = tibble(
    refseq = unlist(map(str_split(genomes$refseq, " "), ~.[[1]][[1]])), 
    seq = str_remove_all(genomes$refseq, clean_seqnames_regex$refseq)),
  genbank  = tibble(
    genbank = unlist(map(str_split(genomes$genbank, " "), ~.[[1]][[1]])), 
    seq = str_remove_all(genomes$genbank, clean_seqnames_regex$genbank))
)

cleanup_genomes_df = function(df){
  df %>%
    mutate(seq = ifelse(str_detect(seq, "mitochondrion|MT"), "M", seq))
}

genomes_df_list = map(genomes_df_list, cleanup_genomes_df)

chr_map = genomes_df_list %>%
  reduce(left_join) %>%
  select(-seq)

write_csv(chr_map, "~/code/nf-core-callingcards/assets/human/chr_map.csv")
```