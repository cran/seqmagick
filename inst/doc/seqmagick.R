## ----style, echo=FALSE, results="asis", message=FALSE-------------------------
knitr::opts_chunk$set(tidy = FALSE,
		   message = FALSE)

## ----echo=FALSE, results="hide", message=FALSE--------------------------------
library(magrittr)
library(Biostrings)
library(muscle)
library("seqmagick")

## -----------------------------------------------------------------------------
tmpgb <- tempfile(fileext = '.gb')
tmpfa <- tempfile(fileext = '.fa')
download_genbank(acc='AB115403', format='genbank', outfile=tmpgb)
download_genbank(acc='AB115403', format='fasta', outfile=tmpfa)
readLines(tmpgb)[1:10]
readLines(tmpfa)

## -----------------------------------------------------------------------------
fa_file <- system.file("extdata/HA.fas", package="seqmagick")
## use the small subset to save compilation time of the vignette
fa2 <- tempfile(fileext = '.fa')
fa_read(fa_file) %>% bs_filter('ATGAAAGTAAAA', by='sequence') %>% fa_write(fa2, type='interleaved')


alnfas <- tempfile(fileext = ".fas")
fa_read(fa2) %>% bs_aln(quiet=TRUE) %>% fa_write(alnfas)

## phylip format is only for aligned sequences
tmpphy <- tempfile(fileext = ".phy")
fas2phy(alnfas, tmpphy, type = 'sequential')
readLines(tmpphy)[1:3]

## -----------------------------------------------------------------------------
phy2fas(tmpphy, alnfas, type = 'interleaved')
readLines(alnfas)[1:4]

## -----------------------------------------------------------------------------
readLines(fa2)[1:10]

## -----------------------------------------------------------------------------
tmpfas <- tempfile(fileext='.fa')
fa_read(fa2) %>% fa_write(tmpfas, type="sequential")
readLines(tmpfas)[1:4]

## -----------------------------------------------------------------------------
tmpphy2 <- tempfile(fileext = '.phy')
phy_read(tmpphy) %>% phy_write(tmpphy2, type="interleaved")
readLines(tmpphy2)[1:4]

## -----------------------------------------------------------------------------
bs <- fa_read(fa_file)
bs_filter(bs, 'ATGAAAGTAAAA', by='sequence')

aln <- bs_filter(bs, 'ATGAAAGTAAAA', by='sequence') %>% bs_aln(quiet=TRUE)

aln

bs_consensus(aln)

## ----echo=FALSE---------------------------------------------------------------
sessionInfo()

