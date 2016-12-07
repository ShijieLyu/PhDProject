###################################################
#
# 	Script for Erhalten von sequenzen aus mm10

# in R 2. version alles am 30.4.2014 aktualisiert, dauert lange: 1h
#      S.Kärst
###################################################

#    source("http://bioconductor.org/biocLite.R")
#    biocLite("BSgenome")
#    biocLite("BSgenome.Mmusculus.UCSC.mm10")
#   biocLite("Mus.musculus")


library(BSgenome.Mmusculus.UCSC.mm10)
#install.packages("biomaRt")
library(biomaRt)
head(listMarts())
ensembl <- useMart("ensembl")
ensembl
listDatasets(ensembl)
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
ensembl
head(listFilters(ensembl))
head(listAttributes(ensembl))
source("http://bioconductor.org/biocLite.R")
biocLite("lumiMouseAll.db")

### ACTION
toString(subseq(Mmusculus$chr3, 36599328,36601675))




##############
suppressMessages(library(BSgenome.Mmusculus.UCSC.mm10))
#untereinheit extraieren (zum Bsp)
subseq(Mmusculus$chr3, 18357968,18360987)
#erhalten der Sequenz
toString(subseq(Mmusculus$chr15, 18357968,18358000))

# primer für gap test in Bbs7
forw <- toString(subseq(Mmusculus$chr3, 36599328,36599928))
# TCATTACTTACCTGAATATGTGGCCAGCACAATTTCATCATAGCCATCTTTTCCCACACAGCCACCCTGGATGGAGGTGACGCTTTCAGACAGCATCTAGAGTAAATGTAAAGTCCACTTTGCTAGCATTGAAACATGTCAAAGACATGTTAAAACAAGCCACAAAACTTTGCTTCTTGTCTCCAGATAGTCAAAACATATTTAGTGTATGTTCATCATGTGCAGACACTGAGGAAAACCTGAATGAATGCTAGAAAATTCCAGCATTTTTAAATATAAGAAAGAAGATTCAAAATAGCTAGTATATTTGTGACAAAACTACATAAAGAAGAGAGGTTTGGTAAATTATATAAACAAAACAGCTAGTAAATTTAGTCTGTAAATCTCTTCAGAGGATTCTTTTTCTCTAATGAAAAAAAAACTATGCTTTCACTGCATAATCTTAAGTTAACGTCTAGAACTCTTATCATAGTCTCATGGAAGTTGTACATAATATTTTAAAGCCATGGCATTTAAAAATCATACCCAGAGCTGGAGAGATGGTTCTGCGGTTAAGAGCGCCCGCTGCTTCTCCAGAGGTCCTGAGTTCAATCCCCAGCAA
rev1 <- toString(subseq(Mmusculus$chr3, 36600040,36600640))
# TGGTGGCACATGCCTTTAATCCCAGCACTCGGGAGGCAGAGGCAGGCAGATTTCTGAGTTCGAGGTCAGCCTGGTCTACAAAGTGAGTGCCAGGACAGTCAGGGCTACACAGAGAAACCCTGTCTCGAAAAACCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAAGAAAGAAAGAAAGAAATCTTAAAAAGAAACATACCCATTGACTCTAAATACCAGGAAAAATTTATATACAAAATCATTGACATTGAAAGTAAGAAACATCAGTTAGAATAACTGAATGCAAGAGAAATAAATGCTACATTGTTCTTTTTCCTCCTTGAACTTTAATTTTCATTTGCAGTGTGTGGGCATGCACACTCCAGTGTGAGTGTGGAGTCAGACAACCTTGAGGAGCTGTCCCTCCTGTTCCACCTTACATGGCTTTGGGGACTGAAAGCAGTCTGCCAGGCTTGTGTGGCACCTCTACACACTAAACCATCTGGCCAGATCTACTGCTTTTAAAAGTCAACATTACTTCTCTCTCTCTCTCTCTCTCTCTCTCATTTGTGCATGTGTACACATACATAGTGTGTATCAACTTTTATTTTCATTGCA
rev2 <- toString(subseq(Mmusculus$chr3, 36601075,36601675))
# GATCTCATTCATATGCAGAACGCCCTCACCTCACGGGGGCTAAGGACAGAATGGCGGGCGCTGCAGTTGTTGAGGAGTAGGGAAGGGATGAGTACGGAAAAGCTGGCTAAGGAGCATCAACATAGTTGCGTATATATGAAGTTCAGACCGTTGCTTTAGTTTTCCTTGAAATTCTGATGGTGACTCATGGCACAGAAAGGAGGAGTGCTTGCCACTGTATAGGGCGACCACAGGTGACAATTATCTAGCACAGGTTTCAAAGGGCTACAAGAAACCATCTTGAGCATTTTCACCCAAGAACTATGTCTACCCTATGTAAACAGTATCCTCAGGGACCACAGGATAACTCCCTGAGTACTGATGCTTGCAGCCAAACTGAGACCTGAGTTTGATTCCTGGGACCCACATGGTGAAACAGGAGAACCAACACCCGGAAGTTATCCTCTGACCTCCGCATGTGTGCCGTGGCGTGCACTTACATACAACATACACACAAATACCACCCCACAAACATTAATATATAAATTTATAAATAAACAATTTAAAACATCTAAAGGTAAAATATTTGGTACCATTATATACTACCATTCATATATAGGTT
toString(subseq(Mmusculus$chr3, 36599328,36601675))

# primer für Ndufs1:171240055

# sequencing für Dmd deleterious
toString(subseq(Mmusculus$chrX, 83728876,83729876))

# 184	184	62	R/C	Cgc/Tgc	-	MT_2934_C/T	YES	deleterious(0)	mt-Nd1
# MT:2934
toString(subseq(Mmusculus$chrM, 2434,3434))
#MT:13781	G	ENSMUSG00000064368	ENSMUST00000082419	Transcript	missense_variant	290	290	97	I/T	aTt/aCt	-	MT_13781_A/G	YES	deleterious(0.04)	mt-Nd6	MT	13552	14070	-1	MGI:102495	mt-Nd6
toString(subseq(Mmusculus$chrM, 13281,14281))
# X:68678857	A	ENSMUSG00000000838	ENSMUST00000114657	Transcript	missense_variant	47	41	14	A/D	gCt/gAt	-	X_68678857_C/A	-	deleterious(0)	Fmr1
toString(subseq(Mmusculus$chrX, 68678357,68679357))
# X:106118378	A	ENSMUSG00000033792	ENSMUST00000055941	Transcript	missense_variant	3837	3761	1254	R/Q	cGg/cAg	-	X_106118378_G/A	YES	deleterious(0.03)	Atp7a
toString(subseq(Mmusculus$chrX, 106117878,106118878))


table(forw)

36885494,36885514
36889868,36889888
36890825,36890845
36895367,36895387
36921287,36921307
36921290,36921310
36940809,36940829
36948332,36948352
36948357,36948377
36948374,36948394
36948409,36948429
37041878,37041898

#geht nicht
start <-c(36885494,36889868,36890825,36895367,36921287,36921290,36940809,36948332,36948357,36948374,36948409,37041878)
end   <-c(36885514,36889888,36890845,36895387,36921307,36921310,36940829,36948352,36948377,36948394,36948429,37041898)
toString(subseq(Mmusculus$chr3,start,end))

#derreihenach
drn <- 1:12

#einzeln
toString(subseq(Mmusculus$chr3,start[1],end[1]))



toString(subseq(Mmusculus$chr3,start[drn],end[drn]))

#[1] "GATAAAACCCGGGAAAATGGG"

toString(subseq(Mmusculus$chr3,start[2],end[2])) # TGATGTTGATCTCTACTATTA
toString(subseq(Mmusculus$chr3,start[3],end[3])) # CAGAAGGGGACATAAGCAGTG
toString(subseq(Mmusculus$chr3,start[4],end[4])) # GAATATTATTACAGATGCTAC
toString(subseq(Mmusculus$chr3,start[5],end[5])) # ...
toString(subseq(Mmusculus$chr3,start[6],end[6]))
toString(subseq(Mmusculus$chr3,start[7],end[7]))
toString(subseq(Mmusculus$chr3,start[8],end[8]))
toString(subseq(Mmusculus$chr3,start[9],end[9]))
toString(subseq(Mmusculus$chr3,start[10],end[10]))
toString(subseq(Mmusculus$chr3,start[11],end[11]))
toString(subseq(Mmusculus$chr3,start[12],end[12]))


GATAAAACCCgGGAAAATGGG
TGATGTTGATcTCTACTATTA
CAGAAGGGGAcATAAGCAGTG
GAATATTATTaCAGATGCTAC...
TCTTCTCCATcTATGTCACCA
TCTCCATCTAtGTCACCACGC
TCGTTTAGCAaTAGATGGATC
CATAGAACATgTGTACATTGT
GGTGAGAAAAgAGGAGAGAGT
GAGTGAGCAGgTCACGTCTCA
GTTGCTATCAcACTTACCTTA
CAAGACAATGaCCAACAACCT


#### Transcription factor tests http://www.bioconductor.org/help/workflows/generegulation/
biocLite(c("MotifDb", "GenomicFeatures", "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene", 
    "org.Sc.sgd.db", "BSgenome.Scerevisiae.UCSC.sacCer3", "motifStack", "seqLogo"))

## See system.file("LICENSE", package="MotifDb") for use restrictions.

biocLite(c("org.Sc.sgd.db","GenomicFeatures"))
library(MotifDb)
library(seqLogo)
library(motifStack)
library(Biostrings)
library(GenomicFeatures)
library(org.Sc.sgd.db)
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)

query(MotifDb, "DAL80")

## MotifDb object of length 2
## | Created from downloaded public sources: 2012-Nov-01
## | 2 position frequency matrices from 2 sources:
## |        JASPAR_CORE:    1
## |             ScerTF:    1
## | 1 organism/s
## |        Scerevisiae:    2
## Scerevisiae-JASPAR_CORE-DAL80-MA0289.1 
## Scerevisiae-ScerTF-DAL80-harbison

pfm.dal80.jaspar <- query(MotifDb, "DAL80")[[1]]
seqLogo(pfm.dal80.jaspar)

dal1 <- "YIR027C"
chromosomal.loc <- transcriptsBy(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene, by = "gene")[dal1]
promoter.dal1 <- getPromoterSeq(chromosomal.loc, Scerevisiae, upstream = 1000, 
    downstream = 0)
pcm.dal80.jaspar <- round(100 * pfm.dal80.jaspar)
matchPWM(pcm.dal80.jaspar, unlist(promoter.dal1)[[1]], "90%")

##   Views on a 1000-letter DNAString subject
## subject: TTGAGGAGTTGTCCACATACACATTAGTGTTG...AAAAAAAAGTGAAATACTGCGAAGAACAAAG
## views:
##     start end width
## [1]   620 626     7 [TGATAAG]
## [2]   637 643     7 [CGATAAG]

query(MotifDb, "DAL80")

## MotifDb object of length 2
## | Created from downloaded public sources: 2012-Nov-01
## | 2 position frequency matrices from 2 sources:
## |        JASPAR_CORE:    1
## |             ScerTF:    1
## | 1 organism/s
## |        Scerevisiae:    2
## Scerevisiae-JASPAR_CORE-DAL80-MA0289.1 
## Scerevisiae-ScerTF-DAL80-harbison

dal80.jaspar <- query(MotifDb, "DAL80")[[1]]
dal80.scertf <- query(MotifDb, "DAL80")[[2]]
seqLogo(dal80.jaspar)

seqLogo(dal80.scertf)

pfm.dal80.jaspar <- new("pfm", mat = query(MotifDb, "dal80")[[1]], name = "DAL80-JASPAR")
pfm.dal80.scertf <- new("pfm", mat = query(MotifDb, "dal80")[[2]], name = "DAL80-ScerTF")
plotMotifLogoStack(DNAmotifAlignment(c(pfm.dal80.scertf, pfm.dal80.jaspar)))

query(MotifDb, "gat1")

## MotifDb object of length 3
## | Created from downloaded public sources: 2012-Nov-01
## | 3 position frequency matrices from 3 sources:
## |        JASPAR_CORE:    1
## |             ScerTF:    1
## |           UniPROBE:    1
## | 1 organism/s
## |        Scerevisiae:    3
## Scerevisiae-JASPAR_CORE-GAT1-MA0300.1 
## Scerevisiae-ScerTF-GAT1-zhu 
## Scerevisiae-UniPROBE-Gat1.UP00287

pfm.gat1.jaspar = new("pfm", mat = query(MotifDb, "gat1")[[1]], name = "GAT1-JASPAR")
pfm.gat1.scertf = new("pfm", mat = query(MotifDb, "gat1")[[2]], name = "GAT1-ScerTF")
pfm.gat1.uniprobe = new("pfm", mat = query(MotifDb, "gat1")[[3]], name = "GAT1-UniPROBE")
plotMotifLogoStack(c(pfm.gat1.uniprobe, pfm.gat1.scertf, pfm.gat1.jaspar))

#...
### TRASH

## 
library(AnnotationDbi)
#help(package="AnnotationDbi")

require(lumiMouseAll.db)
## display the cols
cols(lumiMouseAll.db)
## get the 1st 6 possible keys
keys <- keys(lumiMouseAll.db) 
keys
## lookup gene symbol and unigene ID for the 1st 6 keys
select(lumiMouseAll.db, keys=keys, cols = c("SYMBOL","UNIGENE"))

## get keys based on unigene
keyunis <- head( keys(lumiMouseAll.db, keytype="UNIGENE") )
keyunis
## list supported key types
keytypes(lumiMouseAll.db)
## lookup gene symbol and unigene ID based on unigene IDs by setting
## the keytype to "UNIGENE" and passing in unigene keys:
select(lumiMouseAll.db, keys=keyunis, cols = c("ONTOLOGY","PMID","UNIGENE","ENTREZID","PFAM" ,"IPI"),
       keytype="UNIGENE")



