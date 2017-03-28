## ---- warning=FALSE------------------------------------------------------
library('TreeExp')

## ---- warning=FALSE------------------------------------------------------

readsCount.table = read.table(system.file('extdata/tetraexp.read.counts.raw.txt', 
                                            package='TreeExp'), header = T)
head(readsCount.table[,1:10])

geneInfo.table = read.table(system.file('extdata/tetraexp.length.ortholog.txt',
                                        package='TreeExp'), header = T)
head(geneInfo.table)


## ---- eval=FALSE---------------------------------------------------------
#  taxa.objects = TEconstruct(readCountsFP = system.file('extdata/tetraexp.read.counts.raw.txt', package='TreeExp'),
#    geneInfoFP = system.file('extdata/tetraexp.length.ortholog.txt', package='TreeExp'),
#    taxa = "all", subtaxa = c("Brain", "Cerebellum"), normalize = "TPM")

## ------------------------------------------------------------------------
data(tetraexp)

## ------------------------------------------------------------------------
print(tetraexp.objects, details = TRUE)

## ------------------------------------------------------------------------
print(tetraexp.objects[[1]], printlen = 6)
head(tetraexp.objects[[1]]$normExp.val)

