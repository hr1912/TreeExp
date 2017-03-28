## ---- warning=FALSE------------------------------------------------------
library('TreeExp')

## ------------------------------------------------------------------------
data(tetraexp)

## ---- message=FALSE------------------------------------------------------
dismat <- expdist(tetraexp.objects, taxa = "all",
                 subtaxa = "Brain",
                 method = "pea")
as.dist(dismat)

## ---- message=FALSE, warning=FALSE, results='hide'-----------------------

expression_table <- exptabTE(tetraexp.objects, taxa = "all",
                            subtaxa = "Brain")

dismat <- dist.pea(expression_table)
colnames(dismat) <- colnames(expression_table)
rownames(dismat) <- colnames(dismat)


## ---- eval=FALSE---------------------------------------------------------
#  
#  dismat <- dist.pea(your_own_dataframe)
#  colnames(dismat) <- colnames(your_own_dataframe)
#  rownames(dismat) <- colnames(dismat)
#  

## ---- message=FALSE, warning=FALSE, results='hide', fig.height=4, fig.width=6----
tr <- NJ(dismat)
tr <- root(tr, "Chicken_Brain", resolve.root = T)

exptable <- exptabTE(tetraexp.objects, taxa = "all",
                     subtaxa = "Brain")

f <- function(xx) {
    
     mat <- dist.pea(t(xx))
     # the distance metrics here should be the same as you specified 
     # when you created the expression distance matrix 
    
    colnames(mat) <- rownames(xx)
    rownames(mat) <- colnames(mat)
    
    root(NJ(mat), "Chicken_Brain", resolve.root = T)
    
}

bs <-  boot.phylo(tr, t(exptable), f, B = 100) 

# boot.phylo are sampling in columns of a matrix and we want to sample in rows

tr$node.label = bs
plot(tr, show.node.label = TRUE)

