## ---- warning=FALSE, message=FALSE---------------------------------------
library('TreeExp')

## ---- warning=FALSE, message=FALSE---------------------------------------
data('primatexp')
data('trees')

## ---- warning=FALSE, message=FALSE---------------------------------------

dismat <- expdist(primatexp.objects, taxa = "all", 
                  subtaxa = "brain", method = "sou")

primate_tree <- primatetimetree
primate_tree$tip.label <- colnames(dismat) 
# make sure their names are the same 

exp_tree <- map.ls(primate_tree, dismat) 
# map the expression distance onto the primate time tree

exp_tree <- root(exp_tree, outgroup = "Macaque_Brain", resolve.root = T)
exp_tree <- no0br(exp_tree)
# make a little tweak to the expression tree, 
# and make sure it is rooted and has no zero branch length


## ---- warning=FALSE, message=FALSE---------------------------------------
var_mat <- varMatInv(objects = primatexp.objects,phy = exp_tree,
                     taxa = "all", subtaxa = "Brain")


## ---- warning=FALSE, message=FALSE---------------------------------------
    
exp_table <- exptabTE(primatexp.objects, 
                      taxa = "all", subtaxa = "Brain")

MAG_expression <- exp_table[which(rownames(exp_table) == "ENSG00000105695"),]

## ---- warning=FALSE, message=FALSE---------------------------------------
MAG_anc <- aee(MAG_expression, exp_tree, var_mat, select = "all")


## ---- warning=FALSE, message=FALSE, fig.height=4, fig.width=6------------
primate_tree$node.label <- sprintf("%.4f",MAG_anc$est)
primate_tree$tip.label <- paste0(exp_tree$tip.label, "  ", 
                                 sprintf("%.4f", MAG_expression))

plot(primate_tree, edge.color = "grey80", edge.width = 4, 
     show.node.label = T, align.tip.label = T)

