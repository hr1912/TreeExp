## ---- warning=FALSE, message = FALSE-------------------------------------
library('TreeExp')

## ---- warning = FALSE, message = FALSE-----------------------------------
data('tetraexp')

## ---- warning = FALSE, message = FALSE-----------------------------------

species.group <- c("Human", "Chimpanzee", "Bonobo", "Gorilla", "Orangutan",
                   "Macaque", "Mouse", "Opossum", "Platypus")
### all mammalian species

inv.corr.mat <- corrMatInv(tetraexp.objects, taxa = species.group, subtaxa = "Brain")

## ---- warning = FALSE, message = FALSE-----------------------------------
brain.exptable <- exptabTE(tetraexp.objects, taxa = species.group, subtaxa = "Brain" )
head(brain.exptable)

## ---- warning = FALSE, message = FALSE-----------------------------------
gamma.paras <- estParaGamma(brain.exptable, inv.corr.mat)
cat(gamma.paras)

## ---- warning = FALSE, message = FALSE, fig.height=4, fig.width=6--------
   
brain.Q <- estParaQ(brain.exptable, corrmatinv = inv.corr.mat)
# with prior expression values and inversed correlation matrix
    
brain.post<- estParaWBayesian(brain.Q, gamma.paras)
brain.W <- brain.post$exp # posterior expression values
brain.CI <- brain.post$ci95 # posterior expression 95% confidence interval

names(brain.W) <- rownames(brain.exptable)

head(sort(brain.W, decreasing = T)) #check a few genes with highest seletion pressure

plot(density(brain.W))

