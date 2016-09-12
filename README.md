TreeExp
=======

*TreeExp* is an *R* package that performs analyses of expression
evolution from *RNA-seq* data, including optimized input formatting,
normalization and pair-wise distance evaluation, expression character
tree inference and preliminary phylogenetic network analysis.

*TreeExp* package is under active developing, current stable version 1.0
is available at <https://github.com/hr1912/TreeExp>.

A convenient way to install package from github is through *devtools*
package:

    install.packages('devtools')
    devtools::install_github("hr1912/TreeExp")

Users can also download *TreeExp* package and install locally through:

    install.packages("filePath/TreeExp.1.0.tar.gz", repos = NUll, type = "source")

Load the package in the usual way:

    library('TreeExp')

    ## Loading required package: ape

Expression divergence of two genomes can be considered as concerted
evolution of transcriptome from their common ancestor, which can be
measured by among orthologous genes' variance components. And more
specifically, the expression levels of same tissue across different
species can be treated as taxonomic units, which can be placed at the
tips of character tree that represents expression evolution of these
species.

In here, we will give an example to build a character tree from
expression data (expression phylogeny).

### Input Format:

*TreeExp* package takes in reads count data and gene information file in
certain format:

1.  Gene information file should be a text file in the shape of a
    matrix, in which values are separated by tabs. `Rows` correspond to
    orthologous genes and `columns` correspond to species names. And the
    values in the matrix are in the format of "`GeneId:GeneLength`".

2.  Reads count file should also be a text file in the matrix shape,
    `Rows` correspond to orthologous genes which should be in one-to-one
    correspondence with rows in Gene information file, though gene ids
    are displayed in reads count file. `Columns` correspond to
    sample names. Sample names are in format of
    "`TaxaName_SubtaxaName_ReplicatesName`".

The example files are included in the TreeExp package, which can be
found in `extdata` folder in the package. One can load them in to take a
look:

    readsCount.table = read.table(system.file('extdata/tetraexp.reads.count.raw.txt', 
                                                package='TreeExp'), header = T)
    head(readsCount.table[,1:10])

    ##   homoSapienGeneId Human_Brain_Female Human_Brain_Male1 Human_Brain_Male2
    ## 1  ENSG00000198824              16323             11147             19507
    ## 2  ENSG00000118402              31883             19242             32321
    ## 3  ENSG00000166167             102711             80104            141338
    ## 4  ENSG00000144724              51020             37861             43906
    ## 5  ENSG00000183508                988              1443               760
    ## 6  ENSG00000008086              27839             10208             20666
    ##   Human_Brain_Male3 Human_Brain_Male4 Human_Brain_Male5
    ## 1              4402              6521               683
    ## 2              9841              1164              1362
    ## 3             59947             23104             14613
    ## 4             24645              3968              5653
    ## 5               456               227               456
    ## 6             12970              2569              1940
    ##   Human_Cerebellum_Female Human_Cerebellum_Male Human_Heart_Female
    ## 1                   20032                 43008               3381
    ## 2                   26732                 57081                228
    ## 3                   56103                105866              17022
    ## 4                   26365                 66177              41876
    ## 5                     988                  1215              16913
    ## 6                    4987                 11862                878

    geneInfo.table = read.table(system.file('extdata/tetraexp.gene.length.ortholog.txt',
                                            package='TreeExp'), header = T)
    head(geneInfo.table)

    ##                  Human              Chimpanzee                 Gorilla
    ## 1 ENSG00000198824:3788 ENSPTRG00000023033:3786 ENSGGOG00000009918:3808
    ## 2 ENSG00000118402:3042 ENSPTRG00000018370:2962 ENSGGOG00000009359:2890
    ## 3 ENSG00000166167:6255 ENSPTRG00000002869:6153 ENSGGOG00000005023:5882
    ## 4 ENSG00000144724:9524 ENSPTRG00000015067:6726 ENSGGOG00000006858:6049
    ## 5 ENSG00000183508:5751 ENSPTRG00000001166:5717 ENSGGOG00000004149:5727
    ## 6 ENSG00000008086:3576 ENSPTRG00000021711:2874 ENSGGOG00000010796:3354
    ##                 Orangutan                 Macaque                   Mouse
    ## 1 ENSPPYG00000005537:2442 ENSMMUG00000023599:3781 ENSMUSG00000047710:4040
    ## 2  ENSPPYG00000016794:945 ENSMMUG00000020208:2960 ENSMUSG00000032262:2146
    ## 3 ENSPPYG00000002583:6130 ENSMMUG00000006741:6001 ENSMUSG00000025217:2984
    ## 4 ENSPPYG00000013752:4933 ENSMMUG00000012487:5949 ENSMUSG00000021745:9375
    ## 5 ENSPPYG00000000972:1987 ENSMMUG00000004907:3522 ENSMUSG00000044468:5640
    ## 6 ENSPPYG00000020166:3024 ENSMMUG00000005063:3077 ENSMUSG00000031292:3484
    ##                   Opossum                Platypus                 Chicken
    ## 1 ENSMODG00000003128:2568 ENSOANG00000000610:2487 ENSGALG00000016813:2412
    ## 2  ENSMODG00000018420:945  ENSOANG00000001290:909 ENSGALG00000015876:1322
    ## 3 ENSMODG00000011788:1888 ENSOANG00000007316:1936 ENSGALG00000007820:1962
    ## 4 ENSMODG00000002722:4448 ENSOANG00000006011:3968 ENSGALG00000007177:5396
    ## 5 ENSMODG00000023227:1194 ENSOANG00000003065:1335 ENSGALG00000014453:1179
    ## 6 ENSMODG00000017140:2589 ENSOANG00000004037:2607 ENSGALG00000016529:2712

### Construction:

The construction function `TEconstruct` loads in the reads count data
file as well as a gene information file, and wraps them in a list of
*taxonExp* objects (one *taxaExp* object).

In the package, we include files transformed from six tissues'
expression reads count data of nine tetrapod species. If you want to
transform your own data, a transformation Perl script
`format2treeexp.pl` to format raw outputs of *TopHat2* to "*TreeExp*
compatible" is available at `tools` folder in the package. Or you can
access the script at
<https://github.com/hr1912/TreeExp/blob/master/tools/format2treeexp.pl>

    taxa.objects = TEconstruct(readsCountFP = system.file('extdata/tetraexp.reads.count.raw.txt', package='TreeExp'),
      geneInfoFP = system.file('extdata/tetraexp.gene.length.ortholog.txt', package='TreeExp'), 
      taxa = "all", subtaxa = c("Brain", "Cerebellum"), calRPKM=TRUE, rmOut=TRUE)

The construction process takes **several minutes** on a desktop computer
depending on data size and hardware performance. Specify **"taxa"** and
**"subtaxa"** options in the function when using partial of your data.
The construction process will be faster. If you are hesitated to test
the *TreeExp*, the package has already bundled a constructed object and
you can load the object through:

    data(tetraexp)

You can take a look at what the loaded objects:

    print(tetraexp.objects, details = TRUE)

    ## 
    ##  53 taxonExp objects 
    ## 
    ## object 1 : Human      Brain 
    ## object 2 : Human      Cerebellum 
    ## object 3 : Human      Heart 
    ## object 4 : Human      Kidney 
    ## object 5 : Human      Liver 
    ## object 6 : Human      Testis 
    ## object 7 : Chimpanzee     Brain 
    ## object 8 : Chimpanzee     Cerebellum 
    ## object 9 : Chimpanzee     Heart 
    ## object 10 : Chimpanzee    Kidney 
    ## object 11 : Chimpanzee    Liver 
    ## object 12 : Chimpanzee    Testis 
    ## object 13 : Gorilla   Brain 
    ## object 14 : Gorilla   Cerebellum 
    ## object 15 : Gorilla   Heart 
    ## object 16 : Gorilla   Kidney 
    ## object 17 : Gorilla   Liver 
    ## object 18 : Gorilla   Testis 
    ## object 19 : Orangutan     Brain 
    ## object 20 : Orangutan     Cerebellum 
    ## object 21 : Orangutan     Heart 
    ## object 22 : Orangutan     Kidney 
    ## object 23 : Orangutan     Liver 
    ## object 24 : Macaque   Brain 
    ## object 25 : Macaque   Cerebellum 
    ## object 26 : Macaque   Heart 
    ## object 27 : Macaque   Kidney 
    ## object 28 : Macaque   Liver 
    ## object 29 : Macaque   Testis 
    ## object 30 : Mouse     Brain 
    ## object 31 : Mouse     Cerebellum 
    ## object 32 : Mouse     Heart 
    ## object 33 : Mouse     Kidney 
    ## object 34 : Mouse     Liver 
    ## object 35 : Mouse     Testis 
    ## object 36 : Opossum   Brain 
    ## object 37 : Opossum   Cerebellum 
    ## object 38 : Opossum   Heart 
    ## object 39 : Opossum   Kidney 
    ## object 40 : Opossum   Liver 
    ## object 41 : Opossum   Testis 
    ## object 42 : Platypus      Brain 
    ## object 43 : Platypus      Cerebellum 
    ## object 44 : Platypus      Heart 
    ## object 45 : Platypus      Kidney 
    ## object 46 : Platypus      Liver 
    ## object 47 : Platypus      Testis 
    ## object 48 : Chicken   Brain 
    ## object 49 : Chicken   Cerebellum 
    ## object 50 : Chicken   Heart 
    ## object 51 : Chicken   Kidney 
    ## object 52 : Chicken   Liver 
    ## object 53 : Chicken   Testis

    print(tetraexp.objects[[1]], printlen = 6)

    ## 
    ## One taxonExp object
    ## Taxon name:  Human 
    ## Subtaxon name:  Brain 
    ## Total gene number:  5636 
    ## Total bio replicates number:  6 
    ## Bio replicates titles:
    ## [1] "Female" "Male1"  "Male2"  "Male3"  "Male4"  "Male5" 
    ## Outliers removed
    ## RPKM calculated
    ## Over-dispersion parameter omega NOT calculated

    head(tetraexp.objects[[1]]$rpkm.rmOut)

    ##                 Human_Brain_Female Human_Brain_Male1 Human_Brain_Male2
    ## ENSG00000198824          3.9537403          2.700015         4.7249655
    ## ENSG00000118402          9.6165234          5.803756         9.7486326
    ## ENSG00000166167         15.0663402         11.750193        20.7324083
    ## ENSG00000144724          4.9151773          3.647462         4.2298271
    ## ENSG00000183508          0.1576274          0.230219         0.1212518
    ## ENSG00000008086          7.1428947          2.619155         5.3024556
    ##                 Human_Brain_Male3 Human_Brain_Male4 Human_Brain_Male5
    ## ENSG00000198824        1.06624791        1.57950991        0.16543556
    ## ENSG00000118402        2.96823408        0.35108469        0.41080529
    ## ENSG00000166167        8.79342907        3.38905008        2.14353310
    ## ENSG00000144724        2.37425609        0.38227016        0.54460011
    ## ENSG00000183508        0.07275111        0.03621601        0.07275111
    ## ENSG00000008086        3.32782585        0.65915070        0.49776269

### Distance matrix:

Let us quickly jump to the issue of creating phylogeny from *taxaExp*
object. First, we generate a distance matrix:

    dismat <- expdist(tetraexp.objects, taxa = "all",
                     subtaxa = "Brain",
                     method = "pea")
    dismat

    ##                  Human_Brain Chimpanzee_Brain Gorilla_Brain
    ## Human_Brain       0.00000000       0.00000000    0.00000000
    ## Chimpanzee_Brain  0.04272425       0.00000000    0.00000000
    ## Gorilla_Brain     0.06302204       0.05099573    0.00000000
    ## Orangutan_Brain   0.09460563       0.07342807    0.06543296
    ## Macaque_Brain     0.09178295       0.08056567    0.06950483
    ## Mouse_Brain       0.16133517       0.15274712    0.15187221
    ## Opossum_Brain     0.23696315       0.22046195    0.20750404
    ## Platypus_Brain    0.27608779       0.26075046    0.25261355
    ## Chicken_Brain     0.29903270       0.29092886    0.26943178
    ##                  Orangutan_Brain Macaque_Brain Mouse_Brain Opossum_Brain
    ## Human_Brain           0.00000000     0.0000000   0.0000000     0.0000000
    ## Chimpanzee_Brain      0.00000000     0.0000000   0.0000000     0.0000000
    ## Gorilla_Brain         0.00000000     0.0000000   0.0000000     0.0000000
    ## Orangutan_Brain       0.00000000     0.0000000   0.0000000     0.0000000
    ## Macaque_Brain         0.07304648     0.0000000   0.0000000     0.0000000
    ## Mouse_Brain           0.14767685     0.1219459   0.0000000     0.0000000
    ## Opossum_Brain         0.20089057     0.1860964   0.1752627     0.0000000
    ## Platypus_Brain        0.24856272     0.2385429   0.2341083     0.2024356
    ## Chicken_Brain         0.26518981     0.2420884   0.2394585     0.2593075
    ##                  Platypus_Brain Chicken_Brain
    ## Human_Brain           0.0000000             0
    ## Chimpanzee_Brain      0.0000000             0
    ## Gorilla_Brain         0.0000000             0
    ## Orangutan_Brain       0.0000000             0
    ## Macaque_Brain         0.0000000             0
    ## Mouse_Brain           0.0000000             0
    ## Opossum_Brain         0.0000000             0
    ## Platypus_Brain        0.0000000             0
    ## Chicken_Brain         0.2757437             0

You can specify **"taxa"** and **"subtaxa"** options in the `expdist`
function as well. The default model **"pea"** is to calculate pair-wise
distances by Pearson distance, which equals 1-Pearson’s coefficient of
expression level.

### Expression character tree:

After distance matrix is created, you can construct character tree by
Neighbor-Joining, and bootstrap values can also be generated by
`boot.exphy` function:

    tr <- NJ(dismat)
    tr <- root(tr, "Chicken_Brain")
    bs <- boot.exphy(tr, tetraexp.objects, method = "pea",
                     B = 100, rooted = "Chicken_Brain")
    tr$node.label = bs
    plot(tr, show.node.label = TRUE)

![](README_files/figure-markdown_strict/unnamed-chunk-10-1.png)

By now, an expression tree is successfully constructed. The tree shows
expression patterns' similarities in selected genes of designated
species. The expression tree is largely in accordance with species tree
with minor discrepancy.

Phenomenon of evolutionary history dominates the evolutionary expression
pattern can be described as phylogenetic signals. One way to interpret
highly consistent expression character tree is that expression levels of
transcriptome, representing the regulatory changes, accumulated over
time. Though not as concrete as sequence data, expression levels
generated from transcriptome data across species show strong
phylogenetic signals.
