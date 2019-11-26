# ADnet: MCC-SP: An integration powerful method for Identification of Causal Pathways from genetic variants to complex disease
<br >


## Usage

R code avaiable from https://https://github.com/zhuyuchen95/ADnet/new/master/MCC_SP.r
```
MCC_SP(gene,cg_table1,family_xuhao,sigma, path_test=10000)
```
### Input:

  - gene: Data frame of gene expression data in Alzheimer disease pathway, each row is a sample, the row name is a uniquely identifiable id, each column is a gene, and the column name is the gene name.
  - cg_table1: Covariate data frame, which must include the following 2 columns: 1) apoe_genotype, the APOE genotype of each sample, a numerical variable, the value includes 22, 23, 24, 33, 34, 44; 2) ending: the end of the sample, It is a numerical variable of 0-1, 1 if it is AD patients, and 0 otherwise.
  - sig_gene: It is the return value sig_gene obtained by inputting the data frame gene into the linear_regression function. Can be implemented by R code linear_regression.r, sig_gene contains 3 columns, respectively gene, p, exist. p is the coefficient p value in the linear regression of the gene, exist indicates whether the p value of the gene is less than 0.05, if it is less than 0.05, it is recorded as 1, otherwise it is 0.
  - family_xuhao:It is a data frame of (m-2) * 3, where m is the number of all nodes in the network, and (m-2) represents the starting point (exposure) and end point (end) are removed, that is,the number of each gene selected from sig_gene that can represent the gene family (node). Note that these genes must all come from sig_gene. Node ordering is given for the upstream and downstream order of genes in Alzheimer disease pathway. The three columns of the data frame are: family, gene, xuhao. The first two columns are the gene family and the gene name that represents it; the last column is the corresponding node number, and the values are 2,3, ... m-1.
  - sigma: is the m * m upper triangular matrix representing the network structure, where the diagonal value is 1, and m is the number of all nodes including the start and end points. The row names are the same as the column names, both are 1: m. It is the upstream node and the column is the downstream node. If the two are connected in the network, it is copied to 5, otherwise it is 0. For example, point 1 and point 2 are directedly connected, then network [1,2] = 5, and because of the directionality of the network, network [2,1] = 0.
  - path_test: This variable is a temporary variable that determines how many paths are in the network and will be replaced by the variable 'path' later. Usually, setting 10000 in a small network (about 40 points or less) is sufficient, so the default value = 10000.


### Output

  - result_end: The return value is a list containing 5 components, which sort the paths obtained by 5 methods except Pearson-SP. Since the paths obtained by Pearson-SP name all the paths, the result of Pearson-SP is Path_1: Path_p, p is the total number of paths in the network. The name of each component is the corresponding method name, of which MCC-SP corresponds to the fifth component.
