# MarkovBases25years

Supplementary material for examples in the paper *Markov bases: a 25 year update*. 

-----


This page contains the code necessary to reproduce the few examples in the paper. For sampling fibers of log-linear models for testing goodness of fit, we rely on the   `algstat` package can be  obtained from [github.com/dkahle/algstat](https://github.com/dkahle/algstat). We have installed it using the following commands in `R`: 
```
if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("dkahle/mpoly")
devtools::install_github("coneill-math/m2r")
devtools::install_github("dkahle/latte")
devtools::install_github("dkahle/bertini")
devtools::install_github("dkahle/algstat")
``` 


### Example 2.1: job satisfaction 

The following $4\times4$ contingency table is from David Kahle's algstat github page, Agresti (2002, p.57).
```
Job <- matrix(
  c(1,2,1,0, 3,3,6,1, 10,10,14,9, 6,7,12,11), nrow = 4, ncol = 4,
  dimnames = list(
    "income" = c("< 15k", "15-25k", "25-40k", "> 40k"),
    "satisfaction" = c("VeryD", "LittleD", "ModerateS", "VeryS")
  )
)
```

The design matrix $A$ for the independence model on 4x4 table can be computed as follows:
```
library("algstat")
A <- hmat(c(4,4), 1:2)
```
The vector of sufficient statistics is obtained by computing $Au$ where $u$ is the vectorized version of the table `Job`.
```
b <- t(A %*% c(Job)) 
# Let us find a random move on the fiber: 
m<-rmove(1, A = A, b = b)
# and check that its margins are zero: 
t(A%*%m)
``` 

### Example 3.13: exploring fibers for testing goodness of fit 

If the `algstat` package was installed properly along with all of the algebraic packages that are required, the following command can be run directly from R:
```
loglinear(~ income + satisfaction, data = Job)
```
By default, `loglinear` computes a Markov basis.  Readers who do not have the required software packages installed should follow the next section. 

> Our readers should note that the `algstat` package interfaces with algebraic and combinatorial software packages (`Macaulay2`, `4ti2`, `LattE`) to compute Markov bases and polytopes. However, under standard installations, it is entirely possible that those packages are not avaiable. Therefore, we provide code to compute the necessary bases using an online interface, and then manually import them into the `algstat` package in `R`. 

#### Computing Markov bases using algebraic software directly 

The following code in `Macaulay2` can be run through  [this web interface](https://www.unimelb-macaulay2.cloud.edu.au/) without installing the software:

* Press the orange 'start' button on the website
* once the input line `i1:` appears, copy the following code: 

```
-- define the design matrix A: 
A=matrix{
{1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1},
{1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0},
{0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0},
{0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0},
{0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1}}

-- load the required package for interfacing with 4ti2 for Markov bases computations: 
loadPackage "FourTiTwo"

-- compute the Markov basis; it will be returned as rows of the following matrix: 
toricMarkov A 
```

We have saved the output of the line above to a file and imported it into `R` manually below. 
```
# matrix whose columns are Markov basis elements: 
markovJob <- t(matrix (c( 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, -1, 1,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, -1, 1, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, -1, 0, 1,
  0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, -1, 1, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -1, 0, -1, 0, 1, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, -1, 0, 0, 1,
  0, 0, 0, 0, 0, 0, 1, -1, 0, 0, -1, 1, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, -1, 1,
  0, 0, 0, 0, 0, 1, -1, 0, 0, -1, 1, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, -1, 1, 0,
  0, 0, 0, 0, 0, 1, 0, -1, 0, -1, 0, 1, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, -1, 0, 1,
  0, 0, 0, 0, 1, -1, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0,
  0, 0, 0, 0, 1, 0, -1, 0, -1, 0, 1, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0, -1, 0, 1, 0,
  0, 0, 0, 0, 1, 0, 0, -1, -1, 0, 0, 1, 0, 0, 0, 0,
  0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 1,
  0, 0, 1, -1, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 1, -1, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0,
  0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1,
  0, 1, -1, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 1, -1, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0,
  0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0,
  0, 1, 0, -1, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 1, 0, -1, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0,
  0, 1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1,
  1, -1, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, -1, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0,
  1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0,
  1, 0, -1, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 0, -1, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0,
  1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0,
  1, 0, 0, -1, -1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 1, 0, 0, 0, 0,
  1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 1) ,ncol=16,byrow=TRUE))
# check that the rank of the Markov basis matrix is 9: 
library("Matrix"); 
rankMatrix(markovJob) 
```

Let us verify that each column of `markovJob` is a move for the 4x4 model of independence: 
```
A%*%markovJob # the margins are all zero. 
```

#### Sparse table with low cell counts: histograms of acceptance probabilites 

```
DiaconisSturmfels <- matrix(
  c(1,1,0,0, 0,1,1,0, 0,0,1,1, 1,0,0,1), nrow = 4, ncol = 4,
  dimnames = list(
    "income" = c("< 15k", "15-25k", "25-40k", "> 40k"),
    "satisfaction" = c("VeryD", "LittleD", "ModerateS", "VeryS")
  )
)
latticeBasisDS <- t(matrix( 
  c(1,-1,0,0, -1,1,0,0, 0,0,0,0, 0,0,0,0,
    0,1,-1,0, 0,-1,1,0, 0,0,0,0, 0,0,0,0,
    0,0,1,-1, 0,0,-1,1,  0,0,0,0, 0,0,0,0,  
    0,0,0,0, 1,-1,0,0, -1,1,0,0, 0,0,0,0,
    0,0,0,0, 0,1,-1,0, 0,-1,1,0, 0,0,0,0, 
    0,0,0,0, 0,0,1,-1, 0,0,-1,1,  0,0,0,0, 
    0,0,0,0, 0,0,0,0, 1,-1,0,0, -1,1,0,0,
    0,0,0,0, 0,0,0,0,0,1,-1,0, 0,-1,1,0,
    0,0,0,0, 0,0,0,0,0,0,1,-1, 0,0,-1,1
    ),
  nrow = 16, ncol = 9))
rankMatrix(latticeBasisDS)

# here is one run of the Markov chain fiber sample for each of the bases: 
# Markov basis: 
outMB<- loglinear(~ income + satisfaction, data = DiaconisSturmfels, moves=markovJob)
outMB
# accept_prob: the empirical transition probability of the moves, including the thinned moves. 
outMB$accept_prob
# lattice basis: 
outLB<-loglinear(~ income + satisfaction, data = DiaconisSturmfels, moves=latticeBasisDS)
outLB
outLB$accept_prob
```

The following code was used to produce the histograms in Figure 1: 
```
outLoopLattice<-matrix(c(0),nrow=1,ncol=10)
for(i in c(1:10)) (
  outLoopLattice[,i] <- loglinear(~ income + satisfaction, 
  data = DiaconisSturmfels, moves=latticeBasisDS)$accept_prob 
) 
hist(outLoopLattice[1,],xlim=c(0,1),ylim=c(0,5),main="Lattice basis",xlab="Acceptance probability")
outLoopMarkov<-matrix(c(0),nrow=1,ncol=10)
for(i in c(1:10)) (
  outLoopMarkov[,i] <-loglinear(~ income + satisfaction, 
  data = DiaconisSturmfels, moves=markovJob)$accept_prob 
) 
hist(outLoopMarkov[1,],xlim=c(0,1),ylim=c(0,5), main="Markov basis",xlab="Acceptance probability")
```

For completeness, the reader may wish to check this data set using standard (not algebraic statistics) techniques:
```
fisher.test(DiaconisSturmfels)
MASS::loglm(~ income + satisfaction, data = DiaconisSturmfels)
```


#### Sparse table with larger cell counts: histograms of p-values  

The following is the second data table in the example.  
```
sparseTable <- matrix(
  c(10,0,0,2, 0,3,0,40, 10,0,2,0, 0,3,40,0), nrow = 4, ncol = 4,
  dimnames = list(
    "income" = c("< 15k", "15-25k", "25-40k", "> 40k"),
    "satisfaction" = c("VeryD", "LittleD", "ModerateS", "VeryS")
  )
)
``` 

Unlike the previous table, this one does not fit the model of independence. Again, for completeness, we run the following. 
```
fisher.test(sparseTable)
MASS::loglm(~ income + satisfaction, data = sparseTable)
```
Let us perform one test run using Markov and lattice bases:
```
# Markov basis: 
testrunMB<-loglinear(~ income + satisfaction, data = sparseTable, moves=markovJob)
testrunMB$p.value # the second p-value refers to the chi-square statistic 
testrunMB$accept_prob

# lattice basis: 
testrun<-loglinear(~ income + satisfaction, data = sparseTable, moves=latticeBasisDS)
testrun$p.value # the second p-value refers to the chi-square statistic 
testrun$accept_prob
```

The following code was used to produce the histograms in Figure 2: 
```
outLoopLattice<-matrix(c(0),nrow=6,ncol=100)
for(i in c(1:100)) (
  outLoopLattice[,i] <-loglinear(~ income + satisfaction, 
  data = sparseTable, moves=latticeBasisDS)$p.value
) 
hist(outLoopLattice[2,],xlim=c(0,0.35),main="Lattice bases",xlab="p values for the chi-square statistic")

 
outLoopMarkov<-matrix(c(0),nrow=6,ncol=100)
for(i in c(1:100)) (
  outLoopMarkov[,i] <-loglinear(~ income + satisfaction, 
  data = sparseTable, moves=markovJob)$p.value
) 
hist(outLoopMarkov[2,],xlim=c(0,0.1),main="Markov bases",xlab="p values for the chi-square statistic")
```


### Computing fiber sizes 

The following `python` code was used to compute fiber sizes for the three $4\times 4$ tables in the examples above.
```
import networkx as nx
import numpy as np
from PyNormaliz import *
import time

def get_equalities(A, b):
    """"
    Input
    - A: n x m integer matrix 
    - b: m-dim integer vector
    Output
    - List of equalities Ax=b in Normaliz format
    """
    b = -np.array(b, dtype=int)
    A = np.array(A, dtype=int)
    eqs = np.column_stack((A, b))
    return [list(x) for x in eqs]

if __name__ == "__main__":

    # (Example 2.1) Job satisfaction table
    table_1 = np.array([[1, 3, 10, 6],
                        [2, 3, 10, 7],
                        [1, 6, 14, 12],
                        [0, 1, 9, 11]])

    # (Example 3.13) 
    table_2 = np.array([[1, 0, 0, 1],
                        [1, 1, 0 ,0],
                        [0, 1, 1, 0],
                        [0, 0, 1, 1]])

    # (Example 3.13)
    table_3 = np.array([[10, 0, 10, 0],
                        [0, 3, 0, 3],
                        [0, 0, 2, 40],
                        [2, 40, 0, 0]])

    tables = [table_1, table_2, table_3]

    file_name = "computed_data.txt"  # File to store computed data
    with open(file_name, "w") as file:
        for i in range(len(tables)):
            table = tables[i]
            nrow, ncol = table.shape

            # We use the library NetworkX library in python to recover the design matrix for the
            # independence model with 2 variables and levels (nrow, ncol)
            A = nx.incidence_matrix(nx.complete_bipartite_graph(nrow, ncol)) 
            b = A @ table.flatten()
            eqs = get_equalities(A.toarray(), b)

            start_time = time.time()
            C = Cone(inhom_equations=eqs)
            fiber_size = C.NumberLatticePoints()
            end_time = time.time()
            computation_time = end_time - start_time

            # Prepare the data to store in the file
            data_to_store = f"Table {i+1}:\n {table}\nFiber size: {fiber_size}\nComputation Time: {computation_time}\n\n"
            print(data_to_store)

            file.write(data_to_store)
```

### Example 3.10: example of a non-decomposable model 

If the software package `4ti2` is installed, then the following code illustrates example 3.10 in the paper.

```
library("latte")
library("algstat")
# define the matrices: 
A <- rbind(c(1,1,1,0,0,0),
           c(0,0,0,1,1,1), 
           c(1,0,0,1,0,0), 
           c(0,1,0,0,1,0),
           c(0,0,1,0,0,1))
B <- diag(6)

# To compute g(A,B) we get the Graver basis of B*Gr(A) and retrieve the max 1-norm of its elements 
GrA = graver(A)
BGrA = graver(B %*% GrA)

# Computing g(A, B)
grav_complexity <- max(apply(abs(BGrA), 2, sum))
paste("Graver complexity g(A,B): ", grav_complexity)

# Retrieving the nfold matrix [A,B]^(g) where g is the Graver complexity computed above 
nfold_AB <- rbind(
  cbind(diag(grav_complexity) %x% A),   # Creating a block diagonal matrix with 'grav_complexity' copies of A
  cbind(do.call(cbind, replicate(grav_complexity, B, simplify = FALSE)))  # Concatenating 'grav_complexity' copies of B in the last row block
)

Gr_nfold_AB <- graver(nfold_AB)
paste("|Gr(A,B)^(g)|: ", ncol(Gr_nfold_AB))
```
