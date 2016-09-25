# Measuring performance of eigensolvers for small tridiagonal symmetric eigenproblems.

## Project description

* You are given thousands of small tridiagonal eigenproblems (n<60)           
* Compare the accuracy and speed of three eigensolvers: BX+II (DSTEVX), QR (DSTEQR), MR3 (DSTEMR).
* Compare different matrix types:
   * Random eigenvalues, uniform distribution (0, 1)                            
   * Uniform eigenvalue distribution (page 119 of [paper](http://arxiv.org/pdf/1401.4950v1.pdf))
* Instrument the code to count flops. (use [PAPI](http://icl.cs.utk.edu/papi/))
* Study flops vs accuracy, for different accuracy levels.

## Short description of eigensolvers we test

### QR

The general idea is to iteratively apply QR factorization to the matrix, then apply QR factorization for RxQ (R,Q are obtained from the previous iteration). In practice, an implicit approach is used as described in 7.5 paragraph of [3]. According to [2] the number of operations we need is *3bn^3+O(n^2)* where b is the average number of bulge chases and bulge chase procedure cannot use high level BLAS operations.

### DSTEVX (BX+II)
Uses bisection method to find eigenvalues, based on [Sturm's theorem](https://en.wikipedia.org/wiki/Sturm%27s_theorem) It has *O(nk)* complexity, where k is the number of eigenvalues. If eigenvelues are clustered, it uses Gram-Schmidt orthogonalization --> it is dependent on eigenvalues distribuiton.  Again, no high-level BLAS here.[2]

### MRRR
MRRR (MR3) algorithm is a modification of inverse iteration without Gram-Schmidt orthogonalization --> we can get *O(n^2)* Again, no high-level BLAS routines and the overall complexity depends of the eigenvalues distribution. [2]. 

## Data generation.

The data was generated in [Jupyter notebook](https://github.com/yobibyte/hpmc-eigensolvers/blob/master/hmpc-data-generation.ipynb). Initially the eigenvalues of needed dimension were generated (either sampled from standard random uniform distribution, or taken uniformly from [0,1]. Different distributions were taken, I suppose, to check the dependency of BX+II and MRRR on the eigenvalues clustering. Then, I got random orthogonal matrices Q of the same dimensions and got the problems by Q*M*Q' (for orthogonal matrices the inverse is equivalent to transposed). The equality of initial eigenvalues and solutions to given problems was checked using numpy.testing library.

The generated data was then written to .csv files and read inside C program.

## Code

Each experiment represent a separate run of a program that tests particular data file (10000 problems per file) with different problem dimension (10, 20, 30, 40 and 50), eigenvalue distribution (random uniform and uniform eigenvalue distributions) and algorithm (BX+II, MRRR, QR) to test. Inside each experiment there is a for loop with one iteration per problem. Before start of each iteration the garbage of L3 cache size was created to eliminate cache out of the experiments. PAPI library was used to check the real time of algorithm execution and the number of floating point operation per algorithm call. LAPACKE interface inside OpenBLAS was used, code was compiled with gcc.

## Results

### Accuracy and speed comparison

### Accuracy vs FLOPS (number of floating point operations, **not** flop per second)

* Due to specific of DSTEQR, there is no relative tolerance parameter which defines when the problem is considered as solved. 
* As for DSTEMR, there is TRYRAC input parameter, that (if it is true) will make the procedure to check if the tridiagonal matrix defines the eigenvalues to high relative accuracy ([documentation](http://www.netlib.org/lapack/explore-html/d9/d1e/dstemr_8f_a613f73c16db5b9b111d56fb3e3feff0d.html#a613f73c16db5b9b111d56fb3e3feff0d)]. But using it did not show any difference neither in execution time nor in accuracy obtained.

There also was a possibility to check DSTARRV FORTRAN procedure, but I did not get exactly how should I do it (OpenBLAS has no access to it) and I did not modify FORTRAN code.

So, only DSTEVX that has relative tolerance parameter has adequate results for this experiment.


## Important considerations

* All the experiments carried out for this project use synthetic data. No real problems data as [here](http://www.netlib.org/lapack/lawnspdf/lawn183.pdf) was used.
* only one architecture + tech specs

## References
* [1] http://arxiv.org/pdf/1401.4950v1.pdf
* [2] http://www.netlib.org/lapack/lawnspdf/lawn183.pdf
* [3] Golub, Van Loan, Matrix Computations, [3rd edition](http://web.mit.edu/ehliu/Public/sclark/Golub%20G.H.,%20Van%20Loan%20C.F.-%20Matrix%20Computations.pdf)
* [4] Cholesky decomposition is not needed for this report, but I wrote [this](https://github.com/yobibyte/yobiblog/issues/5) during preparation for the exam, so include it here.
