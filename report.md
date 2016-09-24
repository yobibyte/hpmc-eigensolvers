## Measuring performance of eigensolvers for small tridiagonal symmetric eigenproblems.

### Project **5** description

* You are given thousands of small tridiagonal eigenproblems (n<60)           
* Compare the accuracy and speed of three eigensolvers: BX+II (DSTEVX), QR (DSTEQR), MR3 (DSTEMR).
* Compare different matrix types:
   * Random eigenvalues, uniform distribution (0, 1)                            
   * Uniform eigenvalue distribution (page 119 of [paper](http://arxiv.org/pdf/1401.4950v1.pdf))
* Instrument the code to count flops. (use [PAPI](http://icl.cs.utk.edu/papi/))
* Study flops vs accuracy, for different accuracy levels.

### Short description of eigensolvers we test

#### QR

The general idea is to iteratively apply QR factorization to the matrix, then apply QR factorization for RxQ (R,Q are obtained from the previous iteration). In practice, an implicit approach is used as described in 7.5 paragraph of [3]. According to [2] the number of operations we need is *3bn^3+O(n^2)* where b is the average number of bulge chases and bulge chase procedure cannot use high level BLAS operations.

#### DSTEVX (BX+II)
Uses bisection method to find eigenvalues, based on [Sturm's theorem](https://en.wikipedia.org/wiki/Sturm%27s_theorem) It has *O(nk)* complexity, where k is the number of eigenvalues. If eigenvelues are clustered, it uses Gram-Schmidt orthogonalization --> it is dependent on eigenvalues distribuiton.  Again, no high-level BLAS here.[2]

#### MRRR
MRRR (MR3) algorithm is a modification of inverse iteration without Gram-Schmidt orthogonalization --> we can get *O(n^2)* Again, no high-level BLAS routines and the overall complexity depends of the eigenvalues distribution. [2]. 

### Accuracy and speed comparison

### Accuracy vs FLOPS (number of floating point operations, **not** flop per second)

### Important considerations

* All the experiments carried out for this project use synthetic data. No real problems data as [here](http://www.netlib.org/lapack/lawnspdf/lawn183.pdf) was used.
* only one architecture + tech specs

### References
* [1] http://arxiv.org/pdf/1401.4950v1.pdf
* [2] http://www.netlib.org/lapack/lawnspdf/lawn183.pdf
* [3] Golub, Van Loan, Matrix Computations, [3rd edition](http://web.mit.edu/ehliu/Public/sclark/Golub%20G.H.,%20Van%20Loan%20C.F.-%20Matrix%20Computations.pdf)
* [4] Cholesky decomposition is not needed for this report, but I wrote [this](https://github.com/yobibyte/yobiblog/issues/5) during preparation for the exam, so include it here.
