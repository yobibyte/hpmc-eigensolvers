## Measuring performance of eigensolvers for small tridiagonal symmetric eigenproblems.

### Project **5** description

* You are given thousands of small tridiagonal eigenproblems (n<60)           
* Compare the accuracy and speed of three eigensolvers: BX+II (DSTEVX), QR (DSTEQR), MR3 (DSTEMR).
* Compare different matrix types:
   * Random eigenvalues, uniform distribution (0, 1)                            
   * Uniform eigenvalue distribution (page 119 of [paper](http://arxiv.org/pdf/1401.4950v1.pdf))
* Instrument the code to count flops. (use [PAPI](http://icl.cs.utk.edu/papi/))
* Study flops vs accuracy, for different accuracy levels.

### What eigensolvers we are interested in?

#### QR
#### DSTEVX
#### DSTEMR

### Accuracy and speed comparison

### Accuracy vs FLOPS (number of floating point operations, **not** flop per second)

### Important considerations

* All the experiments carried out for this project use synthetic data. No real problems data as [here](http://www.netlib.org/lapack/lawnspdf/lawn183.pdf) was used.
* only one architecture + tech specs

### References
* http://arxiv.org/pdf/1401.4950v1.pdf
* http://www.netlib.org/lapack/lawnspdf/lawn183.pdf
