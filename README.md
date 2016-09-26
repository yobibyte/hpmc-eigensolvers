# hmpc-eigensolvers
Final project for HMPC course (http://hpac.rwth-aachen.de/teaching/hpmc-16/)

If you just want to read the results, look at the [report](https://github.com/yobibyte/hpmc-eigensolvers/blob/master/report.md).

### How to use this repo
* install [OpenBLAS](https://github.com/xianyi/OpenBLAS) and [PAPI](http://icl.cs.utk.edu/papi/software/index.html). You will also need [numpy](http://www.numpy.org) and [jupyter notebook](http://jupyter.org) for running python code
* modify the [Makefile](https://github.com/yobibyte/hpmc-eigensolvers/blob/master/Makefile) and create 'data' and 'res' folders in this project dir
* generate the data using [this](https://github.com/yobibyte/hpmc-eigensolvers/blob/master/hmpc-data-generation.ipynb) notebook
* ./run.sh (will make the evaluate.c and do all the experiments)
* run [this](https://github.com/yobibyte/hpmc-eigensolvers/blob/master/results_analysis.ipynb) notebook to see the results of the experiments

### misc
* Some useful performance measurement code from the lecture [here](http://hpac.rwth-aachen.de/teaching/hpmc-16/sum.c).
* Timing info from lapack website [here](http://www.netlib.org/lapack/lawn41/node104.html)
* Measuring performance and accuracy of eigensolvers [paper](http://www.netlib.org/lapack/lawnspdf/lawn183.pdf)
