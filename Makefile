evaluate: evaluate.c
	gcc -pg -o evaluate evaluate.c -I /opt/OpenBLAS/include -L/opt/OpenBLAS/lib /usr/lib/ -lopenblas -llapacke -llapack -lgfortran
