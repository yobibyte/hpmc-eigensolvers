evaluate: evaluate.c
	gcc -pg -o evaluate.o evaluate.c -I /opt/OpenBLAS/include -L/opt/OpenBLAS/lib -lopenblas -lgfortran
