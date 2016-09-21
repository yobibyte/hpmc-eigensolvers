evaluate: evaluate.c
	gcc -p -o evaluate.o evaluate.c -I/opt/OpenBLAS/include -I/usr/local/include -L/usr/local/lib -L/opt/OpenBLAS/lib /usr/local/lib/libpapi.a -lopenblas -lgfortran
