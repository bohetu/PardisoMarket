all:
	icpc main.cpp matrixmarket.cpp -o pardiso_solver -mkl -O3
