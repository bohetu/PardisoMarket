#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <cctype>
#include <algorithm>
#include <ctime>
#include <mkl.h>
#include <string>
#include <iostream>
#include "matrixmarket.h"

#define MAX_VALUE 9999

void fillRandom(MKL_REAL *x, MKL_INT n)
{
    for (MKL_INT i = 0; i < n; i++)
        x[i] = i * 1.0 + 1.0;
    //x[i] = rand() % MAX_VALUE + 1.0;
}


bool compareMatrixEntry2(MatrixEntry const &a, MatrixEntry const &b)
{
    if (a.row == b.row)
        return a.col < b.col;
    return a.row < b.row;
}

int writeVector(const char *fname, int N, double *val)
{
    FILE *f;
    MKL_INT i;
    
    if (strcmp(fname, "stdout") == 0)
        f = stdout;
    else
        if ((f = fopen(fname, "w")) == NULL)
            return MM_COULD_NOT_WRITE_FILE;
    
    /* print banner followed by typecode */
    fprintf(f, "%%MatrixMarket matrix coordinate real unsymmetric\n");
    fprintf(f, "%d %d %d\n", N, 1, N);
    
    for (i = 0; i < N; i++)
        fprintf(f, "%lld  1  %llg\n", i + 1, val[i]);

    if (f !=stdout) fclose(f);
    
    return 0;

}

int readMatrix(const char *filename,
                MKL_INT &N_,
                MKL_INT &NZ_,
                MKL_REAL **val_,
                MKL_INT **ia_,
                MKL_INT **ja_)
{
    FILE *f;
    MM_typecode matcode;
    int M, N, nz;
    int i;
    
    if ((f = fopen(filename, "r")) == NULL)
        return -1;
    
    
    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("mm_read_unsymetric: Could not process Matrix Market banner ");
        printf(" in file [%s]\n", filename);
        return -1;
    }
    
    
    
    if ( !(mm_is_real(matcode) && mm_is_matrix(matcode) &&
           mm_is_sparse(matcode)))
    {
        fprintf(stderr, "Sorry, this application does not support ");
        fprintf(stderr, "Market Market type: [%s]\n",
                mm_typecode_to_str(matcode));
        return -1;
    }
    
    /* find out size of sparse matrix: M, N, nz .... */
    
    if (mm_read_mtx_crd_size(f, &M, &N, &nz) !=0)
    {
        fprintf(stderr, "read_unsymmetric_sparse(): could not parse matrix size.\n");
        return -1;
    }
    
    if (nz <= 0) return 0;
    
    std::vector<MatrixEntry> matrix(nz);
    
    /* reseve memory for matrices */
    
    MKL_INT *ia = (MKL_INT *) malloc((N + 1) * sizeof(MKL_INT));
    MKL_INT *ja = (MKL_INT *) malloc(nz * sizeof(MKL_INT));
    MKL_REAL *val = (MKL_REAL *) malloc(nz * sizeof(MKL_REAL));
    
    N_  = N;
    NZ_ = nz;
    *val_ = val;
    *ia_ = ia;
    *ja_ = ja;
 
    
    for (i=0; i < nz; i++)
    {
        fscanf(f, "%lld %lld %lg\n", &matrix[i].row, &matrix[i].col, &matrix[i].val);
    }
    
//    printf("------------------------------------------\n");
    
    
    std::sort(matrix.begin(), matrix.end(), compareMatrixEntry2);
    
//    for (i=0; i < nz; i++)
//    {
//        printf("%lld %lld %lg\n", matrix[i].row, matrix[i].col, matrix[i].val);
//    }
    
    ia[0] = 1;
    MKL_INT lastrow = matrix[0].row;
    MKL_INT pos = 1;
    MKL_INT inserted = 1;
    for (pos = 1; pos < nz; pos++)
    {
        while (pos < nz && matrix[pos].row == lastrow) {
            pos += 1;
        }
//        printf("pos: %lld\n", pos);
        lastrow = matrix[pos].row;
        ia[inserted] = pos + 1;
        inserted  += 1;
    }
    ia[inserted] = pos + 1;
    
//    printf("------------------------------------------\n");
    
    for (i = 0; i < nz; i++)
    {
        ja[i]  = matrix[i].col;
        val[i] = matrix[i].val;
//        printf("%lld %lld %lg\n", matrix[i].row, matrix[i].col, matrix[i].val);
    }
    
//    for (i = 0; i < N + 1; i++)
//        printf("%lld ", ia[i]);
//    printf("\n");
//    
// 
//    
//    printf("------------------------------------------\n");
//    
//
//    for (i = 0; i < nz; i++)
//    {
//        printf("%lld %lg\n", ja[i], val[i]);
//    }
//    
//    printf("------------------------------------------\n");
    fclose(f);
    
    return 0;
}

void usage(char *app)
{
    printf("usage: %s filename verbose\n", app);
}
MKL_INT main(int argc, char *argv[])
{
    if (argc != 3) {
        usage(argv[0]);
        return 1;
    }
    
    char *filename = argv[1];
    int verbose = atoi(argv[2]);
    
    srand(time(NULL));
//    MKL_INT n = 5;
//    MKL_INT ia[ 6] = { 1, 4, 6, 9, 12, 14 };
//    MKL_INT ja[13] = { 1, 2, 4,
//        1, 2,
//        3, 4, 5,
//        1, 3, 4,
//        2, 5 };
//    double a[18] = { 1.0, -1.0, -3.0,
//        -2.0, 5.0,
//        4.0, 6.0, 4.0,
//        -4.0, 2.0, 7.0,
//        8.0, -5.0 };
    
    MKL_INT n, nz;
    MKL_REAL *a;
    MKL_INT *ia, *ja;
    readMatrix(filename, n, nz, &a, &ia, &ja);
    
//    for (MKL_INT i = 0; i < nz; i++)
//        printf("%.2f ", a[i]);
//    printf("\n");
//    
//    for (MKL_INT i = 0; i < nz; i++)
//        printf("%lld ", ja[i]);
//    printf("\n");
//    
//    for (MKL_INT i = 0; i < n + 1; i++)
//        printf("%lld ", ia[i]);
//    printf("\n");
    
    MKL_REAL *b = (MKL_REAL*) malloc (n * sizeof(MKL_REAL));
    MKL_REAL *x = (MKL_REAL*) malloc (n * sizeof(MKL_REAL));
    MKL_REAL *x_rand = (MKL_REAL*) malloc (n * sizeof(MKL_REAL));
    fillRandom(x_rand, n);
    
    MKL_INT mtype = 11;
    MKL_INT nrhs = 1;
    void *pt[64];
    
    MKL_INT iparm[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    MKL_INT i;
    double ddum;
    MKL_INT idum;
    /* -------------------------------------------------------------------- */
    /* .. Setup Pardiso control parameters. */
    /* -------------------------------------------------------------------- */
    for (i = 0; i < 64; i++) {
        iparm[i] = 0;
    }
    iparm[0] = 1;
    iparm[1] = 2;
    iparm[2] = 1;
    iparm[7] = 2;
    iparm[9] = 13;
    iparm[10] = 1;
    iparm[17] = -1;
    iparm[18] = -1;
    maxfct = 1;
    mnum = 1;
    msglvl = verbose;
    error = 0;
    /* -------------------------------------------------------------------- */
    /* .. Initialize the internal solver memory pointer. This is only */
    /* necessary for the FIRST call of the PARDISO solver. */
    /* -------------------------------------------------------------------- */
    for (i = 0; i < 64; i++) {
        pt[i] = 0;
    }
    /* -------------------------------------------------------------------- */
    /* .. Reordering and Symbolic Factorization. This step also allocates */
    /* all memory that is necessary for the factorization. */
    /* -------------------------------------------------------------------- */
    phase = 11;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
            &n, a, ia, ja, &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);
    if (error != 0) {
        printf("\nERROR during symbolic factorization: %d", error);
        exit(1);
    }
    if (verbose > 0)
    {
        printf("\nReordering completed ... ");
        printf("\nNumber of nonzeros in factors = %d", iparm[17]);
        printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
    }
    /* -------------------------------------------------------------------- */
    /* .. Numerical factorization. */
    /* -------------------------------------------------------------------- */
    phase = 22;
    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
            &n, a, ia, ja, &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);
    if (error != 0) {
        printf("\nERROR during numerical factorization: %d", error);
        exit(2);
    }
    if (verbose > 0)
    {
        printf("\nFactorization completed ... ");
    }
    /* -------------------------------------------------------------------- */
    /* .. Back substitution and iterative refinement. */
    /* -------------------------------------------------------------------- */
    phase = 33;
    iparm[7] = 2; /* Max numbers of iterative refinement steps. */
    
    char trans = 'N';
    mkl_dcsrgemv(&trans, &n, a, ia, ja, x_rand, b);
    
    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
            &n, a, ia, ja, &idum, &nrhs,
            iparm, &msglvl, b, x, &error);
    if (error != 0) {
        printf("\nERROR during solution: %d", error);
        exit(3);
    }
    if (verbose > 0)
    {
        printf("\nSolve completed ... ");
        
        printf("\nThe b of the system is: ");
        for (i = 0; i < n; i++) {
            printf("\n b [%d] = % f", i, b[i] );
        }
        
        printf("\nThe solution of the system is: ");
        for (i = 0; i < n; i++) {
            printf("\n x [%d] = % f", i, x[i] );
        }
        printf ("\n");
    }
    
    MKL_REAL *diff = (MKL_REAL*) malloc (n * sizeof(MKL_REAL));
    for (MKL_INT i = 0; i < n; i++)
        diff[i] = x_rand[i] - x[i];
    
    MKL_REAL norm_diff = cblas_dnrm2(n, diff, 1.0);
    MKL_REAL norm_xran = cblas_dnrm2(n, x_rand, 1.0);
    printf("Filename: %s\n", filename);
    printf("Error: %e\n",norm_diff / norm_xran);
    
    
    std::string stdFileName = filename;
    std::string prefix = stdFileName.substr(0, stdFileName.find("A.mtx"));
    std::string xgenName = prefix + "random-x.mtx";
    std::string xsolName = prefix + "solution-x.mtx";
    std::string bVecName = prefix + "rhs-b.mtx";
    
    writeVector(xgenName.c_str(), n, x_rand);
    writeVector(xsolName.c_str(), n, x);
    writeVector(bVecName.c_str(), n, b);
    
    /* -------------------------------------------------------------------- */
    /* .. Termination and release of memory. */
    /* -------------------------------------------------------------------- */
    phase = -1; /* Release internal memory. */
    pardiso(pt, &maxfct, &mnum, &mtype, &phase,
            &n, &ddum, ia, ja, &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);
    
    free(x_rand);
    free(diff);
    return 0;
}
