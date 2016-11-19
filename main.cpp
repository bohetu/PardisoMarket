#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <mkl.h>

#ifndef MKL_INT
#define MKL_INT long long
#endif

#ifndef MKL_REAL
#define MKL_REAL double
#endif

#define MAX_VALUE 9999

void fillRandom(MKL_REAL *x, MKL_INT n)
{
    for (MKL_INT i = 0; i < n; i++)
        x[i] = i * 1.0 + 1.0;
    //x[i] = rand() % MAX_VALUE + 1.0;
}

MKL_INT main( void )
{
    srand(time(NULL));
    int verbose = 0;
    MKL_INT n = 5;
    MKL_INT ia[ 6] = { 1, 4, 6, 9, 12, 14 };
    MKL_INT ja[13] = { 1, 2, 4,
        1, 2,
        3, 4, 5,
        1, 3, 4,
        2, 5 };
    double a[18] = { 1.0, -1.0, -3.0,
        -2.0, 5.0,
        4.0, 6.0, 4.0,
        -4.0, 2.0, 7.0,
        8.0, -5.0 };
    
    MKL_REAL *x_rand = (MKL_REAL*) malloc (n * sizeof(MKL_REAL));
    fillRandom(x_rand, n);
    
    MKL_INT mtype = 11;
    double b[5], x[5];
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
    printf("%e \t %e \t %e\n", norm_diff, norm_xran, norm_diff / norm_xran);
    
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
