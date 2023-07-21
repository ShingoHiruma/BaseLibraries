

#include "MatSolvers.hpp"

#ifdef INTEL_MKL_SOLVER_USING

#include "SparseMatOperators.hpp"
#include "mkl.h"
#include <omp.h>

/* 専用名前空間 */
namespace SRLfem{


/*

//=======================================================
//=======================================================
//=======================================================
Eigenソルバ
//=======================================================
//=======================================================

*/


/*//=======================================================
// ● MKL Pardisoで解く(非対称)
//=======================================================*/
bool MatSolvers::solveMLKpardiso(const slv_int size0, const SparseMat& matA, double* vecB, double *results, int num_para){
    const int n = size0;
    /* スパース形式コピー */
    int *ia = new int[n+1];
    auto row_ptr = matA.matrix->matrix.outerIndexPtr();
    const int total = matA.matrix->matrix.nonZeros();
    int* ja = new int[total];
    auto col_ptr = matA.matrix->getColPtr();
    double* a = new double[total];
    auto val_ptr = matA.matrix->getValuePtr();

    for(int i = 0; i < n; ++i){
        ia[i] = row_ptr[i];
        ia[i]++;
    }
    ia[size0] = total+1;
    for (int i = 0; i < total; ++i){
        ja[i] = col_ptr[i];
        ja[i]++;
        a[i] = val_ptr[i];
    }

    /* 実非対称設定 */
    int mtype = 11;       
    /* 右辺コピー */
    double* b = new double[n];
    for(int i = 0; i < n; ++i){
        b[i] = vecB[i];
    }

    int nrhs = 1; 
    void *pt[64];    
    int iparm[64];
    int maxfct, mnum, phase, error, msglvl;
    double ddum;    /* Double dummy */
    int idum;       /* Integer dummy. */

    omp_set_num_threads(num_para);
    for (int i = 0; i < 64; i++ ){
        iparm[i] = 0;
    }
    iparm[0] = 1;         /* No solver default */
    iparm[1] = 2;         /* Fill-in reordering from METIS */
    iparm[2] = num_para;
    iparm[3] = 0;         /* No iterative-direct algorithm */
    iparm[4] = 0;         /* No user fill-in reducing permutation */
    iparm[5] = 0;         /* Write solution into x */
    iparm[6] = 0;         /* Not in use */
    iparm[7] = 0;         /* Max numbers of iterative refinement steps */
    iparm[8] = 0;         /* Not in use */
    iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;        /* Not in use */
    iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
    iparm[13] = 1;        /* Output: Number of perturbed pivots */
    iparm[14] = 0;        /* Not in use */
    iparm[15] = 0;        /* Not in use */
    iparm[16] = 0;        /* Not in use */
    iparm[17] = 1;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = 1;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */
    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;             /* Which factorization to use. */
    msglvl = 1;           /* Print statistical information in file */
    error = 0;            /* Initialize error flag */    
    for (int i = 0; i < 64; i++ ){
        pt[i] = 0;
    }
    /* 処理１ */
    phase = 11;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 ){
        delete[] ia;
        delete[] ja;
        delete[] a;
        delete[] b;
        return false;
    }
    /* 処理２ */
    phase = 22;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 ){
        delete[] ia;
        delete[] ja;
        delete[] a;
        delete[] b;
        return false;
    }    
    /* 処理３ */
    phase = 33;
    iparm[7] = 2;   /* Max numbers of iterative refinement steps. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, results, &error);
    if ( error != 0 ){
        delete[] ia;
        delete[] ja;
        delete[] a;
        delete[] b;
        return false;
    }
    /* メモリ開放１ */
    phase = -1;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

    delete[] ia;
    delete[] ja;
    delete[] a;
    delete[] b;

    return true;
}


/*//=======================================================
// ● MKL Pardisoで解く(複素/非対称)
//=======================================================*/
bool MatSolvers::solveMLKpardiso(const slv_int size0, const SparseMatC& matA, dcomplex* vecB, dcomplex* results, int num_para){
    const int n = size0;
    /* スパース形式コピー */
    int *ia = new int[n+1];
    auto row_ptr = matA.matrix->matrix.outerIndexPtr();
    const int total = matA.matrix->matrix.nonZeros();
    int* ja = new int[total];
    auto col_ptr = matA.matrix->getColPtr();
    dcomplex* a = new dcomplex[total];
    auto val_ptr = matA.matrix->getValuePtr();

    for(int i = 0; i < n; ++i){
        ia[i] = row_ptr[i];
        ia[i]++;
    }
    ia[size0] = total+1;
    for (int i = 0; i < total; ++i){
        ja[i] = col_ptr[i];
        ja[i]++;
        a[i] = val_ptr[i];
    }

    /* 実非対称設定 */
    int mtype = 13;       
    /* 右辺コピー */
    dcomplex* b = new dcomplex[n];
    for(int i = 0; i < n; ++i){
        b[i] = vecB[i];
    }

    int nrhs = 1; 
    void *pt[64];    
    int iparm[64];
    int maxfct, mnum, phase, error, msglvl;
    double ddum;    /* Double dummy */
    int idum;       /* Integer dummy. */

    omp_set_num_threads(num_para);
    for (int i = 0; i < 64; i++ ){
        iparm[i] = 0;
    }
    iparm[0] = 1;         /* No solver default */
    iparm[1] = 2;         /* Fill-in reordering from METIS */
    iparm[2] = num_para;
    iparm[3] = 0;         /* No iterative-direct algorithm */
    iparm[4] = 0;         /* No user fill-in reducing permutation */
    iparm[5] = 0;         /* Write solution into x */
    iparm[6] = 0;         /* Not in use */
    iparm[7] = 0;         /* Max numbers of iterative refinement steps */
    iparm[8] = 0;         /* Not in use */
    iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;        /* Not in use */
    iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
    iparm[13] = 1;        /* Output: Number of perturbed pivots */
    iparm[14] = 0;        /* Not in use */
    iparm[15] = 0;        /* Not in use */
    iparm[16] = 0;        /* Not in use */
    iparm[17] = 1;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = 1;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */
    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;             /* Which factorization to use. */
    msglvl = 1;           /* Print statistical information in file */
    error = 0;            /* Initialize error flag */    
    for (int i = 0; i < 64; i++ ){
        pt[i] = 0;
    }
    /* 処理１ */
    phase = 11;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 ){
        delete[] ia;
        delete[] ja;
        delete[] a;
        delete[] b;
        return false;
    }
    /* 処理２ */
    phase = 22;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 ){
        delete[] ia;
        delete[] ja;
        delete[] a;
        delete[] b;
        return false;
    }    
    /* 処理３ */
    phase = 33;
    iparm[7] = 2;   /* Max numbers of iterative refinement steps. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, results, &error);
    if ( error != 0 ){
        delete[] ia;
        delete[] ja;
        delete[] a;
        delete[] b;
        return false;
    }
    /* メモリ開放１ */
    phase = -1;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

    delete[] ia;
    delete[] ja;
    delete[] a;
    delete[] b;

    return true;
}


/* end of namespace */
};


#else

namespace SRLfem{

bool MatSolvers::solveMLKpardiso(const slv_int size0, const SparseMat& matA, double* vecB, double *results, int num_para){
    return false;
}
bool MatSolvers::solveMLKpardiso(const slv_int size0, const SparseMatC& matA, dcomplex* vecB, dcomplex* results, int num_para){
    return false;
}
};
#endif

