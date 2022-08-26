#ifndef DEF_MAT_SOLVER
#define DEF_MAT_SOLVER


//#define OMP_USING_MAT_SOL
//#define OMP_USING_ICCG
//#define INTEL_MKL_SOLVER_USING

#include "SparseMat.hpp"
#include "SparseMatC.hpp"
#include <cfloat>

#ifdef OMP_USING_MAT_SOL
#include <omp.h>
#elseifdef OMP_USING_ICCG
#include <omp.h>
#endif

/* Eigenを使ってみる */
#include <000_thirdparty/Eigen/Core>


/* オリジナル名前空間(静止器/回転機FEMライブラリ) */
namespace SRLfem{


/*
//=======================================================
// ■ スパース行列用ソルバ
//=======================================================*/
class MatSolvers{
private:
	/* ICCG専用内部処理 */
	static bool solveICCG(const slv_int size, const double conv_cri, const int max_ite, const double accera, const double normB, 
		const double* diagD, const SparseMatBaseD& matA, const SparseMatBaseD& matL, const SparseMatBaseD& matL_tr, const double *vecB, double *results, bool init=false);
	static void preProcess(const slv_int size0, const SparseMatBaseD& matL, const SparseMatBaseD& matL_tr, const double *diagD, const Eigen::VectorXd& EvecR, Eigen::VectorXd& vec);	
	static void preProcess(const slv_int size0, const SparseMatBaseD& matL, const SparseMatBaseD& matL_tr, const double *diagD, const double *vecR, double *vec);	
	/**/
	/* ICCOCG専用内部処理 */
	static bool solveICCG(const slv_int size, const double conv_cri, const int max_ite, const double accera, const double normB, 
						  const dcomplex* diagD, const SparseMatBaseC& matA, const SparseMatBaseC& matL, const SparseMatBaseC& matL_tr, const dcomplex *vecB, dcomplex *results, bool init=false);
	static void preProcess(const slv_int size0, const SparseMatBaseC& matL, const SparseMatBaseC& matL_tr, const dcomplex *diagD, const dcomplex *vecR, dcomplex *vec);	
	static void preProcess(const slv_int size0, const SparseMatBaseC& matL, const SparseMatBaseC& matL_tr, const dcomplex *diagD, const Eigen::VectorXcd& EvecR, Eigen::VectorXcd& vec);	
	/**/
	/* MKLソルバのラッパー */
	/**/
	template<typename MType, typename VType>
	static bool solveMLKpardisoBase(const slv_int size0, const MType& matA, VType* vecB, VType *results, int mat_mode, int num_para=1);
	/**/
	/*---------------*/
public:
	/*---------------*/
	/* ICCG法 */
	static bool solveICCG(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, const double *vecB, double *results, bool init=false);
	static bool solveICCG(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const double normB, const SparseMat& matA, const double *vecB, double *results, bool init=false);
	static bool solveICCG(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, const Eigen::VectorXd& vecB, double *results, bool init=false);
	/* ICCOCG法 */
	static bool solveICCG(slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMatC& matA, const dcomplex *vecB, dcomplex *results, bool init=false);
	static bool solveICCG(slv_int size0, const double conv_cri, const int max_ite, const double accera, const double normB, const SparseMatC& matA, const dcomplex *vecB, dcomplex *results, bool init=false);
	static bool solveICCG(slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMatC& matA, const Eigen::VectorXcd& vecB, dcomplex *results, bool init=false);
	//
	//
	/* Eigen内部ソルバのラッパーたち */
	//
	static bool solveEigenICCG(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, const Eigen::VectorXd& vecB, Eigen::VectorXd& results, bool init=false);
	static bool solveEigenICCG(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, const double* vecB, double* results, bool init=true);
	static bool solveEigenBiCGstab(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, const Eigen::VectorXd& vecB, Eigen::VectorXd& results, bool init=false);
	static bool solveEigenBiCGstab(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, const double* vecB, double* results, bool init=true);
	static bool solveEigenICCG(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results, bool init=false);
	static bool solveEigenICCG(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA, const dcomplex* vecB, dcomplex* results, bool init=true);
	static bool solveEigenBiCGstab(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results, bool init=false);
	static bool solveEigenBiCGstab(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA,  const dcomplex* vecB, dcomplex* results, bool init=false);
	//
	/* MKLソルバのラッパー */
	//
	static bool solveMLKpardisoSym(const slv_int size0, const SparseMat& matA, double* vecB, double *results, int num_para=1);
	static bool solveMLKpardisoSym(const slv_int size0, const SparseMatC& matA, dcomplex* vecB, dcomplex *results, int num_para=1);
	static bool solveMLKpardiso(const slv_int size0, const SparseMat& matA, double* vecB, double *results, int num_para=1);
	static bool solveMLKpardiso(const slv_int size0, const SparseMatC& matA, dcomplex* vecB, dcomplex *results, int num_para=1);
};

/* end of namespace */
};



#endif