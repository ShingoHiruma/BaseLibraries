#ifndef DEF_MAT_SOLVER_ICCG_MYDEF
#define DEF_MAT_SOLVER_ICCG_MYDEF


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

/* オリジナル名前空間(静止器/回転機FEMライブラリ) */
namespace SRLfem{


/*
//=======================================================
// ■ スパース行列用ソルバ
//=======================================================*/
class MatSolversICCG{
private:
	/* 設定 */
	bool is_diag_scale;								/* 対角スケーリングを行うか */
	bool is_save_best;								/* 最良結果を保存して失敗時に利用するか */
	bool is_save_residual_log;						/* 残差履歴を残すか */
	std::vector<double> residual_log;				/* 残差履歴 */
	/* 発散判定：
	  0=最大反復までやる
	  1：最良×bad_div_valより大きい値がbad_div_countだけ続いたら発散として終わる */
	int diverge_judge_type;
	double bad_div_val;
	int bad_div_count_thres;
	/**/
	/**/
	/* ICCG専用内部処理 */
	bool solveICCG(const slv_int size, const double conv_cri, const int max_ite, const double normB, 
		const double* diagD, const SparseMatBaseD& matA, const SparseMatBaseD& matL, const SparseMatBaseD& matL_tr, const double *vecB, double *results, bool init=false);
	void preProcess(const slv_int size0, const SparseMatBaseD& matL, const SparseMatBaseD& matL_tr, const double *diagD, const Eigen::VectorXd& EvecR, Eigen::VectorXd& vec);	
	void preProcess(const slv_int size0, const SparseMatBaseD& matL, const SparseMatBaseD& matL_tr, const double *diagD, const double *vecR, double *vec);	
	/**/
	/* ICCOCG専用内部処理 */
	bool solveICCG(const slv_int size, const double conv_cri, const int max_ite, const double normB, 
						  const dcomplex* diagD, const SparseMatBaseC& matA, const SparseMatBaseC& matL, const SparseMatBaseC& matL_tr, const dcomplex *vecB, dcomplex *results, bool init=false);
	void preProcess(const slv_int size0, const SparseMatBaseC& matL, const SparseMatBaseC& matL_tr, const dcomplex *diagD, const dcomplex *vecR, dcomplex *vec);	
	void preProcess(const slv_int size0, const SparseMatBaseC& matL, const SparseMatBaseC& matL_tr, const dcomplex *diagD, const Eigen::VectorXcd& EvecR, Eigen::VectorXcd& vec);	
	/**/
	/*---------------*/
public:
	MatSolversICCG();		/* コンストラクタ */

	void setDiagScale(bool bl){ is_diag_scale = bl; };
	void setSaveBest(bool bl){ is_save_best = bl; };
	void setSaveLog(bool bl){ is_save_residual_log = bl; };
	void getResidualLog(std::vector<double>& log);
	void setDirvegeType(int x){ diverge_judge_type = x; };
	void setBadDivVal(double x){bad_div_val=x;};
	void setBadDivCount(int x){bad_div_count_thres=x;};
	/*---------------*/
	/*---------------*/
	/* ICCG法 */
	bool solveICCG(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, const double *vecB, double *results, bool init=false);
	bool solveICCG(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const double normB, const SparseMat& matA, const double *vecB, double *results, bool init=false);
	bool solveICCG(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, const Eigen::VectorXd& vecB, double *results, bool init=false);
	/* ICCOCG法 */
	bool solveICCG(slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMatC& matA, const dcomplex *vecB, dcomplex *results, bool init=false);
	bool solveICCG(slv_int size0, const double conv_cri, const int max_ite, const double accera, const double normB, const SparseMatC& matA, const dcomplex *vecB, dcomplex *results, bool init=false);
	bool solveICCG(slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMatC& matA, const Eigen::VectorXcd& vecB, dcomplex *results, bool init=false);
};

/* end of namespace */
};

#endif
