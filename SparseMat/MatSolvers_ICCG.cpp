
#include "MatSolvers.hpp"
#include "SparseMatOperators.hpp"

/* 専用名前空間 */
namespace SRLfem{


/*

//=======================================================
//=======================================================
//=======================================================
ICCGソルバ
//=======================================================
//=======================================================

*/
/*//=======================================================
// ● ICCGで解く(右辺ノルム内部計算パターン)
//=======================================================*/
bool MatSolvers::solveICCG(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, const double *vecB, double *results){
	const slv_int size = size0;
	double norm = 0;
	for(int i = 0 ; i < size ; i++){
		norm += vecB[i]*vecB[i];
	}
	norm = sqrt(norm);
	return solveICCG(size, conv_cri, max_ite, accera, norm, matA, vecB, results);
}

/*//=======================================================
// ● ICCGで解く・外部実行本体
//=======================================================*/
bool MatSolvers::solveICCG(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const double normB, const SparseMat& matA, const double *vecB, double *results){
	/* コレスキー用スパース行列作成 */
	double *diagD = new double[size0];
	SparseMat matL = matA.IC_decomp(diagD, accera);
	SparseMat matL_tr = matL.trans();

	bool bl= solveICCG(size0, conv_cri, max_ite, accera, normB, diagD, *(matA.matrix), *(matL.matrix), *(matL_tr.matrix), vecB, results);
	delete[] diagD;
	return bl;	
}

/*//=======================================================
// ● ICCGで解く(入力右辺がEigen)
//=======================================================*/
bool MatSolvers::solveICCG(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, const Eigen::VectorXd& vecB, double *results){
	double* vecBa = new double[size0];
	double norm=0;
	for(int i = 0 ; i < size0 ; i++){
		vecBa[i] = vecB(i);
		norm += vecB[i]*vecB[i];
	}
	norm = sqrt(norm);
	bool bl = solveICCG(size0, conv_cri, max_ite, accera, norm, matA, vecBa, results);
	delete[] vecBa;
	return bl;
}

/*========================================*/
/*========================================*/
/*//=======================================================
// ● ICCGで解く（本体）
//=======================================================*/
bool MatSolvers::solveICCG(const slv_int size, const double conv_cri, const int max_ite, const double accera, const double normB, 
	const double* diagD, const SparseMatBaseD& matA, const SparseMatBaseD& matL, const SparseMatBaseD& matL_tr, const double *vecB, double *results){

	/* 要素確保 */
	double alpha;
	double beta;
	Eigen::VectorXd EvecP(size);
	Eigen::VectorXd EvecR(size);
	Eigen::VectorXd EvecLDV(size);
	Eigen::VectorXd EtempAP(size);


	/* 初期設定 */
	/*for(int i = 0 ; i < size ; i++){
		results[i] = 0;
	}*/
		
	slv_int* start_posA = new slv_int[size];
	slv_int* end_posA = new slv_int[size];
	matA.getCols(start_posA, end_posA);
	auto col_ptrA = matA.getColPtr();
	auto val_ptrA = matA.getValuePtr();
#ifdef OMP_USING_ICCG		
	#pragma omp parallel for
#endif
	for(slv_int ii = 0 ; ii < size ; ii++){
		const slv_int c_size = end_posA[ii];
		double ap_temp=0;
		for(slv_int j = start_posA[ii] ; j < c_size ; j++){
			ap_temp += val_ptrA[j] * results[col_ptrA[j]];
		}
		EtempAP(ii) = ap_temp;
		/* 初期残差等計算*/
		EvecR(ii) = vecB[ii] - ap_temp;
	}

	/* 前処理 */
	preProcess(size, matL, matL_tr, diagD, EvecR, EvecP);
	EvecLDV = EvecP;

	bool is_conv = false;
	/* 反復開始 */
	int It = 0;
	for(It = 0 ; It < max_ite ; It++){
		/* AP計算 */

//#ifdef OMP_USING_ICCG		
//		#pragma omp parallel for
//#endif
/*		for(slv_int ii = 0 ; ii < size ; ii++){
			const slv_int c_size = end_posA[ii];
			double tempVal = 0;
			for(slv_int j = start_posA[ii] ; j < c_size ; j++){
				tempVal	+= val_ptrA[j] * EvecP(col_ptrA[j]);
			}
			EtempAP(ii) = tempVal;
		}*/
		EtempAP = matA.matrix * EvecP;
		/* α計算 */
		const double temp = EvecR.dot(EvecLDV);
		const double temp2 = EvecP.dot(EtempAP);
		alpha = temp / temp2;

		/* 解ベクトルと残差計算 */
#ifdef OMP_USING_ICCG		
		#pragma omp parallel for
#endif
		for(slv_int i = 0 ; i < size ; i++){
			results[i] += alpha * EvecP(i);
		}
		EvecR = EvecR - alpha * EtempAP;
		const double normR = EvecR.norm() / normB;
		if(normR < conv_cri){
			//std::cout << "Solved!!! -- " << normR  << " " << It << std::endl;
			is_conv = true;
			break;
			//}else if( It % 100 == 0 ){
			//	std::cout << "NORMR is -- " << normR  << " " << It << std::endl;
		}
		/* v=(LDLT)-1rk　を計算 */
		preProcess(size, matL, matL_tr, diagD, EvecR, EvecLDV);
		/* β計算 */
		beta = EvecR.dot(EvecLDV);
		beta /= temp;
		/* P計算 */
		EvecP = beta*EvecP + EvecLDV;
	}
	delete[] start_posA;
	delete[] end_posA;
	if(It >= max_ite - 1){
		std::cout << "not Convergence!!! " << std::endl;
	}
	return is_conv ;
}

/*//=======================================================
// ● 前処理
//=======================================================*/
void MatSolvers::preProcess(const slv_int size0, const SparseMatBaseD& matL, const SparseMatBaseD& matL_tr, const double *diagD, const double *vecR, double *vec){
	const slv_int size = size0;
	double s;

	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();
	/* 第一方程式計算 */
	for (slv_int i = 0; i < size; i++) {
		s = vecR[i];
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec[j];
		}
		vec[i] = s / val_ptrL1[c_size];
	}
	/* 第二方程式計算 */
	for (slv_int i = size-1; i >= 0; i--) {	
		s = 0;
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s += val_ptrL2[j] * vec[the_row];
		}
		s *= diagD[i];
		vec[i] -= s;
	}
	delete[] start_posL1;
	delete[] end_posL1;
	delete[] start_posL2;
	delete[] end_posL2;
}

/*//=======================================================
// ● 前処理
//=======================================================*/
void MatSolvers::preProcess(const slv_int size0, const SparseMatBaseD& matL, const SparseMatBaseD& matL_tr, const double *diagD, const Eigen::VectorXd& EvecR, Eigen::VectorXd& vec){
	const slv_int size = size0;
	double s;
	slv_int* start_posL1 = new slv_int[size];
	slv_int* end_posL1 = new slv_int[size];
	matL.getCols(start_posL1, end_posL1);
	auto col_ptrL1 = matL.getColPtr();
	auto val_ptrL1 = matL.getValuePtr();

	slv_int* start_posL2 = new slv_int[size];
	slv_int* end_posL2 = new slv_int[size];
	matL_tr.getCols(start_posL2, end_posL2);
	auto col_ptrL2 = matL_tr.getColPtr();
	auto val_ptrL2 = matL_tr.getValuePtr();
	/* 第一方程式計算 */
	for (slv_int i = 0; i < size; i++) {
		s = EvecR(i);
		const slv_int c_size = end_posL1[i] - 1;
		for (slv_int jj = start_posL1[i]; jj < c_size; jj++) {
			slv_int j = col_ptrL1[jj];
			s -= val_ptrL1[jj] * vec(j);
		}
		vec(i) = s / val_ptrL1[c_size];
	}
	/* 第二方程式計算 */
	for (slv_int i = size-1; i >= 0; i--) {	
		s = 0;
		const slv_int c_size = end_posL2[i];
		for (slv_int j = start_posL2[i]+1 ; j < c_size; j++) {
			slv_int the_row = col_ptrL2[j];
			s += val_ptrL2[j] * vec(the_row);
		}
		s *= diagD[i];
		vec(i) -= s;
	}
	delete[] start_posL1;
	delete[] end_posL1;
	delete[] start_posL2;
	delete[] end_posL2;

#ifdef AAAAAAAAAA
	/* 第一方程式計算 */
	for (slv_int i = 0; i < size; i++) {
		s = EvecR(i);
		const slv_int c_size = matL.column_size[i] - 1;
		double *L_ptr = matL.matrix[i];
		slv_int *C_ptr = matL.column[i];
		for (slv_int jj = 0; jj < c_size; jj++) {
			slv_int j = *C_ptr ;
			s -= (*L_ptr) * vec(j);
			L_ptr++;
			C_ptr++;
		}
		vec(i) = s / (*L_ptr);
	}
	/* 第二方程式計算 */
	for (slv_int i = size-1; i >= 0; i--) {	
		s = 0;
		const slv_int c_size = matL_tr.column_size[i];
		double *L_ptr = matL_tr.matrix[i]+1;
		slv_int *C_ptr = matL_tr.column[i]+1;
		for (slv_int j = 1; j < c_size; j++) {
			slv_int the_row = *C_ptr ;
			s += (*L_ptr) * vec(the_row);
			L_ptr++;
			C_ptr++;
		}
		s *= diagD[i];
		vec(i) -= s;
	}
#endif
}

/* end of namespace */
};

