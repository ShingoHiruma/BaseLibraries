#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
using namespace std;


#include <SparseMat/SparseMat.hpp>
#include <SparseMat/SparseMatC.hpp>
#include <SparseMat/SparseMatOperators.hpp>
#include <SparseMat/MatSolvers.hpp>

#ifdef IS_WINDOWS_SISTEM
#ifdef _DEBUG
#pragma comment(lib, "libSparseMat_Deb.lib")
#else
#pragma comment(lib, "libSparseMat.lib")
#endif
#endif


/*//=======================================================
// ●　メイン関数
//=======================================================*/
int main(int argc, char *argv[]){

	/* 実数行列AとB */
	double matA[5][5] ={
		{ 60,    0,   3,  -0.1,   1},
		{  0,   45, 2.0,   0.5,   4},
		{  3,  2.0, 1.5,     0,   0},
		{-0.1, 0.5,   0,   6.4, 1.5},
		{   1,   4,   0,   1.5,  10}
	};
	double matB[5][5] ={
		{1, 2, 3, 0, 5},
		{2, 10, 0, 0, 0},
		{3, 0, 8, 0, 1},
		{0, 0, 0, 4, 0.5},
		{5, 0, 1, 0.5, 1}
	};
	/* 複素行列C */
	SRLfem::dcomplex matC[5][5] ={
		{15.5, 2, 3i, 0, 1.5},
		{2, 10, 0, 0, 0},
		{3i, 0, 8, 0, 1},
		{0, 0, 0, 4.8, 0.5i},
		{1.5, 0, 1, 0.5i, 20.5}
	};

	const int total_size = 5;

	/* 疎行列を定義 */
	SRLfem::SparseMat matAs1(total_size);
	SRLfem::SparseMat matBs1(total_size);
	SRLfem::SparseMatC matCs1(total_size);
	/* 行列A,B,Cの値をセット */
	for(int i = 0 ; i < 5 ; i++){
		for(int j = 0 ; j < 5 ; j++){
			if( fabs(matA[i][j]) > 1.0e-12 ){
				matAs1.add(i, j, matA[i][j]);
			}
			if( fabs(matB[i][j]) > 1.0e-12 ){
				matBs1.add(i, j, matB[i][j]);
			}
			if( fabs(abs(matC[i][j])) > 1.0e-12 ){
				matCs1.add(i, j, matC[i][j]);
			}
		}
	}
	/* 疎行列位置を確定させる */
	matAs1.fix();
	matBs1.fix();
	matCs1.fix();

	/**/
	/**/
	/**/
	/**/

	/* 適当なベクトル */
	double vecB[] = {1, 2, 3, 4, 5};
	SRLfem::dcomplex* vec0C = new SRLfem::dcomplex[total_size];
	for(int i = 0 ; i < total_size ; i++){
		vec0C[i] = i;//rand() / (1.0+RAND_MAX);
	}

	/* 行列・ベクトル積 */
	double* vec01 =matAs1 * vecB;
	for(int i = 0 ; i < total_size ; i++){
		cout << vec01[i] <<endl;
	}
	cout << endl;
	delete[] vec01;


	/* 行列の和や積 */
	SRLfem::SparseMat matX = matAs1 * 10.0;
	SRLfem::SparseMatC matCX = matAs1 * matCs1;
	SRLfem::SparseMatC matCX2 = matAs1 + matBs1 * matCs1;

	SRLfem::SparseMat matZ = SRLfem::SparseMatOperators::plusShift(matAs1, matBs1, 2.0, 0.5, 1, 2);
	matZ.print();



	/**/
	/**/
	/* ソルバ */
	/**/
	/**/

	/* 解を初期化 */
	double* results00 = new double[total_size];
	for(int i = 0 ; i < total_size ; i++){
		results00[i] = 0;
	}
	SRLfem::dcomplex* results01 = new SRLfem::dcomplex[total_size];
	for(int i = 0 ; i < total_size ; i++){
		results01[i] = 0;
	}

	cout <<"start" << endl;
	/* ICCG */
	bool bl1 = SRLfem::MatSolvers::solveICCG(total_size, 1.0e-8, 10000, 1.05, matAs1, vecB, results00);
	/* ICCOCG */
	bool bl2 = SRLfem::MatSolvers::solveICCG(total_size, 1.0e-8, 10000, 1.05, matCs1, vec0C, results01);


	cout << "D norm" << endl;
	double* norm1 = matAs1 * results00;
	for(int i = 0 ; i < total_size ; i++){
		cout << norm1[i] <<endl;
	}
	cout << "C norm" << endl;
	SRLfem::dcomplex* norm2 = matCs1 * results01;
	for(int i = 0 ; i < total_size ; i++){
		cout << norm2[i] <<endl;
	}

	delete[] norm1;
	delete[] norm2;




	return 1;

}
