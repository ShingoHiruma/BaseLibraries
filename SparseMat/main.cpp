
#define _CRT_SECURE_NO_WARNINGS

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
using namespace std;

#include "SparseMatTMPL.hpp"
#include "SparseMatC.hpp"
#include "MatSolvers.hpp"

#include <chrono>

/*//=======================================================
// ●　メイン関数
//=======================================================*/
int main(int argc, char *argv[]){
	/*
	const int test_size = 1500;
	SRLfem::SparseMatBaseD mat1(test_size);
	SRLfem::SparseMatBaseD mat2(test_size);
	mat2.fix();

	auto t_start = std::chrono::system_clock::now();
	for(int i = 0; i < test_size; i++) {
		for(int j = 0; j < test_size; j++) {
			mat1.add(i, j, 1);
		}
	}
	mat1.fix();
	cout << "end1"<<endl;
	auto t_end = std::chrono::system_clock::now();

	for(int i = 0; i < test_size; i++) {
		for(int j = 0; j < test_size; j++) {
			mat2.add(i, j, 1);
		}
	}
	auto t_end2 = std::chrono::system_clock::now();
	cout << "end2"<<endl;

	auto t_dur = t_end - t_start;
	auto t_dur2 = t_end2 - t_end;
	auto t_sec = std::chrono::duration_cast<std::chrono::seconds>(t_dur).count();
	auto t_sec2 = std::chrono::duration_cast<std::chrono::seconds>(t_dur2).count();
	std::cout << "Original insert  time: " << t_sec << " sec \n";
	std::cout << "Eigen insert time: " << t_sec2 << " sec \n";
	*/

	double matA[5][5] ={
		{ 61,    0,   3.1,  -0.1,   1},
		{  0,   42.5, 2.0,   0.5,   4},
		{ 3.1,  2.0, 1.5,     0,   0},
		{-0.1, 0.5,   0,   6.4, 1.5},
		{   1,   4,   0,   1.5,  10}
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
			if( fabs(abs(matC[i][j])) > 1.0e-12 ){
				matCs1.add(i, j, matC[i][j]);
			}
		}
	}
	/* 疎行列位置を確定させる */
	matAs1.refresh();
	matBs1.refresh();
	matCs1.refresh();

	double diagD[100];
	SRLfem::SparseMat matLL = matAs1.IC_decomp(diagD, 1.0);
	//matLL.print();


	/*double vec[5] ={1,2,3,4,5};

	SRLfem::dcomplex* vec2 = matC*vec;
	for(int i = 0 ; i < 5 ; i++){
		cout << vec2[i] <<endl;
	}
	getchar();*/


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

	Eigen::VectorXd vecBe(total_size);
	for(int i = 0 ; i < total_size ; i++){
		vecBe(i) = vecB[i];
	}
	Eigen::VectorXd results_e;

	matAs1.print();
	cout << "---------" << endl;
	cout << "start...." << endl;

	bool bl = SRLfem::MatSolvers::solveEigenICCG(total_size, 1.0e-8, 10000, matAs1, vecBe, results_e);
	cout << "D norm3" << endl;
	Eigen::VectorXd norm1b = matAs1 * results_e;
	cout << norm1b <<endl;
	bool bl2x = SRLfem::MatSolvers::solveEigenICCG(total_size, 1.0e-8, 10000, matAs1, vecBe, results_e, false);


	//matAs1.add(4, 2, 11);
	/*vector<int> aa;
	vector<double> bb;
	matAs1.getTargetColVal(0, aa, bb);
	cout << "results, "  << endl;
	for(auto itr : aa) {
		cout << itr << endl;
	}
	matAs1.print();


	SRLfem::SparseMatBaseD matB2 = matAs1;
	matB2 *= 1.5;
	matB2.print();

	cout << &matAs1 << ", " << &matB2 << endl;*/

	return 1;

}
