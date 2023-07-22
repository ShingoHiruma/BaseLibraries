#ifndef DEF_SPR_MAT_TMPL
#define DEF_SPR_MAT_TMPL

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>


#include <BasicDefines.hpp>
#include <BasicDefinesC.hpp>
#include <vector>
#include <map>
#include <string>
#include <000_thirdparty/Eigen/SparseCore>

/* オリジナル名前空間(静止器/回転機FEMライブラリ) */
namespace SRLfem{


/* ソルバ内のint定義 */
using slv_int = int;

/* プロトタイプ宣言 */
class SparseMat;
class SparseMatC;
class SparseMatOperators;
class MatSolvers;

/*
//=======================================================
// ■ スパース行列テンプレート（自作形式）
//=======================================================*/
template<typename DType> class SparseMatTMPL{
/* フレンドクラス定義 */
friend class SparseMat;
friend class SparseMatC;
friend class SparseMatOperators;
friend class MatSolvers;
private:
	slv_int size;																/* 行列の行数 */
	bool is_fix;																/* スパース位置が確定しているか */
	std::map<slv_int, DType>* tempMat;											/* 一時保存行列（キーが列位置） */
	Eigen::SparseMatrix<DType, Eigen::RowMajor> matrix;											/* Eigenスパース行列 */
	/*-------------------------------------------------------*/
	void add_typeN(slv_int gyo, slv_int retu, DType val);							/* 確定済み行列にpush */
	/*  */
public:
	SparseMatTMPL();
	SparseMatTMPL(slv_int size0);																	/* コンストラクタ */
	SparseMatTMPL(const Eigen::SparseMatrix<DType, Eigen::RowMajor>& Mat0);
	SparseMatTMPL(Eigen::SparseMatrix<DType, Eigen::RowMajor>&& Mat0);
	SparseMatTMPL(const std::vector<Eigen::Triplet<DType>>& sparse_data);						/* Eigenタプルから初期化 */
	~SparseMatTMPL();
	slv_int getSize() const{return size;};
	SparseMatTMPL(const SparseMatTMPL<DType>& Mat);												/* コピーコンストラクタ */
	SparseMatTMPL(SparseMatTMPL<DType>&& Mat) noexcept;											/* ムーブコンストラクタ */
	SparseMatTMPL<DType>& operator=(const SparseMatTMPL<DType>& Mat);							/* 代入オペレータ */
	SparseMatTMPL<DType>& operator=(SparseMatTMPL<DType>&& Mat) noexcept;						/* 代入オペレータ（ムーブ） */
	void FixedCopy(const SparseMatTMPL<DType>& Mat);											/* Fix済みのMatコピー */
	bool isFixed() const{return is_fix;};														/* 確定済みかどうか */
	bool isEmpty() const;																		/* 行列の中身が空かどうか */
	void tempInitialize();																		/* 一時行列を作成 */
	void fixedInitialize(slv_int size0);														/* 空でSparse行列をfix */
	void fix();																					/* 一時配列を確定させる */
	void refresh(){fix();};
	void resetMat();																			/* 確定済み行列の値を０に再セット */
	void getCols(slv_int* cols_data1, slv_int* cols_data2) const;								/* 各行の列のスタート・終了位置を返す */
	auto getColPtr()const{ return matrix.innerIndexPtr(); };									/* 列のポインタを返す */
	auto getValuePtr()const{ return matrix.valuePtr(); };										/* 値のポインタを返す */
	slv_int isInclude(slv_int gyo, slv_int target_r) const;										/* i行目にtarget_r列があるかどうか（あったらそのindexを返す） */	
	slv_int getMaxCol()const;																	/* スパース内の最大の列位置を返す */
	void add(slv_int gyo, slv_int retu, DType val);												/* 一時配列にpush */
	void getTargetRowVal(slv_int target, std::vector<slv_int>& row_pos, std::vector<DType>& row_val)const;	/* 指定した列の非ゼロの行位置と値をvectorに書き出す */
	void getTargetColVal(slv_int target, std::vector<slv_int>& col_pos, std::vector<DType>& col_val)const;	/* 指定した行の非ゼロの列位置と値をvectorに書き出す */
	void operator*=(const DType xval);															/* 自身にスカラ積 */
	/*---------------------------------------------*/
	/*---------------------------------------------*/
	void printMat(const std::string& str="Mat.csv");
	void print() { std::cout << matrix << std::endl;};
};

/*=======================================================================*/
/*=======================================================================*/

/*//=======================================================
// ● コンストラクタ
//=======================================================*/
template<typename DType>
SparseMatTMPL<DType>::SparseMatTMPL(){
	size=0;is_fix=false;
	tempMat=nullptr;
}
template<typename DType>
SparseMatTMPL<DType>::SparseMatTMPL(slv_int size0){
	size=size0;
	is_fix=false;
	tempMat=nullptr;
	/* 一時行列初期化 */
	tempInitialize();
}
template<typename DType>
SparseMatTMPL<DType>::SparseMatTMPL(const Eigen::SparseMatrix<DType, Eigen::RowMajor>& Mat0){
	tempMat = nullptr;
	size = Mat0.rows();
	is_fix = true;
	matrix = Mat0;
}
template<typename DType>
SparseMatTMPL<DType>::SparseMatTMPL(Eigen::SparseMatrix<DType, Eigen::RowMajor>&& Mat0){
	tempMat = nullptr;
	size = Mat0.rows();
	is_fix = true;
	matrix = std::move(Mat0);
}

/*//=======================================================
// ● デストラクタ
//=======================================================*/
template<typename DType>
SparseMatTMPL<DType>::~SparseMatTMPL(){
	if(is_fix){
		matrix.data().squeeze(); 
		matrix.resize(1,1);
	}else{
		if(tempMat != nullptr){
			delete[] tempMat;
		}
	}	
}

/*//=======================================================
// ● コピーコンストラクタ
//=======================================================*/
template<typename DType>
SparseMatTMPL<DType>::SparseMatTMPL(const SparseMatTMPL<DType>& Mat){
	tempMat = nullptr;
	size = 0;
	this->size = Mat.size;
	*this = Mat;	
}
/*//=======================================================
// ● ムーブコンストラクタ
//=======================================================*/
template<typename DType>
SparseMatTMPL<DType>::SparseMatTMPL(SparseMatTMPL<DType>&& Mat) noexcept{
	tempMat = nullptr;
	this->size = Mat.size;
	*this = std::move(Mat);	
}

/*//=======================================================
// ● 代入オペレータ
//=======================================================*/
template<typename DType>
SparseMatTMPL<DType>& SparseMatTMPL<DType>::operator=(const SparseMatTMPL<DType>& Mat){
	bool same = (this->size == Mat.size);
	if( !same ){
		if(tempMat != nullptr){
			delete[] tempMat;
			tempMat = nullptr;
		}
	}
	/* 確定済み部分をコピーするだけ */
	if(Mat.is_fix){
		//this->matrix.swap(Mat.matrix);
		this->matrix =Mat.matrix;
		this->is_fix = Mat.is_fix;
		this->size = Mat.size;
	}else{
		/* コピー先・コピー元が未確定の時 */
		this->is_fix = Mat.is_fix;
		this->size = Mat.size;
		if(tempMat == nullptr){
			tempMat = new std::map<slv_int, DType>[Mat.size];
		}
		for(slv_int i = 0 ; i < size ; i++){
			this->tempMat[i].clear();
			for(auto itr : Mat.tempMat[i]){
				tempMat[i][itr.first] = itr.second;
			}
		}
	}
	return *this;
}

/*//=======================================================
// ● 代入オペレータ（ムーブ）
//=======================================================*/
template<typename DType>
SparseMatTMPL<DType>& SparseMatTMPL<DType>::operator=(SparseMatTMPL<DType>&& Mat) noexcept{
	/* 既存データは削除 */
	if(tempMat != nullptr){
		delete[] tempMat;
		tempMat = nullptr;
	}
	/* スパース位置が確定していたら、確定行列値を移動コピー*/
	if(Mat.is_fix){		
		this->is_fix = Mat.is_fix;
		this->size = Mat.size;
		this->matrix.swap(Mat.matrix);
	/* 未確定なら一時配列を移動コピー */
	}else{
		/* コピー先がまだ空なら全移動コピー */
		this->is_fix = Mat.is_fix;
		this->size = Mat.size;
		tempMat = Mat.tempMat;
		Mat.tempMat = nullptr;
	}
	return *this;
}
/*//=======================================================
// ● Fix済みのMatコピー
//=======================================================*/
template<typename DType>
void SparseMatTMPL<DType>::FixedCopy(const SparseMatTMPL<DType>& Mat){
	//this->matrix.swap(Mat.matrix);
	this->matrix = Mat.matrix;
	this->is_fix = Mat.is_fix;
	this->size = Mat.size;
}


/*//=======================================================
// ● 一時行列の初期化
//=======================================================*/
template<typename DType>
void SparseMatTMPL<DType>::tempInitialize(){
	/* スパース位置が確定していたら、行列値を初期化*/
	if(is_fix){
		matrix.setZero();
		return;
	}
	if(size == 0){
		std::cout << "Mat size is not initialized ! "<< std::endl;
		exit(1);
	}
	/* それ以外は、一時行列を初期化 */
	if(tempMat == nullptr){
		tempMat = new std::map<slv_int, DType>[size];
		for(slv_int i = 0 ; i < size ; i++){
			tempMat[i].clear();
		}
	}
}


/*//=======================================================
// ● 空でSparse行列をfix
//=======================================================*/
template<typename DType>
void SparseMatTMPL<DType>::fixedInitialize(slv_int size0){
	is_fix = true;
	size = size0;
	matrix = Eigen::SparseMatrix<DType>(size, size);
	if(tempMat != nullptr){
		delete[] tempMat;
		tempMat = nullptr;
	}
}

/*//=======================================================
// ● 一時配列を確定させる
//=======================================================*/
template<typename DType>
void SparseMatTMPL<DType>::fix(){
	/* スパース位置が確定していたらなにもしない*/
	if(is_fix){
		std::cout << "already fixed sparse!"<< std::endl;
		return;
	}
	/* Eigenに代入 */
	slv_int max_retu = -1;
	std::vector<Eigen::Triplet<DType>> tripletList;
	for(slv_int i = 0; i < size; i++) {
		const slv_int ss = tempMat[i].size();
		if(ss == 0) {
			continue;
		}
		for(auto& itr : tempMat[i]) {
			const slv_int pos = itr.first;
			const DType val = itr.second;
			tripletList.push_back(Eigen::Triplet<DType>(i, pos, val));
			if(pos > max_retu){
				max_retu = pos;
			}
		}
		std::map<slv_int, DType> empty; empty.clear();
		tempMat[i].swap(empty);
	}
	/* Eigenに確定(列方向の最大値も指定する必要あり)。行数の方が大きいなら正方、列の方が大きいなら長方形サイズにする */
	max_retu++;
	if(max_retu < size){
		max_retu = size;
	}
	matrix = Eigen::SparseMatrix<DType, Eigen::RowMajor>(size, max_retu);
	matrix.setFromTriplets(tripletList.begin(), tripletList.end());

	delete[] tempMat;
	tempMat = nullptr;
	is_fix = true;
}

/*//=======================================================
// ●  Eigenタプルから初期化
//=======================================================*/
template<typename DType>
SparseMatTMPL<DType>::SparseMatTMPL(const std::vector<Eigen::Triplet<DType>>& sparse_data){
	const int data_size = sparse_data.size();
	
	/* 最大列をサーチ */
	slv_int max1 = 0;
	slv_int max2 = 0;
	for(int i = 0 ; i < data_size ; i++){
		slv_int gyo = sparse_data[i].row();
		slv_int retu = sparse_data[i].col();
		if(max1 < gyo) max1 = gyo;
		if(max2 < retu) max2 = retu;
	}
	max1++;
	max2++;
	/* 確定させる */
	this->size = max1;
	matrix = Eigen::SparseMatrix<DType, Eigen::RowMajor>(max1, max2);
	matrix.setFromTriplets(sparse_data.begin(), sparse_data.end());
	tempMat = nullptr;
	is_fix = true;
}


/*//=======================================================
// ● 確定済み行列の値を０に再セット
//=======================================================*/
template<typename DType>
void SparseMatTMPL<DType>::resetMat(){
	/* 確定済みじゃなければ何もしない */
	if(!is_fix){
		return;
	}
	/* ポインタをたどってゼロにする */
	auto val_ptr = matrix.valuePtr();
	const slv_int total_size = matrix.nonZeros();
	slv_int count = 0;
	for(slv_int i = 0; i < total_size; i++) {
		val_ptr[i] = 0.0;
	}
}

/*//=======================================================
// ● 各行の列のスタート・終了位置を返す
//=======================================================*/
template<typename DType>
void SparseMatTMPL<DType>::getCols(slv_int* cols_data1, slv_int* cols_data2) const{
	auto col_ptr = matrix.outerIndexPtr();
	for(slv_int i = 0; i < size-1; i++) {
		cols_data1[i] = col_ptr[i];
		cols_data2[i] = col_ptr[i+1];
	}
	cols_data1[size-1] = col_ptr[size-1];
	cols_data2[size-1] = matrix.nonZeros();
}

/*//=======================================================
// ● 行列の中身が空かどうか(列数がすべての列でゼロなら空でTrue)
//=======================================================*/
template<typename DType>
bool SparseMatTMPL<DType>::isEmpty() const{
	const int nun_zero = matrix.nonZeros();
	return( nun_zero == 0 );
}

/*//=======================================================
// ● i行目にtarget_r列があるかどうか（あったらそのindexを返す）
//=======================================================*/
template<typename DType>
slv_int SparseMatTMPL<DType>::isInclude(slv_int gyo, slv_int target_r) const{
	if(is_fix){
		slv_int count=0;
		/* ポインタをたどって探す */
		auto row_ptr = matrix.innerIndexPtr();
		auto col_ptr = matrix.outerIndexPtr();
		auto val_ptr = matrix.valuePtr();
		const slv_int total_size = matrix.nonZeros();
		for(slv_int i = 0; i < size; i++) {
			const slv_int num = (i==size-1 ? total_size : col_ptr[i+1]);
			for(int j = col_ptr[i]; j < num; j++) {
				if(i == gyo && row_ptr[count] == target_r) {
					return j-col_ptr[i];
					//return count;
				}
				count++;
			}
		}
	}else{
		slv_int counter=0;
		for(const auto& itr : tempMat[gyo]){
			if( itr.first == target_r){
				return( counter );
			}
			counter++;
		}
	}
	return(-1);
}

/*//=======================================================
// ● mapにインサート
//=======================================================*/
template<typename DType>
void SparseMatTMPL<DType>::add(slv_int gyo, slv_int retu, DType val){
	/* スパース位置が確定していたら特別処理へ */
	if(is_fix){
		add_typeN(gyo, retu, val);
		return;
	}
	tempMat[gyo][retu] += val;
}

/*//=======================================================
// ● 行列にインサート
//=======================================================*/
template<typename DType>
void SparseMatTMPL<DType>::add_typeN(slv_int gyo, slv_int retu, DType val){
	matrix.coeffRef(gyo, retu) += val;
}

/*//=======================================================
// ● 指定した列の非ゼロの行位置と値をvectorに書き出す 
//=======================================================*/
template<typename DType>
void SparseMatTMPL<DType>::getTargetRowVal(slv_int target, std::vector<slv_int>& row_pos, std::vector<DType>& row_val) const{
	if( !is_fix ){
		return;
	}
	row_pos.clear();
	row_val.clear();
	slv_int count=0;
	/* ポインタをたどって探す */
	auto row_ptr = matrix.innerIndexPtr();
	auto col_ptr = matrix.outerIndexPtr();
	auto val_ptr = matrix.valuePtr();
	const slv_int total_size = matrix.nonZeros();
	for(slv_int i = 0; i < size; i++) {
		const slv_int num = (i==size-1 ? total_size : col_ptr[i+1]);
		for(int j = col_ptr[i]; j < num; j++) {
			if(row_ptr[count] == target) {
				row_pos.push_back(i);
				row_val.push_back(val_ptr[count]);
			}
			count++;
		}
	}
}
/*//=======================================================
// ● 指定した行の非ゼロの列位置と値をvectorに書き出す
//=======================================================*/
template<typename DType>
void SparseMatTMPL<DType>::getTargetColVal(slv_int target, std::vector<slv_int>& col_pos, std::vector<DType>& col_val)const{
	if( !is_fix ){
		return;
	}
	col_pos.clear();
	col_val.clear();

	
	/* ポインタをたどって探す */
	auto row_ptr = matrix.innerIndexPtr();
	auto col_ptr = matrix.outerIndexPtr();
	auto val_ptr = matrix.valuePtr();
	const slv_int total_size = matrix.nonZeros();
	slv_int i = target;
	slv_int count = col_ptr[i];
	const slv_int num = (i==size-1 ? total_size : col_ptr[i+1]);
	for(int j = col_ptr[i]; j < num; j++) {
		col_pos.push_back(row_ptr[count]);
		col_val.push_back(val_ptr[count]);
		count++;
	}
}

 
/*//=======================================================
// ● スパース内の最大の列位置を返す
//=======================================================*/
template<typename DType>
slv_int SparseMatTMPL<DType>::getMaxCol()const{
	/* 確定前のとき */
	if( !is_fix ){
		slv_int max = 0;
		for(slv_int i = 0; i < size; i++) {
			auto itr = tempMat[i].begin();
			const slv_int the_size = tempMat[i].size();
			for(slv_int j = 0 ; j < the_size ; j++){
				const slv_int tmp = itr->first;
				if(max < tmp){
					max = tmp;
				}
				itr++;
			}
		}
		return max;
	}
	/* 確定済みのとき */
	/* ポインタをたどって探す */
	slv_int count=0;
	auto row_ptr = matrix.innerIndexPtr();
	auto col_ptr = matrix.outerIndexPtr();
	auto val_ptr = matrix.valuePtr();
	const slv_int total_size = matrix.nonZeros();
	slv_int max = 0;
	for(slv_int i = 0; i < size; i++) {
		const slv_int num = (i==size-1 ? total_size : col_ptr[i+1]);
		for(int j = col_ptr[i]; j < num; j++) {
			if(max < row_ptr[count]) {
				max = row_ptr[count];
			}
			count++;
		}
	}
	return max;
}

/*//=======================================================
// ● 自身にスカラ積
//=======================================================*/
template<typename DType>
void SparseMatTMPL<DType>::operator*=(const DType xval){
	this->matrix *= xval;
}


/*//=======================================================
// ● 書き出し
//=======================================================*/
template<typename DType>
void SparseMatTMPL<DType>::printMat(const std::string& str){
	/*std::cout << "????????????????"<<std::endl;
	auto row_ptr = matrix.innerIndexPtr();
	auto col_ptr = matrix.outerIndexPtr();
	auto val_ptr = matrix.valuePtr();
	const slv_int total_size = matrix.nonZeros();
	for(slv_int i = 0 ; i <total_size; i++){
		std::cout << row_ptr[i] << ", " << col_ptr[i] << ", "<< val_ptr[i] << std::endl;
	}*/

	slv_int* start_pos1 = new slv_int[size];
	slv_int* end_pos1 = new slv_int[size];
	this->getCols(start_pos1, end_pos1);
	auto col_ptr1 = matrix.innerIndexPtr();;//matrix.outerIndexPtr();
	auto val_ptr1 = matrix.valuePtr();

	std::fstream fp(str, std::ios::out);
	fp << "Column size : " << size << std::endl;
	for(slv_int i = 0 ; i < size ; i++){
		fp << end_pos1[i]-start_pos1[i] <<  std::endl;
	}

	fp << "Column : "  << std::endl;
	for(slv_int i = 0 ; i < size ; i++){
		slv_int col_size = end_pos1[i]-start_pos1[i];
		if(col_size == 0){
			fp << "-1" << std::endl;
		}else{
			for(slv_int j = start_pos1[i] ; j < end_pos1[i] ; j++){
				fp << col_ptr1[j] << ", ";
			}
			fp << std::endl;
		}
	}
	fp <<std::endl << "Mat value : "  << std::endl;
	for(slv_int i = 0 ; i < size ; i++){
		slv_int col_size = end_pos1[i]-start_pos1[i];
		if(col_size > 0){
			for(slv_int j = start_pos1[i] ; j < end_pos1[i] ; j++){
				fp <<  std::scientific << val_ptr1[j] << ", ";
			}
			fp << std::endl;
		}else {
			fp << std::endl;
		}
	}
	fp.close();
	delete[] start_pos1;
	delete[] end_pos1;

}


/* テンプレート化 */
template class SparseMatTMPL<double>;
template class SparseMatTMPL<dcomplex>;

using SparseMatBaseD = SparseMatTMPL<double>;
using SparseMatBaseC = SparseMatTMPL<dcomplex>;

/* end of namespace */
};

#endif
