
#include "SparseMatOperators.hpp"

#include "SparseMat.hpp"
#include "SparseMatC.hpp"

/* オリジナル名前空間(静止器/回転機FEMライブラリ) */
namespace SRLfem
{

	/*//=======================================================
	// ● 実数→複素化
	//=======================================================*/
	void SparseMatOperators::copyDtoC(SparseMatBaseC &mat_ans, const SparseMatBaseD &mat)
	{
		const slv_int the_size = mat.size;
		SparseMatBaseC *tempB = new SparseMatBaseC(the_size);

		slv_int *start_pos = new slv_int[the_size];
		slv_int *end_pos = new slv_int[the_size];
		mat.getCols(start_pos, end_pos);
		auto col_ptr = mat.getColPtr();	  // matrix.innerIndexPtr();.matrix.innerIndexPtr();
		auto val_ptr = mat.getValuePtr(); //.matrix.valuePtr();
		for (slv_int i = 0; i < the_size; i++)
		{
			const slv_int col_size = end_pos[i];
			for (slv_int j = start_pos[i]; j < col_size; j++)
			{
				const slv_int col = col_ptr[j];
				const double val = val_ptr[j];
				tempB->add(i, col, val);
			}
		}
		delete[] start_pos;
		delete[] end_pos;
		tempB->fix();
		mat_ans = std::move(*tempB);
		delete tempB;
	}

	/*//=======================================================
	// ● (a1*mat1 + a2*mat2)で、mat2に足す位置をpos1,pos2ずらすを計算
	//=======================================================*/
	SparseMat SparseMatOperators::plusShift(const SparseMat &mat1, const SparseMat &mat2, const double a1, const double a2, const slv_int pos1, const slv_int pos2)
	{
		SparseMatBaseD mat_ori;
		SparseMatOperators::plus_operators<SparseMatBaseD, SparseMatBaseD, SparseMatBaseD, double, double>(mat_ori, *(mat1.matrix), *(mat2.matrix), a1, a2, pos1, pos2);
		/* 結果を本体行列にムーブし、終わる */
		SparseMat mat_ans(std::move(mat_ori));
		return mat_ans;
	}
	SparseMatC SparseMatOperators::plusShift(const SparseMatC &mat1, const SparseMatC &mat2, const double a1, const double a2, const slv_int pos1, const slv_int pos2)
	{
		SparseMatBaseC mat_ori;
		SparseMatOperators::plus_operators<SparseMatBaseC, SparseMatBaseC, SparseMatBaseC, dcomplex, dcomplex>(mat_ori, *(mat1.matrix), *(mat2.matrix), a1, a2, pos1, pos2);
		/* 結果を本体行列にムーブし、終わる */
		SparseMatC mat_ans(std::move(mat_ori));
		return mat_ans;
	}

	/*//=======================================================
	// ● (mat1 + mat2)*vecBを計算
	//=======================================================*/
	double *SparseMatOperators::dotVecMat2(const SparseMat &mat1, const SparseMat &mat2, const double *vec1)
	{
		double *ans;
		SparseMatOperators::productVecMat2<SparseMatBaseD, SparseMatBaseD, double, double, double, double>(&ans, *(mat1.matrix), *(mat2.matrix), vec1);
		return ans;
	}
	/**/
	dcomplex *SparseMatOperators::dotVecMat2(const SparseMat &mat1, const SparseMatC &mat2, const double *vec1)
	{
		dcomplex *ans;
		SparseMatOperators::productVecMat2<SparseMatBaseD, SparseMatBaseC, double, dcomplex, double, dcomplex>(&ans, *(mat1.matrix), *(mat2.matrix), vec1);
		return ans;
	}
	/**/
	dcomplex *SparseMatOperators::dotVecMat2(const SparseMatC &mat1, const SparseMat &mat2, const double *vec1)
	{
		dcomplex *ans;
		SparseMatOperators::productVecMat2<SparseMatBaseC, SparseMatBaseD, dcomplex, double, double, dcomplex>(&ans, *(mat1.matrix), *(mat2.matrix), vec1);
		return ans;
	}
	/**/
	dcomplex *SparseMatOperators::dotVecMat2(const SparseMatC &mat1, const SparseMatC &mat2, const double *vec1)
	{
		dcomplex *ans;
		SparseMatOperators::productVecMat2<SparseMatBaseC, SparseMatBaseC, dcomplex, dcomplex, double, dcomplex>(&ans, *(mat1.matrix), *(mat2.matrix), vec1);
		return ans;
	}
	dcomplex *SparseMatOperators::dotVecMat2(const SparseMat &mat1, const SparseMat &mat2, const dcomplex *vec1)
	{
		dcomplex *ans;
		SparseMatOperators::productVecMat2<SparseMatBaseD, SparseMatBaseD, double, double, dcomplex, dcomplex>(&ans, *(mat1.matrix), *(mat2.matrix), vec1);
		return ans;
	}
	/**/
	dcomplex *SparseMatOperators::dotVecMat2(const SparseMat &mat1, const SparseMatC &mat2, const dcomplex *vec1)
	{
		dcomplex *ans;
		SparseMatOperators::productVecMat2<SparseMatBaseD, SparseMatBaseC, double, dcomplex, dcomplex, dcomplex>(&ans, *(mat1.matrix), *(mat2.matrix), vec1);
		return ans;
	}
	/**/
	dcomplex *SparseMatOperators::dotVecMat2(const SparseMatC &mat1, const SparseMat &mat2, const dcomplex *vec1)
	{
		dcomplex *ans;
		SparseMatOperators::productVecMat2<SparseMatBaseC, SparseMatBaseD, dcomplex, double, dcomplex, dcomplex>(&ans, *(mat1.matrix), *(mat2.matrix), vec1);
		return ans;
	}
	/**/
	dcomplex *SparseMatOperators::dotVecMat2(const SparseMatC &mat1, const SparseMatC &mat2, const dcomplex *vec1)
	{
		dcomplex *ans;
		SparseMatOperators::productVecMat2<SparseMatBaseC, SparseMatBaseC, dcomplex, dcomplex, dcomplex, dcomplex>(&ans, *(mat1.matrix), *(mat2.matrix), vec1);
		return ans;
	}

	/*//=======================================================
	// ● 行列AとBとCをかけて結果を返す
	//=======================================================*/
	SparseMat SparseMatOperators::dotMats(const SparseMat &matA, const SparseMat &matB, const SparseMat &matC)
	{
		SparseMatBaseD mat_ori;
		SparseMatOperators::MatProduct<SparseMatBaseD, SparseMatBaseD, SparseMatBaseD, SparseMatBaseD, double, double, double, double>(mat_ori, *(matA.matrix), *(matB.matrix), *(matC.matrix));
		/* 結果を本体行列にムーブし、終わる */
		SparseMat mat_ans(std::move(mat_ori));
		return mat_ans;
	}
	/**/
	SparseMatC SparseMatOperators::dotMats(const SparseMat &matA, const SparseMat &matB, const SparseMatC &matC)
	{
		SparseMatBaseC mat_ori;
		SparseMatOperators::MatProduct<SparseMatBaseD, SparseMatBaseD, SparseMatBaseC, SparseMatBaseC, double, double, dcomplex, dcomplex>(mat_ori, *(matA.matrix), *(matB.matrix), *(matC.matrix));
		/* 結果を本体行列にムーブし、終わる */
		SparseMatC mat_ans(std::move(mat_ori));
		return mat_ans;
	}
	/**/
	SparseMatC SparseMatOperators::dotMats(const SparseMat &matA, const SparseMatC &matB, const SparseMatC &matC)
	{
		SparseMatBaseC mat_ori;
		SparseMatOperators::MatProduct<SparseMatBaseD, SparseMatBaseC, SparseMatBaseC, SparseMatBaseC, double, dcomplex, dcomplex, dcomplex>(mat_ori, *(matA.matrix), *(matB.matrix), *(matC.matrix));
		/* 結果を本体行列にムーブし、終わる */
		SparseMatC mat_ans(std::move(mat_ori));
		return mat_ans;
	}
	/**/
	SparseMatC SparseMatOperators::dotMats(const SparseMatC &matA, const SparseMatC &matB, const SparseMatC &matC)
	{
		SparseMatBaseC mat_ori;
		SparseMatOperators::MatProduct<SparseMatBaseC, SparseMatBaseC, SparseMatBaseC, SparseMatBaseC, dcomplex, dcomplex, dcomplex, dcomplex>(mat_ori, *(matA.matrix), *(matB.matrix), *(matC.matrix));
		/* 結果を本体行列にムーブし、終わる */
		SparseMatC mat_ans(std::move(mat_ori));
		return mat_ans;
	}

	/*//=======================================================
	// ● (配列位置確定時)自身に、行列A±Bを加える(a1*A+a2*B)。Bの開始位置はposだけずらす
	//=======================================================*/
	void SparseMatOperators::plusFix(SparseMat &matAB, const SparseMat &matA, const SparseMat &matB, double a1, double a2, slv_int pos1, slv_int pos2)
	{
		SparseMatOperators::PlusMinusShiftFix<SparseMatBaseD, SparseMatBaseD, SparseMatBaseD, double, double, double>(*(matAB.matrix), *(matA.matrix), *(matB.matrix), a1, a2, pos1, pos2);
	}
	/**/
	void SparseMatOperators::plusFix(SparseMatC &matAB, const SparseMat &matA, const SparseMatC &matB, double a1, double a2, slv_int pos1, slv_int pos2)
	{
		SparseMatOperators::PlusMinusShiftFix<SparseMatBaseD, SparseMatBaseC, SparseMatBaseC, double, dcomplex, dcomplex>(*(matAB.matrix), *(matA.matrix), *(matB.matrix), a1, a2, pos1, pos2);
	}
	/**/
	void SparseMatOperators::plusFix(SparseMatC &matAB, const SparseMatC &matA, const SparseMatC &matB, double a1, double a2, slv_int pos1, slv_int pos2)
	{
		SparseMatOperators::PlusMinusShiftFix<SparseMatBaseC, SparseMatBaseC, SparseMatBaseC, dcomplex, dcomplex, dcomplex>(*(matAB.matrix), *(matA.matrix), *(matB.matrix), a1, a2, pos1, pos2);
	}

	/* (配列位置確定時)自身に、行列A±B±Cを加える(a1*A+a2*B+a3*C) */
	void SparseMatOperators::plusFix(SparseMat &matABC, const SparseMat &matA, const SparseMat &matB, const SparseMat &matC, double a1, double a2, double a3)
	{
		SparseMatOperators::PlusMinusShiftFix<SparseMatBaseD, SparseMatBaseD, SparseMatBaseD, SparseMatBaseD, double, double, double, double>(*(matABC.matrix), *(matA.matrix), *(matB.matrix), *(matC.matrix), a1, a2, a3);
	}
	void SparseMatOperators::plusFix(SparseMatC &matABC, const SparseMatC &matA, const SparseMatC &matB, const SparseMatC &matC, double a1, double a2, double a3)
	{
		SparseMatOperators::PlusMinusShiftFix<SparseMatBaseC, SparseMatBaseC, SparseMatBaseC, SparseMatBaseC, dcomplex, dcomplex, dcomplex, dcomplex>(*(matABC.matrix), *(matA.matrix), *(matB.matrix), *(matC.matrix), a1, a2, a3);
	}

	/*//=======================================================
	// ● (配列位置確定時)行列AとBをかけて、自身ABに格納する
	//=======================================================*/
	void SparseMatOperators::dotFix(SparseMat &matAB, const SparseMat &matA, const SparseMat &matB)
	{
		SparseMatOperators::MatProductFix<SparseMatBaseD, SparseMatBaseD, SparseMatBaseD, double, double, double>(*(matAB.matrix), *(matA.matrix), *(matB.matrix));
	}
	/**/
	void SparseMatOperators::dotFix(SparseMatC &matAB, const SparseMat &matA, const SparseMatC &matB)
	{
		SparseMatOperators::MatProductFix<SparseMatBaseD, SparseMatBaseC, SparseMatBaseC, double, dcomplex, dcomplex>(*(matAB.matrix), *(matA.matrix), *(matB.matrix));
	}
	/**/
	void SparseMatOperators::dotFix(SparseMatC &matAB, const SparseMatC &matA, const SparseMatC &matB)
	{
		SparseMatOperators::MatProductFix<SparseMatBaseC, SparseMatBaseC, SparseMatBaseC, dcomplex, dcomplex, dcomplex>(*(matAB.matrix), *(matA.matrix), *(matB.matrix));
	}

	/*//=======================================================
	// ● (配列位置確定時)行列AとBとCをかけて、自身ABCに格納する
	//=======================================================*/
	void SparseMatOperators::dotFix(SparseMat &matAB, const SparseMat &matA, const SparseMat &matB, const SparseMat &matC)
	{
		SparseMatOperators::MatProductFix<SparseMatBaseD, SparseMatBaseD, SparseMatBaseD, SparseMatBaseD, double, double, double, double>(*(matAB.matrix), *(matA.matrix), *(matB.matrix), *(matC.matrix));
	}
	/**/
	void SparseMatOperators::dotFix(SparseMatC &matAB, const SparseMat &matA, const SparseMat &matB, const SparseMatC &matC)
	{
		SparseMatOperators::MatProductFix<SparseMatBaseD, SparseMatBaseD, SparseMatBaseC, SparseMatBaseC, double, double, dcomplex, dcomplex>(*(matAB.matrix), *(matA.matrix), *(matB.matrix), *(matC.matrix));
	}
	/**/
	void SparseMatOperators::dotFix(SparseMatC &matAB, const SparseMat &matA, const SparseMatC &matB, const SparseMatC &matC)
	{
		SparseMatOperators::MatProductFix<SparseMatBaseD, SparseMatBaseC, SparseMatBaseC, SparseMatBaseC, double, dcomplex, dcomplex, dcomplex>(*(matAB.matrix), *(matA.matrix), *(matB.matrix), *(matC.matrix));
	}
	/**/
	void SparseMatOperators::dotFix(SparseMatC &matAB, const SparseMatC &matA, const SparseMatC &matB, const SparseMatC &matC)
	{
		SparseMatOperators::MatProductFix<SparseMatBaseC, SparseMatBaseC, SparseMatBaseC, SparseMatBaseC, dcomplex, dcomplex, dcomplex, dcomplex>(*(matAB.matrix), *(matA.matrix), *(matB.matrix), *(matC.matrix));
	}

	/* end of namespace */
};
