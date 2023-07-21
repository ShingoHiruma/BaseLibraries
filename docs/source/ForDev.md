# 開発者向け
本ライブラリの内部情報について少し記述しています。
改造案はウェルカムです

## SparseMatについて
### SparseMatの構成
疎行列SparseMatとSparseMatCは、中でSparseMatBaseD/SparseMatBaseCを疎行列の実体として持っています。
SparseMatとSparseMatCはその実体に対する操作を仲介しているだけで、実際のデータはSparseMatBaseD/SparseMatBaseCにあります。
SparseMatBaseD/SparseMatBaseCは、C++テンプレートを使ってEigenの疎行列のラッパを作っています。

* SparseMatTMPL：Eigen疎行列のラッパ。mapを使って仮行列データを作成し、fixが呼ばれるとEigenの疎行列化します。これはヘッダオンリで書いてます。
* SparseMatBaseD：SparseMatTMPL<double>のテンプレート化です。
```cpp:
　template class SparseMatTMPL<double>;
　using SparseMatBaseD = SparseMatTMPL<double>;
```
* SparseMatBaseC：SparseMatTMPL<SRLfem::dcomplex>のテンプレート化です。
```cpp:
　template class SparseMatTMPL<dcomplex>;
　using SparseMatBaseC = SparseMatTMPL<dcomplex>;
```

------------

### 疎行列の要素アクセスについて：

Eigen疎行列は、 SparseMatTMPLのメンバ「Eigen::SparseMatrix<DType, Eigen::RowMajor> matrix」になります。
見ての通り、RowMajorで作成しています。
疎行列の各要素にアクセスしたい場合はEigen疎行列のポインタ情報にアクセスする必要があります。
アクセス用のメソッドをSparseMatTMPLで定義しているので利用してください。

```cpp:

//疎行列の実体を以下とします（実際はSparseMatTMPLのメンバなので個別に作る必要はありません）
SparseMat matA;
//この行列の疎行列実体にアクセス(ポインタです)し、SparseMatTMPLのメソッドを呼び出します
slv_int size = matA.size;
slv_int* start_posL1 = new slv_int[size];
slv_int* end_posL1 = new slv_int[size];
//疎行列のRow方向ポインタを取得します
auto col_ptrL1 = matA.matrix->getColPtr();
//疎行列の値ポインタを取得します
auto val_ptrL1 = matA.matrix->getValuePtr();
//Row方向ポインタの各行の範囲を取得します
matA.matrix->getCols(start_posL1, end_posL1);

//取得した範囲を、各行ごとにループすれば参照できます
for(slv_int ii = 0 ; ii < size ; ii++){
    const slv_int c_size = end_posA[ii];
    //行iiにおける列方向の値があるところは以下のように
    for(slv_int j = start_posA[ii] ; j < c_size ; j++){
        slv_int retu = col_ptrA[j];
        double temp = val_ptrA[j];
        cout << ii <<", " << retu <<  " = " << temp << end;
    }
}
	
```


