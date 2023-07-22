# SparseMat
有限要素法で使う疎行列用ライブラリです。
本ライブラリのstaticライブラリファイルは **「SparseMat.lib」** です。

EigenのEigen::SparseMatrixをベースとして作成されており、Eigenの疎行列ラッパーとして機能します。
double型、complex型の両方に対応し、疎行列の四則演算なども用意されています。
ソルバとしては、オリジナルの加速係数つきICCGソルバが用意されています。
その他、Eigenが用意しているICCGとBiCGstabのラッパも用意しています。

内部のintは、slv_intとして別の型として定義しています。

・using slv_int = int;　(SparseMatTMPL.hppの冒頭)

基本はこのままでいいと思いますが、超大規模にしたい場合はunsingedやlong intに書き換えてください

-------

-------

## SparseMatクラス
実数用の疎行列クラスです。
基本的には、まず初期化した後、add関数で行列要素を追加していきます。
追加が完了したら、fixメソッドで疎行列を確定させてください。
確定後にメソッドや四則演算が行えます。fix前に各種操作を行うとエラーになる恐れがあるので注意してください。
「SparseMat/SparseMat.hpp」をインクルードして使ってください。

### * 基本操作
-------

###### SparseMat();
コンストラクタ：疎行列を初期化します。インスタンスが作成されるだけで内部は何も初期化されません。
この状態でaddをしてもエラーになります。

###### SparseMat(slv_int x);
コンストラクタ：疎行列を初期化します。行数サイズをxで初期化します。
（内部でtempInitializeを呼び、初期化します）

###### SparseMat(const SparseMat& mat);
コンストラクタ：別の疎行列を代入して初期化します。本方法で初期化すると最初からfix状態になります。


###### bool isFixed()
呼び出し元インスタンスが確定済み(fixを呼び出した後か)かどうかを返します。

###### bool isEmpty() const
行列の中身が空かどうかを返します。

###### void tempInitialize()
疎行列の状況を初期化し、fixされていない状態にします。

###### void fix()
疎行列を確定させます。

###### void resetMat()
fix済みの疎行列の０にします。
（疎行列の位置はそのまま、値だけゼロにします）

###### void add(slv_int gyo, slv_int retu, double val)
疎行列の指定の位置に値を加えます。行の位置をgyoで、列の位置をretuで指定し、そこにValを加えます。
fix前の場合は任意の位置にaddできます。
fix後は、すでに値のある位置以外を指定するとエラーになります。

###### slv_int isInclude(slv_int gyo, slv_int target_r)const
疎行列のi行目に、target_r列があるかどうかを判定します。あったらそのindexを返します。

###### void getTargetRowVal(slv_int target, std::vector<slv_int>& row_pos, std::vector<double>& row_val)const
指定した列の非ゼロの行位置と値をvector「row_val」にコピーします。

###### void getTargetColVal(slv_int target, std::vector<slv_int>& col_pos, std::vector<double>& col_val)const
指定した行の非ゼロの列位置と値をvector「col_val」にコピーします。

###### void printMat(const std::string& str="Mat.csv")
疎行列の内容をstrで指定したファイルに書き出します。

###### void print()
疎行列をコンソールに表示します。

-------
### * オペレータ群
-------

###### double* operator*(const double* vec) const;
###### dcomplex* operator*(const dcomplex* vec) const;
疎行列に列ベクトルvecを掛けて、結果ベクトルを返します。
列ベクトルはC++のデフォルトの配列で渡します。
（vecは実数、複素数どちらでも対応）

###### Eigen::VectorXd operator*(const Eigen::VectorXd& vec) const;	
###### Eigen::VectorXcd operator*(const Eigen::VectorXcd& vec) const;	
疎行列に列ベクトルvecを掛けて、結果ベクトルを返します。
列ベクトルはEigenのVectorで渡します。
（vecは実数、複素数どちらでも対応）

###### void operator*=(const double x)
自身の全要素にxを掛けます。

###### SparseMat operator*(const double x) const;
###### SparseMatC operator*(const dcomplex x) const;
自身の全要素にxを掛けた後の疎行列を返します。
（実数、複素数どちらでも対応）

###### SparseMat operator*(const SparseMat& mat) const;
###### SparseMatC operator*(const SparseMatC& mat) const;
自身と行列matの行列積を行い、結果の疎行列を返します。

###### SparseMat operator+(const SparseMat& mat) const;
###### SparseMatC operator+(const SparseMatC& mat) const;
自身と行列matを足し、結果の疎行列を返します。

-------
### * その他の行列オペレータ群
-------

###### SparseMat trans() const
自身の転置行列を返します。

###### SparseMat makeSubMat(slv_int range1a, slv_int range1b, slv_int range2a, slv_int range2b)
範囲[r1a, r1b]×[r2a, r2b]の部分を抜き出して部分行列を作り、得られた部分行列を返します。

###### void MatDiv(SparseMat& matK11, SparseMat& matK12, slv_int rangeA, slv_int rangeB)
列rangeBを境に、行列を2つに分け、行rangeAより下は削除します。

###### SparseMat getMatLower() const;
下三角行列を作成して返します。

###### SparseMat getMatUpper() const
上三角行列を作成して返します。

###### SparseMat inv() const
逆行列（密アルゴリズム・大行列でやるとメモリが吹き飛ぶ！）

###### void round()
行列の要素を四捨五入して丸めます。

-------

## SparseMatCクラス
複素数用の疎行列クラスです。
基本的な操作は全て実数用のSparseMatと同じです。
「SparseMat/SparseMatC.hpp」をインクルードして使ってください。


-------

## SparseMatOperatorsクラス
複数の疎行列の操作を行うstaticなメソッドを管理するクラスです。
「SparseMat/SparseMatOperators.hpp」をインクルードして使ってください。

###### static SparseMat plusShift(const SparseMat& mat1, const SparseMat& mat2, const double a1, const double a2, const slv_int pos1, const slv_int pos2)
###### static SparseMatC plusShift(const SparseMatC& mat1, const SparseMatC& mat2, const double a1, const double a2, const slv_int pos1, const slv_int pos2)
(a1*mat1 + a2*mat2)を計算して、結果の疎行列を返します。
mat2に足す位置をpos1,pos2でずらすことが可能です。

###### static double* dotVecMat2(const SparseMat& mat1, const SparseMat& mat2, const double* vec1);
###### static dcomplex* dotVecMat2(const SparseMat& mat1, const SparseMatC& mat2, const double* vec1);
###### static dcomplex* dotVecMat2(const SparseMatC& mat1, const SparseMat& mat2, const double* vec1);
###### static dcomplex* dotVecMat2(const SparseMatC& mat1, const SparseMatC& mat2, const double* vec1);
###### static dcomplex* dotVecMat2(const SparseMat& mat1, const SparseMat& mat2, const dcomplex* vec1);
###### static dcomplex* dotVecMat2(const SparseMat& mat1, const SparseMatC& mat2, const dcomplex* vec1);
###### static dcomplex* dotVecMat2(const SparseMatC& mat1, const SparseMat& mat2, const dcomplex* vec1);
###### static dcomplex* dotVecMat2(const SparseMatC& mat1, const SparseMatC& mat2, const dcomplex* vec1);
(mat1 + mat2)*vecBを計算し、得られたベクトルを返します。

###### static SparseMat dotMats(const SparseMat& matA, const SparseMat& matB, const SparseMat& matC);
###### static SparseMatC dotMats(const SparseMat& matA, const SparseMat& matB, const SparseMatC& matC);
###### static SparseMatC dotMats(const SparseMat& matA, const SparseMatC& matB, const SparseMatC& matC);
###### static SparseMatC dotMats(const SparseMatC& matA, const SparseMatC& matB, const SparseMatC& matC);
疎行列AとBとCをかけて結果を返します。

###### static void plusFix(SparseMat& matAB, const SparseMat& matA, const SparseMat& matB, double a1=1.0, double a2=1.0, slv_int pos1=0, slv_int pos2=0);
###### static void plusFix(SparseMatC& matAB, const SparseMat& matA, const SparseMatC& matB, double a1=1.0, double a2=1.0, slv_int pos1=0, slv_int pos2=0);
###### static void plusFix(SparseMatC& matAB, const SparseMatC& matA, const SparseMatC& matB, double a1=1.0, double a2=1.0, slv_int pos1=0, slv_int pos2=0);
fix済みの行列matABCに、(a1*matA+a2*matB)の結果を代入します。
matBの開始位置はposだけずらすことが可能です。

###### static void plusFix(SparseMat& matABC, const SparseMat& matA, const SparseMat& matB, const SparseMat& matC, double a1=1.0, double a2=1.0, double a3=1.0);
###### static void plusFix(SparseMatC& matABC, const SparseMatC& matA, const SparseMatC& matB, const SparseMatC& matC, double a1=1.0, double a2=1.0, double a3=1.0);
fix済みの行列matABCに、(a1*A+a2*B+a3*C)の結果を代入します。

###### static void dotFix(SparseMat& matAB, const SparseMat& matA, const SparseMat& matB);
###### static void dotFix(SparseMatC& matAB, const SparseMat& matA, const SparseMatC& matB);
###### static void dotFix(SparseMatC& matAB, const SparseMatC& matA, const SparseMatC& matB);
fix済みの行列matABCに、AとBの行列積の結果を代入します。

###### static void dotFix(SparseMat& matAB, const SparseMat& matA, const SparseMat& matB, const SparseMat& matC);
###### static void dotFix(SparseMatC& matAB, const SparseMat& matA, const SparseMat& matB, const SparseMatC& matC);
###### static void dotFix(SparseMatC& matAB, const SparseMat& matA, const SparseMatC& matB, const SparseMatC& matC);
###### static void dotFix(SparseMatC& matAB, const SparseMatC& matA, const SparseMatC& matB, const SparseMatC& matC);
fix済みの行列matABCに、AとBとCの行列積の結果を代入します。


-------

## MatSolversクラス
疎行列用のソルバを提供するクラスです。全てstaticメソッドです。
「SparseMat/MatSolvers.hpp」をインクルードして使ってください。

###### static SparseMat plusShift(const SparseMat& mat1, const SparseMat& mat2, const double a1, const double a2, const slv_int pos1, const slv_int pos2)
###### static SparseMatC plusShift(const SparseMatC& mat1, const SparseMatC& mat2, const double a1, const double a2, const slv_int pos1, const slv_int pos2)
(a1*mat1 + a2*mat2)を計算して、結果の疎行列を返します。
mat2に足す位置をpos1,pos2でずらすことが可能です。

###### static bool solveICCG(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, const double *vecB, double *results, bool init=false)
###### static bool solveICCG(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const double normB, const SparseMat& matA, const double *vecB, double *results, bool init=false)
###### static bool solveICCG(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMat& matA, const Eigen::VectorXd& vecB, double *results, bool init=false)
自作のICCG法です。右辺をC++配列、Eigen::VectorXdで与えるかで2パターンあります。
また、引く数normBは、右辺bのノルムを事前に与えるか、内部で計算するかです。
（normBは、収束判定で使います。"r/normB < conv_cri"となると収束とします。）

size0：行列の行数です。
conv_cri：収束判定値です。10^-6など。
max_ite：最大反復数です。
accera：加速係数です。
matA：Ax=bのAです。対称行列としてください。
vecB:Ax=bのbです。size0のベクトルです。
results：解ベクトルです（Ax=bのx）。メモリは事前にsize0だけ確保しておいてください。
init：解ベクトルxをこの処理の前にゼロ初期化するかです。過去の値がxに入っていて流用したい場合はfalseにしてください。それ以外はtrueにすれば内部でゼロ初期化してから始めます。

###### static bool solveICCG(slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMatC& matA, const dcomplex *vecB, dcomplex *results, bool init=false);
###### static bool solveICCG(slv_int size0, const double conv_cri, const int max_ite, const double accera, const double normB, const SparseMatC& matA, const dcomplex *vecB, dcomplex *results, bool init=false);
###### static bool solveICCG(slv_int size0, const double conv_cri, const int max_ite, const double accera, const SparseMatC& matA, const Eigen::VectorXcd& vecB, dcomplex *results, bool init=false);
上記ICCG法の複素数版です。使い方は同じです。
ただしmatAは対称行列にのみ対応しています（エルミートなどは未対応）

###### static bool solveEigenICCG(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, const Eigen::VectorXd& vecB, Eigen::VectorXd& results, bool init=false);
###### static bool solveEigenICCG(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, const double* vecB, double* results, bool init=true);
###### static bool solveEigenBiCGstab(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, const Eigen::VectorXd& vecB, Eigen::VectorXd& results, bool init=false);
###### static bool solveEigenBiCGstab(const slv_int size0, const double conv_cri, const int max_ite, const SparseMat& matA, const double* vecB, double* results, bool init=true);
###### static bool solveEigenICCG(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results, bool init=false);
###### static bool solveEigenICCG(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA, const dcomplex* vecB, dcomplex* results, bool init=true);
###### static bool solveEigenBiCGstab(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA, const Eigen::VectorXcd& vecB, Eigen::VectorXcd& results, bool init=false);
###### static bool solveEigenBiCGstab(const slv_int size0, const double conv_cri, const int max_ite, const SparseMatC& matA,  const dcomplex* vecB, dcomplex* results, bool init=false);
Eigenに付属している疎行列ソルバのラッパーです。
使い方は上記ICCG法の引数と同じです。
