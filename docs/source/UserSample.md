# 利用例
本ライブラリの利用方法の例を示します。
githubのexampleにも利用例があるので参照ください。

## BasicFuncsの利用例
### Mt乱数の利用例
Mt.hppで定義されているメルセンヌ乱数の利用例です。
まずはシードを初期化します。

```cpp:

#include <BasicFuncs/Mt.hpp>

int main(int argc, char *argv[]){
    //10のシードを作成
    int seed = 10;
    //初期化する（一応C++のデフォルト乱数も初期化）
    srand(seed);
    SRLfem::Mt::init_rand(seed);
}

```

あとはMt.hpp内に定義している乱数メソッドを自由に呼び出して利用してください。
例えば正規乱数なら以下のようになります。

```cpp:

/* 平均u、分散sig*sigの正規乱数の生成 */
double u = 5.0;
double sig = 100;
double rand　= SRLfem::Mt::normal_rand(u, sig);
cout << rand << endl;

```

------

### 基本関数の利用例
CommonFuncs.hppで定義されている基本関数の利用例です。
基本的には、やりたい処理をhppの定義から探して呼び出すだけです。
例えば、ガウス積分点[-1,+1]を生成したいなら以下のようになります。

```cpp:

int num_g = 6;
double* tg = new double[num_g];
double* wg = new double[num_g];

//ガウス積分点[-1, 1]の生成
SRLfem::CommonFuncs::setGauss(tg, wg, num_g);
for(int i = 0 ; i < num_g ; i++){
    cout << i << "th point = " << tg[i] << " and it's w = " << wg[i] << endl;
}

delete[] tg;
delete[] wg;

```


------

## SparseMatの利用例
### 疎行列の作成
SparseMat.hppもしくはSparseMatC.hppで定義している疎行列の利用例です。
基本的には、行数を決めてインスタンスを作成し、addで要素を付け足していくだけです。
行列が完成したらfixで固定します。


```cpp:

#include <SparseMat/SparseMat.hpp>
#include <SparseMat/SparseMatC.hpp>

int main(int argc, char *argv[]){
    int size = 5;
    /* 行数5の疎行列を作成して初期化 */
    SRLfem::SparseMat matA1(size);

    //適当に要素を付け足していく
    matA1.add(0, 0, 5.0);
    matA1.add(1, 2, -1.0);
    matA1.add(4, 5, 1.0);
    //matA1.add(5, 1, 1.0); //サイズが5をオーバーしたのでエラーで落ちます

    //同じ位置にaddしてもOKです
    matA1.add(0, 0, 1.0);

    //終わったら確定します
    matA1.fix();

    matA1.print();
    //(0, 0)=6
    //(1, 2)=-1
    //(4, 5)=1
    //と、addしたところだけ値を持った疎行列が表示されます

}

```

------


### 疎行列の四則演算
fix済みの疎行列は四則演算が可能です

```cpp:

#include <SparseMat/SparseMat.hpp>
#include <SparseMat/SparseMatC.hpp>
#include <SparseMat/SparseMatOperators.hpp>

int main(int argc, char *argv[]){
    //
    //～適当に行列を作成してfixさせます
    //  matA, matBを作っているとしましょう
    //

    /* 適当なベクトル */
    double vecB[] = {1, 2, 3, 4, 5};
    /* 行列・ベクトル積 */
    double* vec01 =matA * vecB;
    for(int i = 0 ; i < total_size ; i++){
        cout << vec01[i] <<endl;
    }
    cout << endl;
    delete[] vec01;//ベクトル積で作った結果は動的ポインタで保持しているので、使い終わったらdeleteして下さい

    /* 行列の和や積も自由にできます */
    SRLfem::SparseMat matX = matA * 10.0;
    SRLfem::SparseMat matCX = matA * matB;
    //複数の操作をまとめて書いても動きます
    SRLfem::SparseMat matCX2 = matA + matB * matCX;

    //SparseMatOperatorsで定義されている少し複雑な処理は、個別にstaticメソッドを読んで実行します
    SRLfem::SparseMat matZ = SRLfem::SparseMatOperators::plusShift(matAs1, matBs1, 2.0, 0.5, 1, 2);
    matZ.print();

}

```



------


### 疎行列ソルバ
fix済みの疎行列に対して、疎行列用ソルバを利用できます。



```cpp:

#include <SparseMat/SparseMat.hpp>
#include <SparseMat/SparseMatC.hpp>
#include <SparseMat/MatSolvers.hpp>

int main(int argc, char *argv[]){
    int total_siz = 5:

    //
    //～適当に行列を作成してfixさせます
    //  matAを作っているとしましょう
    //


	/* 解を初期化 */
    double* results00 = new double[total_size];
    for(int i = 0 ; i < total_size ; i++){
        results00[i] = 0;
    }    
    //右辺b
    double vecB[] = {1, 2, 3, 4, 5};
    /* ICCGでAx=bを解きます（収束判定10^-8、加速係数1.05） */
    bool bl1 = SRLfem::MatSolvers::solveICCG(total_size, 1.0e-8, 10000, 1.05, matAs1, vecB, results00);;


}

```

