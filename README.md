# BasicLibralies
## 電磁界有限要素法用基本ライブラリ
電磁界有限要素法で使いそうな基本関数と疎行列ライブラリです。
なお、内部で他のオープンソースEigenとnlohmann-jsonを利用しています（000_thirdparty内に配置）。

### BasicFuncsライブラリ
ガウスの積分点やC++11のメルセンヌ乱数のラッパーなどを提供する、基本関数ライブラリです

### SparseMat
疎行列とソルバのライブラリです。
Eigen疎行列のラッパー的に動作します。
ソルバは、オリジナル実装のICCG法と、Eigen内部のソルバのラッパーからなります。
