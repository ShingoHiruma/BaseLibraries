.. BaseLibraries documentation master file, created by sphinx-quickstart
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=====================================
BaseLibraries解説用ページ
=====================================

| このページは、電磁界有限要素法関連で使うC++用基本ライブラリの解説ページです。 
| ソースコードは下記にアップされています。 
https://github.com/JP-MARs/BaseLibraries

=====================================

ライブラリ内容
------------------


本ライブラリは複数のライブラリ群をまとめたもので、現在、以下のライブラリから成ります。

* **基本乱数・基本処理ライブラリ：BasicFuncs**　～メルセンヌ乱数やガウス積分用積分点の生成など、数値解析/最適化などで使いそうな基本処理をまとめたライブラリ
* **疎行列ライブラリ：SparseMat**　～有限要素法で用いる疎行列とソルバのライブラリ。基本的にはEigenライブラリ のラッパーとして機能します。実数と複素数の両方に対応します。

注釈

* 本ライブラリは全て名前空間 **「SRLfem」** 内に定義しています。呼び出す際は名前空間を付けるか、using namespaceを使ってください。
* 複素数の専用名として、SRLfem::dcomplexを定義しています。std::complex<double>として定義しているだけです。
* 行列のintはslv_intとして定義しています。intの別名として定義しているだけです。もしlong intなどに変えたければusingの項を変更してください。（詳細はSparseMatの説明のページで）

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Basic.md
   BaseLibraries.md
   SparseMat.md
   UserSample.md
   ForDev.md


