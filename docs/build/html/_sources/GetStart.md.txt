# [基本的な使い方]

## Get Start
利用方法：CMakeの簡易サンプルを置いています。コンパイルしてexampleを利用するだけなら、cmakeするだけです。

```

mkdir build
cd build
cmake ..

```

これでbuildフォルダにexampleの実行ファイルができるので、それを実行すればテストできます。
また、本ライブラリの元フォルダ直下に、Staticライブラリ（.libファイル）が出来ますので、自分のプログラムにリンクして使ってください。

## 自作プログラムへのリンク
cmakeでbuildしたlibファイルは、通常のstaticライブラリと同様の方法で自作プログラムにリンクできます。
ライブラリを置いた場所を「/BaseLibraries」とし、自作プログラムを「/TestProg」に置いたとします。
このときに、例えば本ライブラリの疎行列ライブラリを使う場合は、makefileに以下を記述します。

```

CXX = icpx
CFLAGS = -I/BaseLibraries
OPTS = -qmkl -ipo
LIBS = -L/BaseLibraries -lSparseMat

```


