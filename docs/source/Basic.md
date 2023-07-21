# 基本ヘッダ
本ライブラリで利用する基本ヘッダは、下記の3つです。
* BasicDefines.hpp
* BasicDefinesC.hpp
* BasicDefinesJSON.hpp

-------

## BasicDefines.hpp
名前空間SRLfemの定義と、基本定数を「namespace CommonDef」で定義しています。
* PI = 3.14159265358979323846	    円周率
* PI_2 = 1.57079632679489661923		円周率/2
* PI_4 = 0.785398163397448309616	円周率/4
* L_SPEED = 299792458.0				光の速度
* MYU0 = 4.0*PI*1.0e-7				真空透磁率
* VNYU0 = 1.0/MYU0					真空磁気抵抗率
* EPS0 = 8.8541878128*1.0e-12         真空誘電率
* NORMB_EPS = 1.0e-10		        微小ゼロとみなす磁束密度ノルム
* MATH_E = 2.71828182845904523536;  eの値
* SQRT2 = 1.41421356237309504880    sqrt(2)
* SQRT1_2 = 1.41421356237309504880  1/sqrt(2)
* DegToRad = PI/180.0   [deg]から[rad]への変換係数 
* RadToDeg = 180.0/PI	[rad]から[deg]への変換係数 

## BasicDefinesC.hpp
名前空間SRLfem内に、dcomplexを定義しています
（using dcomplex = std::complex<double>　をやっているだけ）

## BasicDefinesJSON.hpp
名前空間SRLfem内に、JSON型c_jsonを定義しています。
外部ライブラリnlohmann/json.hppを読み込んで別名を付けています。
（using c_json = nlohmann::json　をやっているだけ）