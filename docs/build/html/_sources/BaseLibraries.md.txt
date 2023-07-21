# BasicFuncs
有限要素法や最適化など、数値計算や進化計算関連で使う基本操作をライブラリ化したものです。

-------

## Mtクラス
C++11のメルセンヌツイスタのラッパとして働く、乱数関連のライブラリです。「BasicFuncs/Mt.hpp」をインクルードして使ってください。
全てStatic関数です。インスタンスを作成せずに呼び出してください。

###### void init_rand()
乱数シード初期化の初期化メソッドです。ランダムに初期化します

###### void init_rand(int xx)
乱数シード初期化メソッドです。引数で指定した整数を乱数シードにして初期化します。

###### int genrand_int32()
整数乱数をメルセンヌ・ツイスタで生成します。範囲は「0～int最大まで」です。

###### int mrand(int m)
メルセンヌ・ツイスタで、整数ｍ未満の整数値をランダムに返します。

###### int mrand(int a, int b)
メルセンヌ・ツイスタで、a<=x<=bの範囲の一様整数乱数を返します。

###### double genrand_real1()
メルセンヌ・ツイスタで、0<= x <= 1の乱数を発生します。

###### double unif_rand(double a, double b)
メルセンヌ・ツイスタで、a<=x<=bの範囲の一様乱数を発生します。

###### void make_diff_rands(const int n, const int max, int* rands)
n個の異なる乱数列(0～max未満)を整数で生成します。

###### void make_diff_rands(const int n, double* rands)
n個の異なる乱数列(0～1)を実数で生成します。

###### double normal_rand(double u, double sigma, bool limit_sig=false)
正規分布に従う乱数を生成します(limit_sig=true: 絶対値４sigmaまでの範囲でしか出力しない)

###### double lognorm_rand(double u, double sigma)
対数正規分布に従う乱数を生成します。

###### double norm_CDF(double u, double sigma, double x)
正規分布の累積分布関数を返します。

###### double norm_CDF(double x);
標準正規分布の累積分布関数を返します。

###### double exp_rand(double lam)
指数分布に従う乱数を生成します。

###### double chai_rand(double n)
カイ2乗分布に従う乱数を生成します（n：自由度）

###### double fisher_rand(double m, double n)
フィッシャーのF分布（m,n：自由度） に従う乱数を生成します。
###### double student_t_rand(double n)
studentのｔ分布生成（n：自由度） に従う乱数を生成します。

###### double poisson_rand(double u)
ポアソン分布に従う乱数を生成します。

###### double gamma_rand(double a, double b)
ガンマ分布に従う乱数を生成します。

###### double weibull_rand(double a, double b)
ワイブル分布に従う乱数を生成します。

###### bool bernoulli_rand(double a)
ベルヌーイ分布(確率aでtrueを生成)に従う乱数を生成します。

###### int binomial_rand(double a, int n)
二項分布(成功確率a5の事象をn回施行し、成功した回数を返す)に従う乱数を生成します。

-------

## CommonFuncsクラス
数値解析に使いそうな操作をまとめたライブラリです。ガウス積分点の生成などがあります。「BasicFuncs/CommonFuncs.hpp」をインクルードして使ってください。
全てStatic関数です。インスタンスを作成せずに呼び出してください。

###### int round(double m)
四捨五入メソッドです。いまはstdで普通にありますが、一応。

###### double sigmoid(double x, double constK)
入力ｘをシグモイド関数に入力した結果を返します。Kは係数です。
f(x) = 1 / (1+exp(-K*x))

###### void setGauss(double* tg, double* wg, int numP )
ガウス積分点を作成します。積分点数はnumPで、積分点をtgに、積分重みをwgに格納します。
numPは2~12で指定してください。
ポインタは事前にnumPのサイズで確保しておいてください。

###### void setGaussHEX(double* tg_x, double* tg_y, double* tg_z, double* wg, int numP )
六面体用のガウス積分点を作成します。積分点数はnumPです。
積分点のx,y,z成分ををtg_x、_y、_zに格納し、その積分重みをwgに格納します。
ポインタは事前にnumPのサイズで確保しておいてください。

###### void setGaussPRI(double* tg_x, double* tg_y, double* tg_z, double* wg_xy, double* wg_z, int numP_xy, int numP_z)
プリズム用のガウス積分点を作成します。xy面(三角形面)の積分点はnumP_xy、z方向はnumP_zです。
積分点のx,y,z成分ををtg_x、_y、_zに格納し、その積分重みをwgに格納します。
ポインタは事前にnumPのサイズで確保しておいてください。

###### void calcGaussPoint(double* gauss_p, double* gauss_w, const int num_gauss)
ガウス積分点を作成します。積分点数はnumPで、積分点をtgに、積分重みをwgに格納します。
ガウスの積分点数に制限はありません。
ポインタは事前にnumPのサイズで確保しておいてください。

###### void make_three_phase(double* three_rad, double u_phase_rad, bool minus=true)
三相の位相差を計算して返す（u_radから計算）。デフォルト相順はU=0，V=-120，W=-240。
つまり、
three_rad[0]=u_phase_rad
three_rad[1]=u_phase_rad-120deg*PI/180
three_rad[2]=u_phase_rad-240deg*PI/180

###### void dft(int size, double *wave_data, double *dft_r, double *dft_im)
波形wave_dataの離散フーリエ変換を行います。wave_dataのデータ数は「size」です。
dft_rに変換後の実部、dft_imに虚部を代入します。
なお、変換結果のsize/2以降はエイリアシングで意味がない結果です。
各データのポインタは事前に確保しておいてください。

###### void dft(int size, double *wave_data, double *dft_r, double *dft_im, double *ampl, double *phase)
波形wave_dataの離散フーリエ変換を行います。wave_dataのデータ数は「size」です。
dft_rに変換後の実部、dft_imに虚部を代入します。
また、各次数の大きさと位相[rad]を、amplとphaseに代入します。
各データのポインタは事前に確保しておいてください。

###### void idft(int size, double *dft_r, double *dft_im, double *wave_data)
離散フーリエ後の結果dft_rとdft_imから、逆変換を行った波形を作成します。
サイズはsizeで、結果波形はwave_dataに格納されます。
各データのポインタは事前に確保しておいてください。

###### void dft_filtering(int size, int target, double *dft_r, double *dft_im)
離散フーリエ後の結果dft_rとdft_imの、target次数以外をゼロにします。

###### double calc_phase_deg(int size, double* wave1, double* wave2)
２つの波形の基本波の位相差をdeg単位で返します。
両波形のサイズは共にsizeです。

###### void ThreeFullWaveMaker(int sub_wave_count, int wave_points, int& full_points, double** sub_waves, double*** full_wave)
3相対称の部分波形からフル波形を作ります。
sub_wave_countは部分波形の電気角の角度で、60，120，180のいずれかが入ります。
wave_pointsは部分波形の点数で、値はsub_wavesで渡します。(waveはU[0deg]V[-120deg],W[-240deg]の順)
部分波形をフル波形にした結果は、full_waveに格納されます。ポインタは事前に確保せずにOKです。
（sub_wave_countが60でwave_pointsが50なら、フルの波形サイズは6*50＝300になり、full_fave[3][300]なるサイズの波形が出来ます）

