#ifndef DEF_MT_C11_DEFTYPE0_BASICS
#define DEF_MT_C11_DEFTYPE0_BASICS

#include <cstdlib>
#include <cmath>
#include <random>

//64bitのMTを使う場合にon
#define MT_TYPE_64

/* オリジナル名前空間(静止器/回転機FEMライブラリ) */
namespace SRLfem{


/*
//=======================================================
// ■ メルセンヌ乱数使用したメソッド管理staticクラス
//=======================================================
//  メルセンヌ乱数を用いた処理を提供するクラス
//=======================================================*/
class Mt{
private:
#ifdef MT_TYPE_64
	static std::mt19937_64 mt_r;
#else
	static std::mt19937 mt_r;
#endif
public:	
	static void init_rand();														/* 乱数シード初期化 */ 
	static void init_rand(int xx);													/* 乱数シード初期化(外部指定) */ 
	static void reset(){init_rand();};
	static int genrand_int32();														/* 整数乱数デフォルト生成（0～最大まで） */
	static int mrand(int m);														/* ｍ未満の整数を返す */
	static int mrand(int a, int b);													/* a<=x<=bの範囲の一様整数乱数を返す */
	static double genrand_real1(){return(unif_rand(0.0,1.0));};						/* 0<= x <= 1の乱数発生 */
	static double unif_rand(double a, double b);									/* a<=x<=bの範囲の一様乱数を返す */
	static void make_diff_rands(const int n, const int max, int* rands);			/* n個の異なる乱数列(0～max未満)を生成する(整数) */
	static void make_diff_rands(const int n, double* rands);						/* n個の異なる乱数列(0～1)を生成する(実数) */
	/**/
	static double normal_rand(double u, double sigma, bool limit_sig=false);		/* 正規分布生成(limit_sig=true: 絶対値４sigmaまでの範囲でしか出力しない) */
	static double lognorm_rand(double u, double sigma);								/* 対数正規分布生成 */
	static double norm_CDF(double u, double sigma, double x);						/* 正規分布の累積分布関数を返す */
	static double norm_CDF(double x);												/* 標準正規分布の累積分布関数を返す */
	/**/
	static double exp_rand(double lam);												/* 指数分布 */
	static double chai_rand(double n);												/* カイ2乗分布（n：自由度） */
	static double fisher_rand(double m, double n);									/* フィッシャーのF分布（m,n：自由度） */
	static double student_t_rand(double n);											/* studentのｔ分布生成（n：自由度） */
	/**/
	static double poisson_rand(double u);											/* ポアソン分布生成 */
	static double gamma_rand(double a, double b);									/* ガンマ分布生成 */
	static double weibull_rand(double a, double b);									/* ワイブル分布生成 */
	/**/
	static bool bernoulli_rand(double a);											/* ベルヌーイ分布生成(確率aでtrueを生成) */
	static int binomial_rand(double a, int n);										/* 二項分布生成(成功確率a5の事象をn回施行し、成功した回数を返す) */
};


/* end of namespace */
}

#endif
