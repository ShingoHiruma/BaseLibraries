#ifndef DEF_MY_HEADER_1_COMMONFUNCS_BASICS
#define DEF_MY_HEADER_1_COMMONFUNCS_BASICS

#include <cmath>
#include <iostream>
#include <BasicDefines.hpp>

/* 簡略化の記号 */


/* オリジナル名前空間(静止器/回転機FEMライブラリ) */
namespace SRLfem{


/*
//=======================================================
// ■ CommonFuncs:共有用関数たち管理クラス
//=======================================================
// 共有的に使う関数を提供
//=======================================================*/
class CommonFuncs{
public:
	static int round(double m);
	inline static double sigmoid(double x, double constK){return( 1.0/(1.0 + exp(-1.0*constK*x))  );};	/* シグモイド */
	/**/
	static void setGauss(double* tg, double* wg, int numP );																	/* ガウス積分点[-1, 1] */
	static void setGaussHEX(double* tg_x, double* tg_y, double* tg_z, double* wg, int numP );									/* ガウス積分点(六面体用) */
	static void setGaussPRI(double* tg_x, double* tg_y, double* tg_z, double* wg_xy, double* wg_z, int numP_xy, int numP_z);	/* ガウス積分点(プリズム用) */
	/**/
	static void calcGaussPoint(double* gauss_p, double* gauss_w, const int num_gauss);											/* ガウス積分点計算（積分点任意ver） */
	/**/
	/**/
	inline static void make_three_phase(double* three_rad, double u_phase_rad, bool minus=true);								/* 三相の位相差を計算して返す（u_radから計算）。デフォルト相順はU=0，V=-120，W=-240 */
	/**/
	/**/
	/*---------------------------*/
	/*--離散フーリエ変換系-------*/
	static void dft(int size, double *wave_data, double *dft_r, double *dft_im);									/* 通常離散フーリエ */
	static void dft(int size, double *wave_data, double *dft_r, double *dft_im, double *ampl, double *phase);		/* 通常離散フーリエ・各次数の振幅と位相も計算 */
	static void idft(int size, double *dft_r, double *dft_im, double *wave_data);									/* 通常離散逆フーリエ */
	static void dft_filtering(int size, int target, double *dft_r, double *dft_im);									/* target次数の成分のみ抽出するフィルタ */
	/*---------------------------*/
	static double calc_phase_deg(int size, double* wave1, double* wave2);											/* ２つの波形の基本波の位相差をdeg単位で返す */
	static void ThreeFullWaveMaker(int sub_wave_count, int wave_points, int& full_points, double** sub_waves, double*** full_wave);	/* 3相対称の部分波形からフル波形を作る。sub_wave_countは部分波形の電気角の角度(waveはU[0deg],V[-120deg],W[-240deg]の順、countは60,120,180のどれか) */
};


/* end of namespace */
}

#endif
