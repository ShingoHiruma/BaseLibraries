#include "CommonFuncs.hpp"

#include <BasicDefinesc.hpp>


/* オリジナル名前空間(静止器/回転機FEMライブラリ) */
namespace SRLfem {


/*//=======================================================
// ●　通常離散フーリエ
//=======================================================*/
void CommonFuncs::dft(int size, double* wave_data, double* dft_r, double* dft_im) {
	/* まず元波形を複素数化 */
	dcomplex* wave = new dcomplex[size];
	dcomplex* wave_dft = new dcomplex[size];
	for (int i = 0; i < size; i++) {
		wave[i].real(wave_data[i]);
		wave[i].imag(0);
	}

	const dcomplex j2pn = dcomplex(0.0, -2.0 * CommonDef::PI / size);
	for (int k = 0; k < size; k++) {
		wave_dft[k].real(0.0);
		wave_dft[k].imag(0.0);
		for (int i = 0; i < size; i++) {
			const double coef = k * i;
			wave_dft[k] += wave[i] * exp(j2pn * coef);
		}
		dft_r[k] = wave_dft[k].real();
		dft_im[k] = wave_dft[k].imag();
	}
	delete[] wave;
	delete[] wave_dft;
}

/*//=======================================================
// ●　通常離散フーリエ・各次数の振幅と位相も計算
//=======================================================*/
void CommonFuncs::dft(int size, double* wave_data, double* dft_r, double* dft_im, double* ampl, double* phase) {
	/* まず変換 */
	CommonFuncs::dft(size, wave_data, dft_r, dft_im);

	/* 出てきたスペクトルから振幅と位相を計算する */
	for (int k = 0; k < size / 2; k++) {
		dcomplex wave_dft = dcomplex(dft_r[k], dft_im[k]);
		wave_dft /= (double)size;
		ampl[k] = abs(wave_dft);
		/* エイリアシングで対象側に半分あるので、直流分以外は2倍する */
		if (k > 0) {
			ampl[k] *= 2;
		}
		phase[k] = atan2(dft_im[k], dft_r[k]);
		phase[k] *= 180.0 / CommonDef::PI;
	}
}

/*//=======================================================
// ●　通常離散逆フーリエ
//=======================================================*/
void CommonFuncs::idft(int size, double* dft_r, double* dft_im, double* wave_data) {
	/* まず元波形を複素数化 */
	dcomplex* wave = new dcomplex[size];
	dcomplex* wave_dft = new dcomplex[size];
	for (int i = 0; i < size; i++) {
		wave_dft[i].real(dft_r[i]);
		wave_dft[i].imag(dft_im[i]);
	}

	const dcomplex j2pn = dcomplex(0.0, 2.0 * CommonDef::PI / size);
	for (int k = 0; k < size; k++) {
		wave[k].real(0.0);
		wave[k].imag(0.0);
		for (int i = 0; i < size; i++) {
			const double coef = k * i;
			wave[k] += wave_dft[i] * exp(j2pn * coef);
		}
		wave[k] /= (double)(size);
		wave_data[k] = wave[k].real();
	}
	delete[] wave;
	delete[] wave_dft;
}

/*//=======================================================
// ●　target次数の成分のみ抽出するフィルタ
//=======================================================*/
void CommonFuncs::dft_filtering(int size, int target, double* dft_r, double* dft_im) {
	for (int k = 0; k < size; k++) {
		/* エイリアシングで、対象側にもあるのでそこも保存 */
		bool is_target = (k == target || k == size - target);
		if (!is_target) {
			dft_r[k] = 0;
			dft_im[k] = 0;
		}
	}
}

/*//=======================================================
// ●　２つの波形の基本波の位相差をdeg単位で返す
//=======================================================*/
double CommonFuncs::calc_phase_deg(int size, double* wave1, double* wave2){
	double* dft_r = new double[size];
	double* dft_im = new double[size];
	double* ampl = new double[size];
	double* phase = new double[size];
	/* まず変換 */
	CommonFuncs::dft(size, wave1, dft_r, dft_im, ampl, phase);
	double base_rad1 = phase[1];
	/* 位相を引く */ 
	CommonFuncs::dft(size, wave2, dft_r, dft_im, ampl, phase);
	base_rad1 -= phase[1];

	delete[] dft_r;
	delete[] dft_im;
	delete[] ampl;
	delete[] phase;
	/* degに直して返す */
	return( base_rad1*180.0/CommonDef::PI  );
}

/*//=======================================================
// ●　3相対称の部分波形からフル波形を作る(waveはU[0deg],V[-120deg],W[-240deg]の順)
//     sub_wave_countは部分波形の電気角の角度
//     countは60,120,180のどれか
//=======================================================*/
void CommonFuncs::ThreeFullWaveMaker(int sub_wave_count, int wave_points, int& full_points, double** sub_waves, double*** full_wave) {
	int full_wave_count;
	double** temp_full_wave = new double* [3];
	/* 電気角1/6からフル波形にする */
	if (sub_wave_count == 60) {
		full_wave_count = 6 * wave_points;
		for (int i = 0; i < 3; i++) {
			temp_full_wave[i] = new double[full_wave_count];
			for (int j = 0; j < wave_points; j++) {
				temp_full_wave[i][j] = sub_waves[i][j];
			}
		}
		for (int j = wave_points; j < full_wave_count; j++) {
			temp_full_wave[0][j] = -1.0 * temp_full_wave[1][j - wave_points];
			temp_full_wave[1][j] = -1.0 * temp_full_wave[2][j - wave_points];
			temp_full_wave[2][j] = -1.0 * temp_full_wave[0][j - wave_points];
		}
		/* 電気角1/3からフル波形にする */
	}else if (sub_wave_count == 120) {
		full_wave_count = 3 * wave_points;
		for (int i = 0; i < 3; i++) {
			temp_full_wave[i] = new double[full_wave_count];
			for (int j = 0; j < wave_points; j++) {
				temp_full_wave[i][j] = sub_waves[i][j];
			}
		}
		for (int j = wave_points; j < full_wave_count; j++) {
			temp_full_wave[0][j] = temp_full_wave[2][j - wave_points];
			temp_full_wave[1][j] = temp_full_wave[0][j - wave_points];
			temp_full_wave[2][j] = temp_full_wave[1][j - wave_points];
		}
		/* 電気角1/2からフル波形にする */
	}else if (sub_wave_count == 180) {
		full_wave_count = 2 * wave_points;
		for (int i = 0; i < 3; i++) {
			temp_full_wave[i] = new double[full_wave_count];
			for (int j = 0; j < wave_points; j++) {
				temp_full_wave[i][j] = sub_waves[i][j];
			}
		}
		for (int j = wave_points; j < full_wave_count; j++) {
			temp_full_wave[0][j] = -1.0 * temp_full_wave[0][j - wave_points];
			temp_full_wave[1][j] = -1.0 * temp_full_wave[1][j - wave_points];
			temp_full_wave[2][j] = -1.0 * temp_full_wave[2][j - wave_points];
		}
	}else {
		std::cout << "Three sub wave is not correct ! " << std::endl;
		*full_wave = nullptr;
		delete[] temp_full_wave;
		return;
	}
	/* 実体を渡して終わる */
	full_points = full_wave_count;
	*full_wave = std::move(temp_full_wave);
	temp_full_wave = nullptr;
}


/*end of namespace */
};
