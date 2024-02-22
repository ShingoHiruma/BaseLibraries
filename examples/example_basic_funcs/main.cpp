#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
using namespace std;


#include <BasicFuncs/Mt.hpp>
#include <BasicFuncs/CommonFuncs.hpp>



/*//=======================================================
// ●　メイン関数
//=======================================================*/
int main(int argc, char *argv[]){

	int seed = (unsigned)time(NULL);
	cout << "Seed: "<< seed << endl;
	/* 乱数シードのセット */
	srand(seed);
	SRLfem::Mt::init_rand(seed);

	/* ０～１乱数の生成 */
	double rand1 = SRLfem::Mt::genrand_real1();
	cout << rand1 << endl;

	/* a～bの一様乱数の生成 */
	double a = 10;
	double b = 50;
	double rand2 = SRLfem::Mt::unif_rand(a, b);
	cout << rand2 << endl;


	/* 平均u、分散sig*sigの正規乱数の生成 */
	double u = 5.0;
	double sig = 100;
	double rand3 = SRLfem::Mt::normal_rand(u, sig);
	cout << rand3 << endl;


	/*====================================*/
	/*====================================*/
	/*====================================*/
	cout << "------------------------" << endl;
	cout << "------------------------" << endl;
	cout << "------------------------" << endl;


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

	return 1;

}
