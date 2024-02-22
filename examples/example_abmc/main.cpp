#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
using namespace std;

#include <SparseMat/SparseMat.hpp>
#include <SparseMat/SparseMatC.hpp>
#include <SparseMat/SparseMatOperators.hpp>
#include <SparseMat/MatSolvers.hpp>
#include <SparseMat/MatSolvers_ABMCICCG.hpp>

#ifdef IS_WINDOWS_SISTEM
#ifdef _DEBUG
#pragma comment(lib, "libSparseMat_Deb.lib")
#else
#pragma comment(lib, "libSparseMat.lib")
#endif
#endif

/*//=======================================================
// ●　メイン関数
//=======================================================*/
int main(int argc, char *argv[])
{
	std::cout << "Read Matrix Market format" << std::endl;
	std::string filename = "./data/mat_large.dat";
	std::string filename2 = "./data/vec_b_large.dat";
	std::ifstream file(filename);
	std::ifstream file2(filename2);

	if (!file.is_open())
	{
		std::cerr << "Failed to open file: " << filename << std::endl;
		return false;
	}

	std::string line;
	// Skip comments and empty lines
	while (std::getline(file, line))
	{
		if (line[0] != '%')
			break;
	}

	// Read matrix dimensions and number of non-zero entries
	std::istringstream iss(line);
	int numRows, numCols, numNZ;
	iss >> numRows >> numCols >> numNZ;
	SRLfem::SparseMat matrix = SRLfem::SparseMat(numRows);

	double *vec_b = new double[numRows];
	for (int i = 0; i < numRows; i++)
	{
		file2 >> vec_b[i];
	}
	int currentRow = 0;

	for (int i = 0; i < numNZ; i++)
	{
		int row, col;
		double value;
		file >> row >> col >> value;
		row--;
		col--; // Adjust from 1-based to 0-based indexing

		matrix.add(row, col, value);

		if (row > currentRow)
		{
			currentRow = row;
			// std::cout << row << "," << col << "," << value << std::endl;
		}
	}

	matrix.fix();

	/* initialization */
	double *results00 = new double[numRows];
	for (int i = 0; i < numRows; i++)
	{
		results00[i] = 0;
	}
	double *results01 = new double[numRows];
	for (int i = 0; i < numRows; i++)
	{
		results01[i] = 0;
	}

	cout << "start" << endl;

	/*ICCG with ABMC ordering*/
	double epsilon = 1.0e-10;
	int max_itr = numRows;
	double accera = 1.1;
	double normB = 0.0;
	for (int i = 0; i < numRows; i++)
	{
		normB += vec_b[i] * vec_b[i];
	}
	normB = sqrt(normB);

	std::cout << numRows << std::endl;

	auto start = std::chrono::high_resolution_clock::now(); // record start time
	omp_set_num_threads(20);
	// bool bl1 = SRLfem::MatSolvers::solveICCG(numRows, epsilon, max_itr, accera, matrix, vec_b, results00);
	bool bl2 = SRLfem::ABMCICCG::solveICCGwithABMC(numRows, epsilon, max_itr, accera, normB, matrix, vec_b, results00, 512, 40);
	auto end = std::chrono::high_resolution_clock::now();			  // record end time
	std::chrono::duration<double, std::milli> duration = end - start; // calculate duration (ms)
	std::cout << "Total: " << duration.count() << "ms" << std::endl;

	return 1;
}
