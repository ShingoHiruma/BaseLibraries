#include <vector>
#include <queue>
#include "MatSolvers.hpp"
#include <chrono>
#include <omp.h>

/*
    This file is part of JP_Mars, a C++ template library for linear algebra.

    Copyright (C) 2024 Shingo Hiruma <hiruma.shingo.3w@kyoto-u.ac.jp>

    This Source Code Form is subject to the terms of the Mozilla
    Public License v. 2.0. If a copy of the MPL was not distributed
    with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

    This is a parallel version of the ICCG solver using the Algebraic Block Multi-Color (ABMC) ordering.
    The solver is based on the following paper:
    "Algebraic Block Multi-Color Ordering Method for Parallel Multi-Threaded Sparse Triangular Solver in ICCG Method"
    proposed by T. Iwashita, H. Nakashima, and Y. Takahashi, in the Conference: Parallel & Distributed Processing Symposium (IPDPS), 2012 IEEE 26th International.
*/

/* 専用名前空間 */
namespace SRLfem
{
    class MatSolversABMCICCG : public MatSolvers
    {
    private:
        // Sort by Algebraic Multi-Color Ordering
        void sortAlgebraicMultiColor(int size0,
                                     const SparseMatBaseD &A,
                                     SparseMatBaseD &PA,
                                     const double *b,
                                     double *Pb,
                                     std::vector<std::vector<int>> &color_list,
                                     std::vector<int> &ordering,
                                     std::vector<int> &reverse_ordering,
                                     int num_colors);

        // Make Algebraic Block
        void makeAlgebraicBlock(int size0,
                                const SparseMatBaseD &A,
                                SparseMatBaseD &Mb,
                                std::vector<std::vector<int>> &block_list,
                                int num_blocks);

        // Sort by Algebraic Block Multi-Color Ordering
        void sortAlgebraicBlockMultiColor(int size0,
                                          const SparseMatBaseD &A,
                                          SparseMatBaseD &PA,
                                          const double *b,
                                          double *Pb,
                                          std::vector<std::vector<int>> &block_list,
                                          std::vector<std::vector<int>> &color_list,
                                          std::vector<int> &ordering,
                                          std::vector<int> &reverse_ordering,
                                          int num_blocks,
                                          int num_colors);

        // ICCG Method using ABMC ordering
        bool parallelIccgSolv(int size0,
                              const SparseMatBaseD &A,
                              const SparseMatBaseD &L,
                              const SparseMatBaseD &Lt,
                              std::vector<std::vector<int>> &block_list,
                              std::vector<std::vector<int>> &color_list,
                              const double *diagD,
                              const double *b,
                              double *x,
                              double eps,
                              int numItr);

        // solve [L][D][Lt]{x} = {b} for ICCG solver using ABMC ordering
        void parallelIcSolv(int size0,
                            const SparseMatBaseD &L,
                            const SparseMatBaseD &Lt,
                            std::vector<std::vector<int>> &block_list,
                            std::vector<std::vector<int>> &color_list,
                            const double *diagD,
                            const Eigen::VectorXd &b, Eigen::VectorXd &x);

    public:
        MatSolversABMCICCG();
        ~MatSolversABMCICCG();

        void setDiagScale(bool bl) { is_diag_scale = bl; };
        void setSaveBest(bool bl) { is_save_best = bl; };
        void setSaveLog(bool bl) { is_save_residual_log = bl; };
        void getResidualLog(std::vector<double> &log);
        void setDirvegeType(int x) { diverge_judge_type = x; };
        void setBadDivVal(double x) { bad_div_val = x; };
        void setBadDivCount(int x) { bad_div_count_thres = x; };

        bool solveICCGwithABMC(const slv_int size0,
                               const double conv_cri,
                               const int max_ite,
                               const double accera,
                               const double normB,
                               const SparseMat &matA,
                               double *vec_b,
                               double *vec_x,
                               int num_blocks,
                               int num_colors);
    };
} // namespace SRLfem
