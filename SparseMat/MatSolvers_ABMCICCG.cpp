#include "MatSolversICCG.hpp"
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

namespace SRLfem
{
    void MatSolvers::sortAlgebraicMultiColor(int size0, const SparseMatBaseD &A, SparseMatBaseD &PA, const double *b, double *Pb, std::vector<std::vector<int>> &color_list, std::vector<int> &ordering, std::vector<int> &reverse_ordering, int num_colors)
    {
        const slv_int numRow = size0;
        int num_colors_ = num_colors;
        // Check maximum NZ element in each row
        int k = 0;
        auto col_ptr = A.getColPtr();
        auto val_ptr = A.getValuePtr();
        slv_int *start_pos = new slv_int[numRow];
        slv_int *end_pos = new slv_int[numRow];
        A.getCols(start_pos, end_pos);

        for (int i = 0; i < numRow; i++)
        {
            int mu = start_pos[i];
            while (col_ptr[mu] != i)
            {
                mu++;
            }
            if (num_colors_ < mu - start_pos[i])
            {
                num_colors_ = mu - start_pos[i] + 1;
            }
        }

        int icolor = 0;
        ordering.resize(numRow);
        reverse_ordering.resize(numRow);
        color_list.resize(num_colors_);
        std::vector<int> color(numRow);
        for (int i = 0; i < numRow; i++)
        {
            color[i] = -1;
            int mu = start_pos[i];
            while (col_ptr[mu] != i)
            {
                if (color[col_ptr[mu]] == icolor)
                {
                    icolor = (icolor + 1) % num_colors_;
                    mu = start_pos[i] - 1;
                }
                mu++;
            }

            color[i] = icolor;
            color_list[icolor].push_back(i);
            icolor = (icolor + 1) % num_colors_;
        }

        k = 0;
        for (int i = 0; i < color_list.size(); i++)
        {
            for (int j = 0; j < color_list[i].size(); j++)
            {
                ordering[color_list[i][j]] = k;
                reverse_ordering[k] = color_list[i][j];
                color_list[i][j] = k;
                k++;
            }
        }

        for (int i = 0; i < numRow; i++)
        {
            int numNZ = end_pos[i] - start_pos[i];
            for (int mu = start_pos[i]; mu < end_pos[i]; mu++)
            {
                PA.add(ordering[i], ordering[col_ptr[mu]], val_ptr[mu]);
            }
            Pb[ordering[i]] = b[i];
        }
        delete[] start_pos;
        delete[] end_pos;
    }

    void MatSolvers::makeAlgebraicBlock(int size0, const SparseMatBaseD &A,
                                        SparseMatBaseD &Mb,
                                        std::vector<std::vector<int>> &block_list,
                                        int num_blocks)
    {
        const slv_int numRow = size0;
        int max_block = (numRow + num_blocks - 1) / num_blocks;
        block_list.resize(max_block);

        auto col_ptr = A.getColPtr();
        auto val_ptr = A.getValuePtr();
        slv_int *start_pos = new slv_int[numRow];
        slv_int *end_pos = new slv_int[numRow];
        A.getCols(start_pos, end_pos);

        std::queue<int> que; // que.front(), que.pop(), que.push(i)
        std::vector<int> block_assign(numRow, -1);
        for (int i = 0; i < max_block; i++)
        {
            for (int j = 0; j < numRow; j++)
            {
                if (block_assign[j] == -1)
                {
                    block_assign[j] = i;
                    block_list[i].push_back(j);
                    for (int k = start_pos[j]; k < end_pos[j]; k++)
                    {
                        if (block_assign[col_ptr[k]] == -1)
                        {
                            que.push(col_ptr[k]);
                        }
                    }
                    while (block_list[i].size() < num_blocks && que.size() != 0)
                    {
                        if (block_assign[que.front()] == -1)
                        {
                            block_assign[que.front()] = i;
                            block_list[i].push_back(que.front());
                            for (int k = start_pos[que.front()]; k < end_pos[que.front()]; k++)
                            {
                                if (block_assign[col_ptr[k]] == -1)
                                {
                                    que.push(col_ptr[k]);
                                }
                            }
                        }
                        que.pop();
                    }
                }
                if (block_list[i].size() == num_blocks)
                {
                    std::queue<int>().swap(que);
                    break;
                }
            }
        }

        for (int i = 0; i < numRow; i++)
        {
            for (int mu = start_pos[i]; mu < end_pos[i]; mu++)
            {
                Mb.add(block_assign[i], block_assign[col_ptr[mu]], 1.0);
            }
        }
        Mb.fix();
        delete[] start_pos;
        delete[] end_pos;
    }

    void MatSolvers::sortAlgebraicBlockMultiColor(int size0,
                                                  const SparseMatBaseD &A,
                                                  SparseMatBaseD &PA,
                                                  const double *b,
                                                  double *Pb,
                                                  std::vector<std::vector<int>> &block_list,
                                                  std::vector<std::vector<int>> &color_list,
                                                  std::vector<int> &ordering,
                                                  std::vector<int> &reverse_ordering,
                                                  int num_blocks,
                                                  int num_colors)
    {
        int numRow = size0;

        auto col_ptr = A.getColPtr();
        auto val_ptr = A.getValuePtr();
        slv_int *start_pos = new slv_int[numRow];
        slv_int *end_pos = new slv_int[numRow];
        A.getCols(start_pos, end_pos);

        int max_block = (numRow + num_blocks - 1) / num_blocks;
        SparseMatBaseD Mb(max_block);
        std::vector<std::vector<int>> reverse_block_list;
        auto start = std::chrono::high_resolution_clock::now(); // record start time
        makeAlgebraicBlock(numRow, A, Mb, reverse_block_list, num_blocks);
        auto end = std::chrono::high_resolution_clock::now();             // record stop time
        std::chrono::duration<double, std::milli> duration = end - start; // calculate duration (ms)
        std::cout << "assign block: " << duration.count() << "ms" << std::endl;

        SparseMatBaseD PMb(max_block);
        double *y = new double[max_block];
        double *z = new double[max_block];
        std::vector<int> block_ordering;
        std::vector<int> reverse_block_ordering;
        start = std::chrono::high_resolution_clock::now(); // record start time
        sortAlgebraicMultiColor(max_block, Mb, PMb, y, z, color_list, block_ordering, reverse_block_ordering, num_colors);
        end = std::chrono::high_resolution_clock::now(); // record stop time
        duration = end - start;                          // calculate duration (ms)
        std::cout << "AMC: " << duration.count() << "ms" << std::endl;
        // color_list: soted by new block index
        // block_ordering: old block index -> new block index
        // reverse_block_ordering: new block index -> old block index
        delete[] y;
        delete[] z;
        int k = 0;

        ordering.resize(numRow);
        reverse_ordering.resize(numRow);

        // size of block_list is equivalent to reverse_block_list
        block_list.resize(reverse_block_list.size());

        start = std::chrono::high_resolution_clock::now(); // record start time

        // make block_list: new block index order
        for (int i = 0; i < reverse_block_ordering.size(); i++)
        {
            int bl = reverse_block_list[reverse_block_ordering[i]].size();
            block_list[i].resize(bl);

            for (int j = 0; j < bl; j++)
            {
                ordering[reverse_block_list[reverse_block_ordering[i]][j]] = k;
                reverse_ordering[k] = reverse_block_list[reverse_block_ordering[i]][j];
                block_list[i][j] = k;
                k++;
            }
        }

        end = std::chrono::high_resolution_clock::now(); // record stop time
        duration = end - start;                          // calculate duration (ms)
        std::cout << "sort by block: " << duration.count() << "ms" << std::endl;

        start = std::chrono::high_resolution_clock::now(); // record start time
        for (int i = 0; i < numRow; i++)
        {
            for (int mu = start_pos[i]; mu < end_pos[i]; mu++)
            {
                PA.add(ordering[i], ordering[col_ptr[mu]], val_ptr[mu]);
            }
            Pb[ordering[i]] = b[i];
        }

        end = std::chrono::high_resolution_clock::now(); // record stop time
        duration = end - start;                          // calculate duration (ms)
        std::cout << "make permutation matrix: " << duration.count() << "ms" << std::endl;
        delete[] start_pos;
        delete[] end_pos;
    }

    /* Conjugate gradient method for ICCG parallelized by ABMC ordering*/
    bool MatSolvers::parallelIccgSolv(int size0,
                                      const SparseMatBaseD &A,
                                      const SparseMatBaseD &L,
                                      const SparseMatBaseD &Lt,
                                      std::vector<std::vector<int>> &block_list,
                                      std::vector<std::vector<int>> &color_list,
                                      const double *diagD,
                                      const double *b,
                                      double *x,
                                      double eps,
                                      int numItr)
    {
        /*     Solution of system of linear equation                     */
        /*             A*X = B (A:symmetric)                             */
        /*      Conjugate gradient method for ICCG                       */
        /*                                                               */
        /*     ver. 1.00 2024.2.21  S. Hiruma                           */
        /* ************************************************************* */
        //    ==== input ====
        //     A......Matrix
        //     L......Matrix [L][D][Lt] = A

        //     b[i]...right hand side vector(i=0,1,2,...n-1)
        //      eps...fabsolute tolerance for residual vector
        //
        //    ==== in/output ====
        //     x[i].....input : initial guess for solution(i=0,1,2,...n-1)
        //             output : solution
        //
        //    ==== output ====
        //    r(i)...residual vector (i=0,1,2,...n-1)
        //
        //    ==== return value ====
        //        converged : Number of iteration
        //    not converged : -1

        int i, k;
        bool is_conv = false;
        double eps2;
        double r2sum;
        double rur0, rur1, pap, alpha, beta;

        int n = size0;
        Eigen::VectorXd EvecX(n);

        eps2 = eps * eps;

        try
        {
            double firstr2sum;

            // ap = new double[n];
            // p = new double[n];
            // ru = new double[n];

            Eigen::VectorXd p(n);
            Eigen::VectorXd ap(n);
            Eigen::VectorXd r(n);
            Eigen::VectorXd ru(n);
            for (int i = 0; i < n; i++)
            {
                p(i) = 0.0;
                ap(i) = 0.0;
                EvecX(i) = 0.0;
            }

            /***Cal. {AP} = [A]{x0}***/
            ap = A.matrix * EvecX;

            /***Cal. residual vector {r0} = {b} - [A]{X0} ***/
            r2sum = 0.0;
#pragma omp parallel for
            for (i = 0; i < n; i++)
            {
                r(i) = b[i] - ap(i);
                // r2sum += r(i) * r(i);
            }
            r2sum = r.dot(r);
            firstr2sum = fabs(r2sum);

            // std::cout << "Initial residual = " << sqrt(firstr2sum) << std::endl;

            /***solve [L][D][Lt]{ru} = {r}***/
            parallelIcSolv(size0, L, Lt, block_list, color_list, diagD, r, ru);

            rur0 = 0.0;
#pragma omp parallel for
            for (i = 0; i < n; i++)
            {
                p(i) = ru(i);
                // rur0 += r(i) * ru(i);
            }
            rur0 = r.dot(ru);

            /* 最良結果の保存用（フラグがonなら） */
            double *best_results = nullptr;
            double best_resi_value = 1.0e+6;
            if (is_save_best)
            {
                best_results = new double[size0];
            }
            int bad_counter = 0;

            /**** iteration loop ****/
            for (k = 0; k < numItr; k++)
            {

                /* *** {ap} = [A]{p} *** */
                ap = A.matrix * p;

                pap = p.dot(ap);

                alpha = rur0 / pap;
#pragma omp parallel for
                for (i = 0; i < n; i++)
                {
                    x[i] += alpha * p(i);
                    r(i) -= alpha * ap(i);
                }
                r2sum = r.dot(r);

                /* フラグがonなら、残差保存 */
                if (is_save_residual_log)
                {
                    residual_log.push_back(sqrt(r2sum));
                }

                /*** check convergence**/
                if (fabs(r2sum) <= eps2)
                {
                    // std::cout << "---- converged ------" << std::endl;
                    // std::cout << "Number of iteration = " << k + 1 << std::endl;
                    numItr = k + 1;
                    is_conv = true;
                    break;
                }

                if (normR < best_resi_value)
                {
                    best_resi_value = normR;
                    /* 最良値の更新(フラグがonなら) */
                    if (is_save_best)
                    {
                        for (slv_int i = 0; i < size; i++)
                        {
                            best_results[i] = results[i];
                        }
                    }
                }

                /* check divergence */
                if (diverge_judge_type == 1)
                {
                    /* 最良値×val以下なら、発散カウント初期化 */
                    if (normR < best_resi_value * bad_div_val)
                    {
                        bad_counter = 0;
                    }
                    /* 最良値×val以上なら、発散カウント＋ */
                    if (normR >= best_resi_value * bad_div_val)
                    {
                        bad_counter++;
                    }
                    /* 発散カウントが閾値オーバー＝発散扱いで終わる */
                    if (bad_counter >= bad_div_count_thres)
                    {
                        is_conv = false;
                        break;
                    }
                }

                /***solve [L][D][Lt]{ru} = {r}***/
                parallelIcSolv(size0, L, Lt, block_list, color_list, diagD, r, ru);

                rur1 = 0.0;
                // for (i = 0; i < n; i++)
                // {
                //     rur1 += r(i) * ru(i);
                // }
                rur1 = r.dot(ru);
                beta = rur1 / rur0;
                rur0 = rur1;
#pragma omp parallel for
                for (i = 0; i < n; i++)
                {
                    p(i) = ru(i) + beta * p(i);
                }
            }

            /*** not convergence**/
            if (!is_conv)
            {
                numItr = -1;
                std::cout << "--- not converged ---" << std::endl;
                std::cout << "Number of iteration = " << k << std::endl;

                /* 最良値を代入(フラグがonなら) */
                if (is_save_best)
                {
                    for (slv_int i = 0; i < size; i++)
                    {
                        results[i] = best_results[i];
                    }
                    delete[] best_results;
                }
            }
        }
        catch (...)
        {
            throw;
        }
        /***  termination ***/
        return is_conv;
    }

    /* *********************************************************** */
    void MatSolvers::parallelIcSolv(int size0, const SparseMatBaseD &L, const SparseMatBaseD &Lt, std::vector<std::vector<int>> &block_list, std::vector<std::vector<int>> &color_list,
                                    const double *diagD, const Eigen::VectorXd &b, Eigen::VectorXd &x)
    {
        //
        //              Solve system of linear equations
        //                  (L*D*Lt)*X = B
        //             in ICCG method (Lt = transpose of L)
        //             LDLt is incomplete Cholesky decomposition
        //             arranged by ICDCMP
        //
        /*     ver. 1.00 2024.2.21  S. Hiruma
        /* *********************************************************** */
        int n, k, l, mu, nu;
        int col;
        int bl;
        int idx;
        double t;
        double *y;
        int *columnP;
        double *valueP;
        int numNZ;

        n = size0;

        try
        {
            auto col_ptr = L.getColPtr();
            auto val_ptr = L.getValuePtr();
            slv_int *start_pos = new slv_int[n];
            slv_int *end_pos = new slv_int[n];
            L.getCols(start_pos, end_pos);

            auto col_ptr2 = Lt.getColPtr();
            auto val_ptr2 = Lt.getValuePtr();
            slv_int *start_pos2 = new slv_int[n];
            slv_int *end_pos2 = new slv_int[n];
            Lt.getCols(start_pos2, end_pos2);

            y = new double[n];
            /*     ---- forward substitution ----*/
            /* solve [L]{y} = {b} */
            for (col = 0; col < color_list.size(); col++)
            {
#pragma omp parallel for private(bl, l, mu, t, idx)
                for (k = 0; k < color_list[col].size(); k++)
                {
                    bl = color_list[col][k];
                    for (l = 0; l < block_list[bl].size(); l++)
                    {
                        idx = block_list[bl][l];
                        t = b(idx);
                        for (mu = start_pos[idx]; mu < end_pos[idx] - 1; mu++)
                        { // numNZ-1 is the diagnal !!
                            t -= val_ptr[mu] * y[col_ptr[mu]];
                        }
                        y[idx] = t / val_ptr[end_pos[idx] - 1]; // SparseMat IC version
                    }
                }
            }
            /* init x[] */
            for (k = 0; k < n; k++)
            {
                x(k) = y[k];
            }

            /*     ---- backward substitution ----*/
            /* solve [D][L]{x} = {y} */
            for (col = color_list.size() - 1; col >= 0; col--)
            {
#pragma omp parallel for private(bl, l, nu, t, idx)
                for (k = color_list[col].size() - 1; k >= 0; k--)
                {
                    bl = color_list[col][k];
                    for (l = block_list[bl].size() - 1; l >= 0; l--)
                    {
                        idx = block_list[bl][l];
                        t = 0.0;
                        for (nu = start_pos2[idx] + 1; nu < end_pos2[idx]; nu++)
                        {
                            // t -= 1.0 / val_ptr2[start_pos2[idx]] * val_ptr2[nu] * x(col_ptr2[nu]);
                            t -= val_ptr2[nu] * x(col_ptr2[nu]);
                        }
                        // t /= val_ptr2[start_pos2[idx]];
                        t *= diagD[idx];
                        x(idx) += t;
                    }
                }
            }
            delete[] y;
            delete[] start_pos;
            delete[] end_pos;
            delete[] start_pos2;
            delete[] end_pos2;
        }
        catch (...)
        {
            throw;
        }
    }

    bool MatSolvers::solveICCGwithABMC(const slv_int size0, const double conv_cri, const int max_ite, const double accera, const double normB, const SparseMat &matA, double *vec_b, double *vec_x, int num_blocks, int num_colors)
    {
        bool is_conv = false;

        double epsilon = conv_cri;
        int n = size0;

        // 残差ベクトル
        double *vec_r = new double[n];

        double sum = 0.0;
        for (int i = 0; i < n; i++)
        {
            vec_x[i] = 0.0;
            vec_r[i] = 0.0;
        }
        epsilon *= normB;

        SparseMat PA(n);
        double *vec_Pb = new double[n];

        std::vector<std::vector<int>> block_list;
        std::vector<std::vector<int>> color_list;
        std::vector<int> ordering;
        std::vector<int> reverse_ordering;

        std::cout << "start ABMC" << std::endl;
        auto start = std::chrono::high_resolution_clock::now(); // record start time
        MatSolvers::sortAlgebraicBlockMultiColor(n, *(matA.matrix), *(PA.matrix), vec_b, vec_Pb, block_list, color_list, ordering, reverse_ordering, num_blocks, num_colors);
        auto end = std::chrono::high_resolution_clock::now();             // record stop time
        std::chrono::duration<double, std::milli> duration = end - start; // calculate duration (ms)
        std::cout << "ABMC(total): " << duration.count() << "ms" << std::endl;
        PA.fix();

        double ga = accera;
        double *diagD = new double[n];
        start = std::chrono::high_resolution_clock::now(); // record start time
        SparseMat matL = PA.IC_decomp(diagD, accera);
        end = std::chrono::high_resolution_clock::now(); // record stop time
        duration = end - start;                          // calculate duration (ms)
        std::cout << "IC_decomp: " << duration.count() << "ms" << std::endl;

        start = std::chrono::high_resolution_clock::now(); // record start time
        SparseMat matL_tr = matL.trans();
        end = std::chrono::high_resolution_clock::now(); // record stop time
        duration = end - start;                          // calculate duration (ms)
        std::cout << "Transpose matrix: " << duration.count() << "ms" << std::endl;

        start = std::chrono::high_resolution_clock::now(); // record start time
        is_conv = MatSolvers::parallelIccgSolv(n, *(PA.matrix), *(matL.matrix), *(matL_tr.matrix), block_list, color_list, diagD, vec_Pb, vec_x, epsilon, max_ite);
        end = std::chrono::high_resolution_clock::now(); // record stop time
        duration = end - start;                          // calculate duration (ms)
        std::cout << "Iteration: " << duration.count() << "ms" << std::endl;

        for (int i = 0; i < n; i++)
        {
            vec_r[i] = vec_x[i];
        }

        double normX = 0.0;
        for (int i = 0; i < n; i++)
        {
            vec_x[reverse_ordering[i]] = vec_r[i];
            normX += vec_x[reverse_ordering[i]] * vec_x[reverse_ordering[i]];
        }
        std::cout << "normX = " << sqrt(normX) << std::endl;

        delete[] vec_r;
        delete[] vec_Pb;
        delete[] diagD;

        return is_conv;
    }
} // namespace SRLfem
