/*
Copyright (c) 2023, University of Luxembourg
All rights reserved.

Redistributions and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS'
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
*/

#ifndef G2O_LINEAR_SOLVER_EIGEN_PARTIAL_LU_H
#define G2O_LINEAR_SOLVER_EIGEN_PARTIAL_LU_H

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/src/Core/util/Constants.h>
#include <cassert>

#include "g2o/core/batch_stats.h"
#include "g2o/core/linear_solver.h"
#include "g2o/core/marginal_covariance_cholesky.h"
#include "g2o/stuff/logger.h"
#include "g2o/stuff/timeutil.h"

namespace g2o {

/**
 * \brief Linear solver which uses the sparse LU solver from Eigen
 *
 * Implements a robust linear solver using Eigen's SparseLU decomposition.
 */
template <typename MatrixType>
class LinearSolverEigenLU : public LinearSolverCCS<MatrixType> {
public:
  typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SparseMatrix;
  typedef Eigen::Triplet<double> Triplet;

  using luDecomposition =
      Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>>;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  LinearSolverEigenLU() : LinearSolverCCS<MatrixType>(), _init(true) {}

  virtual bool init() override {
    _init = true;
    return true;
  }

  virtual bool solve(const SparseBlockMatrix<MatrixType> &A, double *x,
                     double *b) override {
    double t;
    if (!computeLU(A, t))
      return false;

    // Solving the linear system using LU decomposition
    VectorX::MapType xx(x, _sparseMatrix.cols());
    VectorX::ConstMapType bb(b, _sparseMatrix.cols());
    xx = _lu.solve(bb);

    if (_lu.info() != Eigen::Success) {
      G2O_ERROR("LU solve failed.");
      return false;
    }

    // Collect statistics
    G2OBatchStatistics *globalStats = G2OBatchStatistics::globalStats();
    if (globalStats) {
      globalStats->timeNumericDecomposition = get_monotonic_time() - t;
      globalStats->choleskyNNZ = 0;
    }

    return true;
  }

protected:
  bool _init;
  SparseMatrix _sparseMatrix;
  luDecomposition _lu;

  // Compute LU decomposition
  bool computeLU(const SparseBlockMatrix<MatrixType> &A, double &t) {
    // perform some operations only once upon init
    if (_init)
      _sparseMatrix.resize(A.rows(), A.cols());
    fillSparseMatrix(A, !_init);
    if (_init)
      computeSymbolicDecomposition();
    _init = false;

    t = get_monotonic_time();
    _lu.factorize(
        _sparseMatrix
            .selfadjointView<Eigen::Upper>()); // Perform Lu decomposition
    if (_lu.info() != Eigen::Success) {
      if (this->writeDebug()) {
        G2O_ERROR("Lu decomposition failed, writing debug.txt");
        A.writeOctave("debug.txt");
      } else {
        G2O_DEBUG("Lu decomposition failed");
      }
      return false;
    }

    return true;
  }

  /**
   * compute the symbolic decomposition of the matrix only once.
   * Since A has the same pattern in all the iterations, we only
   * compute the fill-in reducing ordering once and re-use for all
   * the following iterations.
   */
  void computeSymbolicDecomposition() {
    double t = get_monotonic_time();
    _lu.analyzePattern(_sparseMatrix.selfadjointView<Eigen::Upper>());
    _init = false; // Mark as initialized

    G2OBatchStatistics *globalStats = G2OBatchStatistics::globalStats();
    if (globalStats)
      globalStats->timeSymbolicDecomposition = get_monotonic_time() - t;
  }

  // Fill the sparse matrix with values from A
  void fillSparseMatrix(const SparseBlockMatrix<MatrixType> &A,
                        bool onlyValues) {
    if (onlyValues) {
      this->_ccsMatrix->fillCCS(_sparseMatrix.valuePtr(), true);
      return;
    }
    this->initMatrixStructure(A);
    _sparseMatrix.resizeNonZeros(A.nonZeros());
    int nz = this->_ccsMatrix->fillCCS(_sparseMatrix.outerIndexPtr(),
                                       _sparseMatrix.innerIndexPtr(),
                                       _sparseMatrix.valuePtr(), true);
    (void)nz;
    assert(nz <= static_cast<int>(_sparseMatrix.data().size()));
  }

  virtual bool solveBlocks_impl(
      const SparseBlockMatrix<MatrixType> &A,
      std::function<void(MarginalCovarianceCholesky &)> /*compute*/) override {
    double t;
    if (!computeLU(A, t))
      return false;

    MarginalCovarianceCholesky mcc;

    // mcc.setCholeskyFactor(
    //    _lu.matrixL().rows(),
    //   const_cast<int*>(_lu.matrixL().outerIndexPtr()),
    //  const_cast<int*>(_lu.matrixL().innerIndexPtr()),
    // const_cast<double*>(_lu.matrixL().valuePtr()),
    //    nullptr);  // No explicit permutation available
    //  compute(mcc);
    return true;
  }
};

} // namespace g2o

#include "solver_eigen_lu.hpp"

#endif
