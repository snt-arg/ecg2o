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

#ifndef G2O_LINEAR_SOLVER_EIGEN_QR_H
#define G2O_LINEAR_SOLVER_EIGEN_QR_H

#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <Eigen/src/Core/util/Constants.h>
#include <cassert>

#include "g2o/core/batch_stats.h"
#include "g2o/core/linear_solver.h"
#include "g2o/core/marginal_covariance_cholesky.h"
#include "g2o/stuff/logger.h"
#include "g2o/stuff/timeutil.h"

namespace g2o {

/**
 * \brief linear solver which uses the sparse QR solver from Eigen
 *
 * Has no dependencies except Eigen. Hence, should compile almost everywhere
 * without too much issues. Performance should be similar to CSparse.
 */
template <typename MatrixType>
class LinearSolverEigenQR : public LinearSolverCCS<MatrixType> {
public:
  typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SparseMatrix;
  typedef Eigen::Triplet<double> Triplet;

  using sparseQRDecomposition =
      Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>>;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
  LinearSolverEigenQR() : LinearSolverCCS<MatrixType>(), _init(true) {}

  virtual bool init() {
    _init = true;
    return true;
  }

  bool solve(const SparseBlockMatrix<MatrixType> &A, double *x, double *b) {
    double t;
    bool qrState = computeSparseQR(A, t);
    if (!qrState)
      return false;

    // Solving the system
    Eigen::VectorXd::MapType xx(x, _sparseMatrix.cols());
    Eigen::VectorXd::ConstMapType bb(b, _sparseMatrix.cols());

    xx = _sparseQR.solve(bb);
    G2OBatchStatistics *globalStats = G2OBatchStatistics::globalStats();
    if (globalStats) {
      globalStats->timeNumericDecomposition = get_monotonic_time() - t;
      globalStats->choleskyNNZ = _sparseQR.matrixR().nonZeros();
    }
    return true;
  }

protected:
  bool _init;
  SparseMatrix _sparseMatrix;
  sparseQRDecomposition _sparseQR;

  // compute the QR factorization
  bool computeSparseQR(const SparseBlockMatrix<MatrixType> &A, double &t) {
    // perform some operations only once upon init
    if (_init)
      _sparseMatrix.resize(A.rows(), A.cols());
    fillSparseMatrix(A, !_init);
    if (_init)
      computeSymbolicDecomposition();
    _init = false;

    t = get_monotonic_time();
    _sparseQR.factorize(_sparseMatrix.selfadjointView<Eigen::Upper>());
    if (_sparseQR.info() != Eigen::Success) {
      if (this->writeDebug()) {
        G2O_ERROR("Sparse QR failure, writing debug.txt (Hessian loadable by "
                  "Octave)");
        A.writeOctave("debug.txt");
      } else {
        G2O_DEBUG("QR failure");
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
    _sparseQR.analyzePattern(_sparseMatrix.selfadjointView<Eigen::Upper>());
    _init = false; // Mark as initialized

    G2OBatchStatistics *globalStats = G2OBatchStatistics::globalStats();
    if (globalStats)
      globalStats->timeSymbolicDecomposition = get_monotonic_time() - t;
  }

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

  /**
   * Implementation of the general parts for computing the inverse blocks of the
   * linear system matrix. Here we call a function to do the underlying
   * computation.
   */
  bool solveBlocks_impl(
      const SparseBlockMatrix<MatrixType> &A,
      [[maybe_unused]] std::function<void(MarginalCovarianceCholesky &)>
          compute) {
    // compute the QR factor
    double t;
    bool sparseQRState = computeSparseQR(A, t);
    if (!sparseQRState)
      return false;

    // book keeping statistics
    G2OBatchStatistics *globalStats = G2OBatchStatistics::globalStats();
    if (globalStats) {
      globalStats->choleskyNNZ = _sparseQR.matrixR().nonZeros();
    }
    return true;
  }
};

} // namespace g2o

#include "solver_eigen_qr.hpp"
#endif
