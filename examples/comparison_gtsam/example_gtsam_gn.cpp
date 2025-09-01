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

/*
     minimize     0.5*||x1 + e^(-x2)||^2 + 0.5*||x1^2 + 2*x2 + 1||^2
    subject to    x1 + x1^3 + x2 + x2^2  = 0

*/

#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_dogleg.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/solvers/cholmod/linear_solver_cholmod.h>

#include <Eigen/Core>
#include <iostream>

#include "ecg2o/linear_solver_eigen_lu.h" // for using the Eigen LU solver
#include "ecg2o/linear_solver_eigen_qr.h" // for using the Eigen QR solver
#include "ecg2o/sparse_optimizer_eq.h" // for utilizing the factor graph solver
#include "example_gtsam_vertices_edges.h" // for defining the edges

int main(int argc, char **argv) {
  std::cout
      << "OP: min  0.5*||x1 + e^(-x2)||^2 + 0.5*||x1^2 + 2*x2 + 1||^2 s.t.    "
         " x1 + x1^3 + x2 + x2^2  = 0"
      << std::endl;
  std::cout << "Usage: " << argv[0]
            << "x_start<int> y_start<int> solverType:<0=GN,1=LV,2=DL> "
               "linearSolverType:<0=CHOLMOD,1=QR,2=LU> iterations<int> "
            << std::endl;
  std::cout << " " << std::endl;
  std::cout << "./example -0.2 -0.2 0 0 150" << std::endl;

  std::cout << " KKT-based GN method is Used\n";

  int argCount = 1;
  double xStart = (argc > argCount) ? std::atof(argv[argCount]) : -1.2;

  argCount++;
  double yStart = (argc > argCount) ? std::atof(argv[argCount]) : -1;
  argCount++;
  int solverType = (argc > argCount) ? std::atoi(argv[argCount]) : 0;

  argCount++;
  int linearSolverType = (argc > argCount) ? std::atoi(argv[argCount]) : 0;

  argCount++;
  int numberOfIterations = (argc > argCount) ? std::atoi(argv[argCount]) : 150;

  // Initialize optimizer
  std::cout << " KKT-based GN method is Used\n";
  g2o::SparseOptimizerEq optimizer;
  optimizer.setVerbose(true);

  std::unique_ptr<g2o::BlockSolverX> blockSolver;
  std::unique_ptr<g2o::LinearSolver<g2o::BlockSolverX::PoseMatrixType>>
      linearSolver;
  std::unique_ptr<g2o::OptimizationAlgorithm> algorithm;

  switch (linearSolverType) {
  case 1:
    std::cout << "Using Eigen QR linear solver" << std::endl;
    linearSolver = std::make_unique<
        g2o::LinearSolverEigenQR<g2o::BlockSolverX::PoseMatrixType>>();
    break;
  case 2:
    std::cout << "Using Eigen LU linear solver" << std::endl;
    linearSolver = std::make_unique<
        g2o::LinearSolverEigenLU<g2o::BlockSolverX::PoseMatrixType>>();
    break;

  default:
    std::cout << "Using CHOLMOD linear solver" << std::endl;
    linearSolver = std::make_unique<
        g2o::LinearSolverCholmod<g2o::BlockSolverX::PoseMatrixType>>();
    break;
  }

  blockSolver = std::make_unique<g2o::BlockSolverX>(std::move(linearSolver));

  // Set up the optimization algorithm
  switch (solverType) {
  case 1:
    std::cout << "Using Levenberg-Marquardt solver" << std::endl;
    algorithm = std::make_unique<g2o::OptimizationAlgorithmLevenberg>(
        std::move(blockSolver));
    break;
  case 2:
    std::cout << "Using Dogleg solver" << std::endl;
    algorithm = std::make_unique<g2o::OptimizationAlgorithmDogleg>(
        std::move(blockSolver));
    break;
  default:
    std::cout << "Using Gauss-Newton solver" << std::endl;
    algorithm = std::make_unique<g2o::OptimizationAlgorithmGaussNewton>(
        std::move(blockSolver));
    break;
  }

  optimizer.setAlgorithm(algorithm.get());

  // Add the vertex to the optimizer
  auto xy = std::make_shared<VertexXY>();
  xy->setId(0);
  xy->setEstimate(Eigen::Vector2d(xStart, yStart)); // initial estimate
  optimizer.addVertex(xy.get());

  // Add the cost edgeXY to the optimizer: 0.5*||x1 + e^(-x2)||^2 + 0.5*||x1^2 +
  // 2*x2 + 1||^2
  auto edgeXY = std::make_shared<EdgeXY>();
  edgeXY->setVertex(0, xy.get());
  edgeXY->setMeasurement(Eigen::Vector2d(0, 0)); //
  edgeXY->setInformation(Eigen::Matrix2d::Identity() * 0.5);
  optimizer.addEdge(edgeXY.get());

  auto *edgeEq = new EdgeEq(); // add x1 + x1^3 + x2 + x2^2 = 0
  edgeEq->setVertexLagrangeMultiplierId(10);
  edgeEq->setVertex(0, xy.get());
  optimizer.addEdgeEq(edgeEq);

  // Optimize
  optimizer.initializeOptimization();
  class terminationCriterionType;

  optimizer.epsilon_stop_threshold = 1e-3; // Stopping criterion for update norm

  optimizer.optimize(numberOfIterations);

  // Output the results
  std::cout << "Optimized x: " << xy->estimate().transpose() << std::endl;

  // std::cout << "Optimized nu: " << nu->estimate() << std::endl;

  return 0;
}
