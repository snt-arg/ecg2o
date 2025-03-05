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
     minimize     ||x_N + r_N ||^2_P + sum{k=1 .. N-1}(||x_k + r_k ||^2_Q +||u_k||^2_R)
     subject to     x_{k+1} - (x_k + delta_t/m(u_k -F_{resist})) = 0 , k= 0.. N-1
                       

    ||e||^2_X= e'Xe 
*/




#include "oc_parameters.h"
#include "oc_formulation.h"
#include <Eigen/src/Core/products/Parallelizer.h>
#include <g2o/core/base_fixed_sized_edge.h>
#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/optimization_algorithm_dogleg.h>
#include <g2o/solvers/cholmod/linear_solver_cholmod.h>


#include <Eigen/Core>
#include <iostream>
#include <memory>
#include <vector>

#include "ecg2o/sparse_optimizer_eq.h" // for utilizing the factor graph solver
#include "ecg2o/sparse_optimizer_al.h" // using the AL solver
 


 
 

int main(int argc, char** argv) {
       std::cout << "MPC-based ACC" << std::endl;
    std::cout << "Usage: " << argv[0] << std::endl;
  

    // Parse command-line arguments
    int argCount = 1;
    int Horizon = (argc > argCount) ? std::atoi(argv[argCount]) : 385; //  horizon length
    argCount++;
    int numofRuns = (argc > argCount) ? std::atoi(argv[argCount]) : 1; 
    argCount++;
    int numberOfIterations = (argc > argCount) ? std::atoi(argv[argCount]) : 50;
    argCount++;
    int solverType = (argc > argCount) ? std::atoi(argv[argCount]) :0;

    std::unique_ptr<g2o::BlockSolverX> blockSolver;
    std::unique_ptr<g2o::LinearSolver<g2o::BlockSolverX::PoseMatrixType>> linearSolver;
    std::unique_ptr<g2o::OptimizationAlgorithm> algorithm;


    std::cout << "Using CHOLMOD linear solver" << std::endl;
    linearSolver = std::make_unique<g2o::LinearSolverCholmod<g2o::BlockSolverX::PoseMatrixType>>();
  
    blockSolver = std::make_unique<g2o::BlockSolverX>(std::move(linearSolver));
 
    // Set up the optimization algorithm
    switch (solverType) {
        case 1:
            algorithm = std::make_unique<g2o::OptimizationAlgorithmLevenberg>(std::move(blockSolver));
            break;
        case 2:
            algorithm = std::make_unique<g2o::OptimizationAlgorithmDogleg>(std::move(blockSolver));
            break;
        default:
            algorithm = std::make_unique<g2o::OptimizationAlgorithmGaussNewton>(std::move(blockSolver));
            break;
    }

    // Initialize MPC parameters
    auto param = std::make_shared<g2o::oc::OCParameters>();

//---------------------------------------------------------------------------------------------

#define USE_Eq
    // Initialize optimizer
     #ifdef USE_Eq
    auto optimizer = std::make_shared<g2o::SparseOptimizerEq>();
    auto oc = std::make_shared<g2o::oc::OCFormulation<g2o::SparseOptimizerEq>>(Horizon, optimizer, param);
    #else
    auto optimizer = std::make_shared<g2o::SparseOptimizerAL>();      
    auto oc = std::make_shared<g2o::oc::OCFormulation<g2o::SparseOptimizerAL>>(Horizon, optimizer, param);
    std::vector<double> settings = {1, 0.5, 50000, 1, 10, 1};     // [0] max_num_inner_iterations, [1] rho_min, [2] rho_max, [3] init_rho, [4] rho_upate_factor, [5] init_eq_lagrange_multiplier
    optimizer->algSettings(settings);   
    #endif
    


    optimizer->setVerbose(true);
    optimizer->setAlgorithm(algorithm.get());
    optimizer->setVLagrangianInitial(0); // the initial value of the Lagrangian vertex

          
    oc->setupOC();
 
    // Update parameters
    param->set_initial_v_h_0(0);
    param->set_initial_f(0);
    param->set_linearDynamics(false);

    // Set the initial guess
    oc->setInitialGuess(true);

    std::vector<double> initial = oc->getResults();
    std::cout << "Initial Guess:\n"
              << Eigen::Map<Eigen::VectorXd>(initial.data(), initial.size()).transpose() << std::endl;


    optimizer->initializeOptimization();

    auto start_time = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < numofRuns; i++){
          oc->setInitialGuess(true);
          optimizer->optimize(numberOfIterations);
    }
 

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = end_time - start_time;




    // Retrieve and print results
    std::vector<double> results = oc->getResults();
    double u_input = oc->getForceInput(0);

   // std::cout << "Optimization done.\nResults:\n"
          //    << Eigen::Map<Eigen::VectorXd>(results.data(), results.size()).transpose() << std::endl;
    std::cout << "u_input: " << u_input << std::endl;
    
    std::cout << "Elapsed time: " << elapsed_seconds.count()/numofRuns << "s\n";

    return 0;
    
}
