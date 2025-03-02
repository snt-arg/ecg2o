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
    int numofRuns = (argc > argCount) ? std::atoi(argv[argCount]) : 1; // v_lagrange_initial
    argCount++;
    int numberOfIterations = (argc > argCount) ? std::atoi(argv[argCount]) : 500;
    argCount++;
    int linearSolverType = (argc > argCount) ? std::atoi(argv[argCount]) :13;
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
    #endif
    
    optimizer->setVerbose(true);
    optimizer->setAlgorithm(algorithm.get());
    optimizer->setVLagrangianInitial(0); // the initial value of the Lagrangian vertex

          
 
    // Update parameters
    param->set_initial_v_h_0(0);
    param->set_initial_d_h_0(10.2);
    param->set_TracForce_Prev(111.0);
    param->set_BrakeForce_Prev(222.0);
    //std::vector<double> vp_prediction = {0.5, 1, 1.5, 1.7, 1.9, 2.1,2.2,2.3,2.4,2.5 ,2.6};
    //param->set_vp_prediction(vp_prediction);
   oc->setupOC();

    // Set the initial guess
    oc->setInitialGuess(true);

    std::vector<double> initial = oc->getResults();
    std::cout << "Initial Guess:\n"
              << Eigen::Map<Eigen::VectorXd>(initial.data(), initial.size()).transpose() << std::endl;


    optimizer->initializeOptimization();

    auto start_time = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < numofRuns; i++){
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
