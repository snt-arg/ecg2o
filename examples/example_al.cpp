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
     minimize     (xy(1)-2)^2 +  (xy(2)-9)^2 + (z-50)^2
     subject to
                      xy(2) + z    == 3, 
                      xy(1) + xy(2)    == 2.5,
%  Note that x = xy(0) and y = xy(1)               
% The answer is x = 15, y = -12.5, z = 15.5

OP: min ||xy(1)-a||^2 +  ||xy(2)-b||^2 + ||z-c||^2 s.t.    xy(2) + z == 3, xy(1) + xy(2) == 2.5
*/


#include <g2o/core/block_solver.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/core/optimization_algorithm_levenberg.h>
#include <g2o/core/optimization_algorithm_dogleg.h>
#include <g2o/solvers/cholmod/linear_solver_cholmod.h>

#include <Eigen/Core>
#include <iostream>


#include "ecg2o/sparse_optimizer_eq.h" // for utilizing the factor graph solver
#include "ecg2o/sparse_optimizer_al.h" // using the AL solver
 
#include "include/example_vertices_edges.h" // for defining the edges

 
 

 
 


int main(int argc, char** argv) {
    std::cout << "OP: min ||xy(1)-a||^2 +  ||xy(2)-b||^2 + ||z-c||^2 s.t.    xy(2) + z == 3, xy(1) + xy(2) == 2.5," << std::endl;
    std::cout << "Usage: " << argv[0] << "solverType:<0=GN,1=LV,2=DL> iterations<int> a<int> b<int> c<int>"<< std::endl;
    std::cout << "x_start<int> y_start<int> z_start<int> " << std::endl;
    std::cout << "setting <{[0] max_num_inner_iterations, [1] rho_min, [2] rho_max, [3] init_rho, [4] rho_upate_factor, [5] init_eq_lagrange_multiplier}>" << std::endl;
    std::cout << "./example_al 0 150 2 9 50 10 10 10 1 0.5 5 1 1.2 0"<<std::endl;

    
    int argCount = 1;
    int solverType = (argc > argCount) ? std::atoi(argv[argCount]) : 0;
    
    argCount++;
    int numberOfIterations = (argc > argCount) ? std::atoi(argv[argCount]) : 150;
    
    argCount++;
    int a = (argc > argCount) ? std::atoi(argv[argCount]) : 2;
    
    argCount++;
    int b = (argc > argCount) ? std::atoi(argv[argCount]) : 9;
    
    argCount++;
    int c = (argc > argCount) ? std::atoi(argv[argCount]) : 50;
    
    argCount++;
    int xStart = (argc > argCount) ? std::atoi(argv[argCount]) : 10;
    
    argCount++;
    int yStart = (argc > argCount) ? std::atoi(argv[argCount]) : 10;
    
    argCount++;
    int zStart = (argc > argCount) ? std::atoi(argv[argCount]) : 10;
   
    // [0] max_num_inner_iterations, [1] rho_min, [2] rho_max, [3] init_rho, [4] rho_upate_factor, [5] init_eq_lagrange_multiplier
    std::vector<double> settings = {1, 0.5, 5, 1, 1.4, 0};
    
    argCount++;
    if (argc > argCount) settings[0] = std::atoi(argv[argCount]);

    argCount++;
    if (argc > argCount) settings[1] = std::atof(argv[argCount]);

    argCount++;
    if (argc > argCount) settings[2] = std::atof(argv[argCount]);

    argCount++;
    if (argc > argCount) settings[3] = std::atof(argv[argCount]);

    argCount++;
    if (argc > argCount) settings[4] = std::atof(argv[argCount]);

    argCount++;
    if (argc > argCount) settings[5] = std::atof(argv[argCount]);



    // Initialize optimizer
    g2o::SparseOptimizerAL optimizer;
    optimizer.setVerbose(true);

    std::unique_ptr<g2o::BlockSolverX> blockSolver;
    std::unique_ptr<g2o::LinearSolver<g2o::BlockSolverX::PoseMatrixType>> linearSolver;
    std::unique_ptr<g2o::OptimizationAlgorithm> algorithm;

  

    std::cout << "Using CHOLMOD linear solver" << std::endl;
    linearSolver = std::make_unique<g2o::LinearSolverCholmod<g2o::BlockSolverX::PoseMatrixType>>();
  

    blockSolver = std::make_unique<g2o::BlockSolverX>(std::move(linearSolver));

    // Set up the optimization algorithm
    switch (solverType) {
        case 1:
            std::cout << "Using Levenberg-Marquardt solver" << std::endl;
            algorithm = std::make_unique<g2o::OptimizationAlgorithmLevenberg>(std::move(blockSolver));
            break;
        case 2:
            std::cout << "Using Dogleg solver" << std::endl;
            algorithm = std::make_unique<g2o::OptimizationAlgorithmDogleg>(std::move(blockSolver));
            break;
        default:
            std::cout << "Using Gauss-Newton solver" << std::endl;
            algorithm = std::make_unique<g2o::OptimizationAlgorithmGaussNewton>(std::move(blockSolver));
            break;
    }

    optimizer.setAlgorithm(algorithm.get()); 

    //auto* nu = new g2o::mpc::VertexMPC<2>();
    
    
    // Add the vertex to the optimizer
    auto* xy = new VertexXY();
    xy->setId(0);
    xy->setEstimate(Eigen::Vector2d(xStart, yStart));  // initial estimate
    optimizer.addVertex(xy);

    // Create and add vertices
    auto* z = new VertexZ();
    z->setId(1);
    z->setEstimate(zStart);  // initial estimate
    optimizer.addVertex(z);
    
    // Add the EdgeZ to the optimizer
    auto edgeXY = std::make_shared<EdgeXY>();
    edgeXY->setVertex(0, xy);
    edgeXY->setMeasurement(Eigen::Vector2d(a, b));  // target value (2, 9)
    edgeXY->setInformation(Eigen::Matrix2d::Identity());
    optimizer.addEdge(edgeXY.get());
 

    auto edgeZ = std::make_shared<EdgeZ>(); // (z- a)^2  
    edgeZ->setVertex(0, z);
    edgeZ->setMeasurement(c);  // target value 50
    edgeZ->setInformation(Eigen::Matrix<double, 1, 1>::Identity());
    optimizer.addEdge(edgeZ.get());     
 



    auto* edgeEq = new  EdgeEq();  
    edgeEq->setVertexLagrangeMultiplierId(10);
    edgeEq->setVertex(0, xy);
    edgeEq->setVertex(1, z);
    optimizer.addEdgeEq(edgeEq);



 
    // Optimize
    optimizer.initializeOptimization();
    
  
    optimizer.algSettings(settings); 
    optimizer.optimize(numberOfIterations);

    // Output the results
    std::cout << "Optimized xy: " << xy->estimate().transpose() << std::endl;
    std::cout << "Optimized z: " << z->estimate() << std::endl;
    std::cout << "Optimized Lagrange Multiplier: " << edgeEq->lagrangeMultiplier().transpose() << std::endl;


    //std::cout << "Optimized nu: " << nu->estimate() << std::endl;

    std::cout << "Settings: ";
    for (const auto& setting : settings) {
        std::cout << setting << " ";
    }
    std::cout << std::endl;
    
    std::cout << "-------- al " << std::endl;
    return 0;

}
