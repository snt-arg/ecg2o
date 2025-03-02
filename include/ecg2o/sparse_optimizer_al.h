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


#ifndef G2O_GRAPH_OPTIMIZER_AL_H_
#define G2O_GRAPH_OPTIMIZER_AL_H_

#include "g2o/core/g2o_core_api.h"
#include "vertex_lagrange_multiplier.h"
#include <iostream>
#include <memory>
#include <optimizable_graph.h>
#include <sparse_optimizer.h>
#include <vector>
#include "g2o/core/optimization_algorithm_with_hessian.h"
#include "g2o/core/sparse_optimizer.h"
#include "optimization_algorithm.h"
#include "g2o/stuff/timeutil.h"
#include "g2o/stuff/logger.h"
#include <Eigen/Core>
#include <iostream>
#include <vector>
#include "solver.h"



namespace g2o {

//forwad declaration
template <int D, typename E, typename... VertexTypes>
class BaseFixedSizedEdgeEq;
 
// Specialized optimizer inheriting from SparseOptimizer
class G2O_CORE_API SparseOptimizerAL : public SparseOptimizer{
  
     
   
  public:
  double epsilon_stop_threshold = 1e-4;
  double epsilon_eq_threshold = 1e-4;


  protected: 
  std::set<OptimizableGraph::Edge*> _edgeEqSet; // Set of equality edges
  std::unordered_map<void*, std::function<void()>> updateMultipliersEqFunctionMap;

  bool _ALTerminationFlag = false;
  int _max_num_inner_iterations = 1; // Maximum number of inner iterations for the AL algorithm
  double _rho_min = .5;  // Minimum value for the penalty parameter
  double _rho_max = 50000;  // Maximum value for the penalty parameter
  double _init_rho = 1;
  double _rho_upate_factor = 10;
  double _init_eq_lagrange_multiplier = 0;
  std::vector<std::shared_ptr<g2o::OptimizableGraph::Vertex>> _vEqLagrangeMultiplier; // Set of Lagrangian vertices


 public:
   SparseOptimizerAL() = default;  // Default constructor
   ~SparseOptimizerAL() = default;  // Virtual destructor 
  
  
   template <int D, typename E, typename... VertexTypes>
   bool addEdgeEq(BaseFixedSizedEdgeEq<D,E, VertexTypes...>* e);

   template <int D, typename E, typename... VertexTypes, typename EdgeType>
   void updateMultipliers(EdgeType& edge);

   template <int D, typename E, typename... VertexTypes>
   void constructQuadraticFormEq(BaseFixedSizedEdgeEq<D, E, VertexTypes...>& edge);
     
   // solver parameter


   // Setter function
   bool algSettings(std::vector<double>& settings){  
    int number_of_settings = 6;
    if (settings.size() != number_of_settings) { 
     std::cerr << "[Error] Invalid settings for the AL algorithm \n"
               << "Setting vector size must be <<" << number_of_settings << "\n"
               <<  " max_num_inner_iterations (int), rho_min, rho_max, init_rho, init_eq_lagrange_multiplier" << std::endl;
               return false;
    }
    _max_num_inner_iterations = (int)settings[0]; // Maximum number of inner iterations for the AL algorithm
    _rho_min = settings[1];  // Minimum value for the penalty parameter
    _rho_max = settings[2];  // Maximum value for the penalty parameter
    _init_rho = settings[3]; // Initial value for the penalty parameter
    _rho_upate_factor = settings[4];
    _init_eq_lagrange_multiplier = settings[5]; // Initial value for the Lagrangian vertex
    return true;
  
  }

   void setPenaltyParameter(double rho){ _init_rho = rho; };
   void setMaxNumInnerIterations(int numIterations) { _max_num_inner_iterations = numIterations; }
   void setVLagrangianInitial(double vInitial) {_init_eq_lagrange_multiplier = vInitial; }
   void resetVLagrangian(){ 
    for (auto& vertex : _vEqLagrangeMultiplier) {
      std::cout << "Resetting the Lagrangian vertex to the initial value" << std::endl;
       vertex->setToOrigin();    
     }
   }
  
   //getter for the solver parameters
    double getVLagrangianInitial();
    int getMaxNumInnerIterations() { return _max_num_inner_iterations; }
    double getPenaltyParameter();
    double getMultiplierUpdateFactor();

    virtual int optimize(int iterations, bool online = false) override {

      if (_ivMap.size() == 0) {
        G2O_WARN(
            "0 vertices to optimize, maybe forgot to call "
            "initializeOptimization()");
        return -1;
      }
     
      // reset _t to the initial value
    
      int cjIterations = 0;
      double cumTime = 0;
      bool innerLoopStop = false;
      bool outerLoopStop = false;
      bool ok = true;
      resetVLagrangian();
    
    
      ok = _algorithm->init(online);
      if (!ok) {
        G2O_ERROR("Error while initializing");
        return -1;
      }
    
      _batchStatistics.clear();
      if (_computeBatchStatistics) _batchStatistics.resize(iterations);
      
      OptimizationAlgorithm::SolverResult result = OptimizationAlgorithm::OK;
    
     
    
       while (!outerLoopStop) {
      
        innerLoopStop = false;
         
        int i;
        while (!innerLoopStop) {
          i = 0;
          i++; 
          preIteration(cjIterations);
              
          if (_computeBatchStatistics) {
            G2OBatchStatistics& cstat = _batchStatistics[cjIterations];
            G2OBatchStatistics::setGlobalStats(&cstat);
            cstat.iteration = cjIterations;
            cstat.numEdges = _activeEdges.size();
            cstat.numVertices = _activeVertices.size();
          }
    
          double ts = get_monotonic_time();
          result = _algorithm->solve(cjIterations, online);
          ok = (result == OptimizationAlgorithm::OK);
    
          bool errorComputed = false;
          if (_computeBatchStatistics) {
            computeActiveErrors();
            errorComputed = true;
            _batchStatistics[cjIterations].chi2 = activeRobustChi2();
            _batchStatistics[cjIterations].timeIteration = get_monotonic_time() - ts;
          }
    
          if (verbose()) {
            double dts = get_monotonic_time() - ts;
            cumTime += dts;
            if (!errorComputed) computeActiveErrors();
            std::cerr << "iteration= " << cjIterations << "\t chi2= " << FIXED(activeRobustChi2())
                << "\t time= " << dts << "\t cumTime= " << cumTime
                << "\t edges= " << _activeEdges.size();
            _algorithm->printVerbose(std::cerr);
            std::cerr << std::endl;
          }     
          ++cjIterations;
          postIteration(cjIterations);
          
          // termination criteria
          innerLoopStop = i >= _max_num_inner_iterations || i >=  iterations || !ok || verifyTermination(); // -1 is the defualt and here 10 time the defualt
          
        }
    
   
        //update the multiplier for Equality constraints
        for (auto& edge : _edgeEqSet) {      
          executeEqMultiplierUpdate(edge); 
        }
    
    
        if (result == OptimizationAlgorithm::Fail) {
        return 0;
        }
    
        if (cjIterations >= iterations ||  terminate()  ) { // check if the maximum number of iterations is reached
          break;
        }
        // check if the termination condition is satisfied
        outerLoopStop = verifyTermination() && verifyEqFeasibility(_edgeEqSet);
      
       
        
     
      } // End of inner loop
    
       return cjIterations; 
  
     }


    bool verifyTermination(){  
      OptimizationAlgorithmWithHessian* algorithm = static_cast<OptimizationAlgorithmWithHessian*>(_algorithm);
      const double* update = algorithm->solver().x();
      size_t xSize = algorithm->solver().vectorSize();
      Eigen::Map<const Eigen::VectorXd> updateVec(update, xSize);
      auto epsilon = epsilon_stop_threshold;     
      return updateVec.norm() < epsilon;
    }

    bool verifyEqFeasibility(const std::set<OptimizableGraph::Edge*>& edgeSet) {
      auto epsilon = epsilon_eq_threshold;    

      for (auto& edge : edgeSet) {
          edge->computeError();             
          const double* error = edge->errorData();   
          int edgeDimension = edge->dimension() /2 ; // divide by 2 to account for equality constraints excluding the lagrange multiplier
          for (int i = 0; i < edgeDimension; ++i) {
              if (std::abs(error[i]) > epsilon) {
                  return false;  // If any error exceeds epsilon, return false
              }
          }
      }
      
      return true; // Return true only if all edges satisfy the condition
  }


  // Call the stored function later
void executeEqMultiplierUpdate(void* edgePtr) {
  auto it = updateMultipliersEqFunctionMap.find(edgePtr);
  if (it != updateMultipliersEqFunctionMap.end()) {
      it->second(); // Calls updateMultipliers on the correct edge
  } else {
      std::cerr << "[Error] Edge not found in function map" << std::endl;
  }
}


};





template<int D, typename E, typename... VertexTypes>
bool SparseOptimizerAL::addEdgeEq(BaseFixedSizedEdgeEq<D, E, VertexTypes...>* e) {
      // create new vertex for the lagrangian
      auto vLagrangian = std::make_shared<VertexLagrangeMultiplier<D>>();
      //auto* nu = new VertexLagrangeMultiplier<D>();
      vLagrangian->setId(e->getVertexLagrangeMultiplierId());
      vLagrangian->setInitialValue(_init_eq_lagrange_multiplier);
      vLagrangian->setFixed(true);
      
      // Add the Nu vertex to the optimizer                
      this->addVertex(vLagrangian.get());

      //assign the vertex to the edge
      e->setVertex(e->vertices().size() -1, vLagrangian.get());
      
      bool eresult = OptimizableGraph::addEdge(e);
      if (!eresult) {
        std::cerr << "[Error] adding Eq edge to the optimizer" << std::endl;
        return false;
      } 
      
       Eigen::Matrix<double, D, 1> rho_bar = Eigen::Matrix<double, D, 1>::Constant(_init_rho);
 
      // In equality constraints, the information matrix 2D is used to store the constraint violation and the rho_bar
      // we initialize the the constraint violation to error and rho_bar to as initial value 
      // Assign first D diagonal entries from `error.segment(0, D)  
      e->setConstraintViolationPrev(e->error().segment(0, D).cwiseAbs());
      e->setConstraintViolationPrev(std::move(rho_bar), D);

      e->setRho(rho_bar);

      //set the ConstructQuadraticFormImpl function
      auto lambda = [this](BaseFixedSizedEdgeEq<D, E, VertexTypes...>& edge) {
      this->constructQuadraticFormEq(edge);
      };    
      e->setConstructQuadraticFormImpl(lambda);

      // Store the function in a map instead
      updateMultipliersEqFunctionMap[e] = [this, e]() {
      this->updateMultipliers<D, E, VertexTypes...>(*e); 
      };
 
     //add teh edge to the set of inequality set
     _edgeEqSet.insert(e);
     _vEqLagrangeMultiplier.push_back(vLagrangian);

     return true;
}
 



template <int D, typename E, typename... VertexTypes, typename EdgeType>
void SparseOptimizerAL::updateMultipliers(EdgeType& edge) {
    auto error = edge.error();
    auto rho = edge.rho();
    auto multiplier = edge.lagrangeMultiplier();
    Eigen::Matrix<double, D, 1> rho_bar; 
    auto constraints_violation_prev = edge.information().diagonal().segment(0, D);
    auto constraints_violation = error.segment(0, D);

        

    for (int i = 0; i < D; ++i) {
      double x_d = 0, x_i = 0;
      
      switch (2) {
        case 1:       
           constraints_violation = constraints_violation.cwiseAbs();           
           rho_bar = constraints_violation_prev.segment(D, D);
       
          if (constraints_violation_prev[i] > constraints_violation[i])
              x_d = 1 - constraints_violation[i] / constraints_violation_prev[i];

          if (constraints_violation[i] > constraints_violation_prev[i])
              x_i = 1 - constraints_violation_prev[i] / constraints_violation[i];          
                  rho[i] = rho_bar[i] + x_d * (_rho_max - rho_bar[i]);
                  rho_bar[i] = rho_bar[i] + x_d * (_rho_max - rho_bar[i]);            
                  edge.setConstraintViolationPrev(constraints_violation);
                  edge.setConstraintViolationPrev(std::move(rho_bar), D);  
          break;
        case 2: 
                rho[i] = _rho_upate_factor * rho[i];
                if (rho[i] < _rho_min) rho[i] = _rho_min;
                if (rho[i] > _rho_max) rho[i] = _rho_max;
                break;
               
          }

        //update mulitplier  
         multiplier[i] = multiplier[i] + 2 * rho[i] * error[i];
    }
       
    edge.setLagrangeMultiplier(multiplier);
    edge.setRho(rho);
      
}



// Implementation of the quadratic form construction
template <int D, typename E, typename... VertexTypes>
void SparseOptimizerAL::constructQuadraticFormEq(BaseFixedSizedEdgeEq<D, E, VertexTypes...>& edge) { 
  auto error = edge.error();
  auto gamma = edge.lagrangeMultiplier();
   auto rho = edge.rho();  
  Eigen::Matrix<double, 2 * D, 1> weightedError = Eigen::Matrix<double, 2 * D, 1>::Zero();
  Eigen::DiagonalMatrix<double, 2 * D> omega_matrix = Eigen::DiagonalMatrix<double, 2 * D>(Eigen::Matrix<double, 2 * D, 1>::Zero());
 for (int i = 0; i < D; ++i) {    
    omega_matrix.diagonal()[i] =    rho[i]; 
    weightedError[i] =  -(omega_matrix.diagonal()[i] * error[i]+  0.5 * gamma[i]);
  }
  static const std::size_t _nr_of_vertices = sizeof...(VertexTypes) +1; // add the lagrange multiplier vertex
  edge.constructQuadraticFormNs(omega_matrix, weightedError, std::make_index_sequence<_nr_of_vertices>());    
} 

} // namespace g2o

#endif  // G2O_GRAPH_OPTIMIZER_AL_H_
