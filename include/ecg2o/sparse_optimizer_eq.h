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



#ifndef G2O_GRAPH_OPTIMIZER_EQ_H
#define G2O_GRAPH_OPTIMIZER_EQ_H

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
template <int D, typename E, typename... VertexTypes>
class BaseFixedSizedEdgeEq;



// Specialized optimizer inheriting from SparseOptimizer
class G2O_CORE_API SparseOptimizerEq : public SparseOptimizer{

 public:
   SparseOptimizerEq() = default;;  // Default constructor
   ~SparseOptimizerEq() = default;  // Virtual destructor 

   template <int D, typename E, typename... VertexTypes>
   bool addEdgeEq(BaseFixedSizedEdgeEq<D,E, VertexTypes...>* e);

   int optimize(int iterations, bool online = false) override {
    if (_ivMap.size() == 0) {
      G2O_WARN(
          "0 vertices to optimize, maybe forgot to call "
          "initializeOptimization()");
      return -1;
    }
  
    int cjIterations = 0;
    double cumTime = 0;
    bool ok = true;
    bool stop = false;
  
    ok = _algorithm->init(online);
    if (!ok) {
      G2O_ERROR("Error while initializing");
      return -1;
    }
    resetVLagrangian();

  
    _batchStatistics.clear();
    if (_computeBatchStatistics) _batchStatistics.resize(iterations);
  
    OptimizationAlgorithm::SolverResult result = OptimizationAlgorithm::OK;
    for (int i = 0; i < iterations && !terminate() && !stop && ok; i++) {
      preIteration(i);
  
      if (_computeBatchStatistics) {
        G2OBatchStatistics& cstat = _batchStatistics[i];
        G2OBatchStatistics::setGlobalStats(&cstat);
        cstat.iteration = i;
        cstat.numEdges = _activeEdges.size();
        cstat.numVertices = _activeVertices.size();
      }
  
      double ts = get_monotonic_time();
      result = _algorithm->solve(i, online);
      ok = (result == OptimizationAlgorithm::OK);
  
      bool errorComputed = false;
      if (_computeBatchStatistics) {
        computeActiveErrors();
        errorComputed = true;
        _batchStatistics[i].chi2 = activeRobustChi2();
        _batchStatistics[i].timeIteration = get_monotonic_time() - ts;
      }
  
      if (verbose()) {
        double dts = get_monotonic_time() - ts;
        cumTime += dts;
        if (!errorComputed) computeActiveErrors();
        std::cerr << "iteration= " << i << "\t chi2= " << FIXED(activeRobustChi2())
             << "\t time= " << dts << "\t cumTime= " << cumTime
             << "\t edges= " << _activeEdges.size();
        _algorithm->printVerbose(std::cerr);
        std::cerr << std::endl;
      }
      ++cjIterations;
      postIteration(i);
      stop = verifyTermination();
    }
    if (result == OptimizationAlgorithm::Fail) {
      return 0;
    }
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

  void resetVLagrangian() { 
    for (auto& vertex : _vEqLagrangeMultiplier) {
      std::cout << "Resetting the Lagrangian vertex to the initial value" << std::endl;
       vertex->setToOrigin();    
    }
  }

  public:
    void setVLagrangianInitial(double vInitial) {_init_eq_lagrange_multiplier = vInitial; }
  
    double _init_eq_lagrange_multiplier = 0.0; // Initial value for the Lagrangian vertex
    double epsilon_stop_threshold = 1e-3;

    protected: 
    std::vector<std::shared_ptr<g2o::OptimizableGraph::Vertex>> _vEqLagrangeMultiplier; // Set of Lagrangian vertices
    bool _TerminationFlag;
    
    
  };
  template <int D, typename E, typename... VertexTypes>
  bool SparseOptimizerEq::addEdgeEq(BaseFixedSizedEdgeEq<D, E, VertexTypes...>* e){
    // create new vertex for the lagrangian
      auto vLagrangian = std::make_shared<VertexLagrangeMultiplier<D>>();
      //auto* nu = new VertexLagrangeMultiplier<D>();
      vLagrangian->setId(e->getVertexLagrangeMultiplierId());
      vLagrangian->setInitialValue(_init_eq_lagrange_multiplier);
      
    // Add the Nu vertex to the optimizer                
      this->addVertex(vLagrangian.get());

      //assign the vertex to the edge
      e->setVertex(e->vertices().size() -1, vLagrangian.get());
      e->initializeInformationMatrix();
      
      bool eresult = OptimizableGraph::addEdge(e);
      if (!eresult) {
        std::cerr << "[Error] adding Eq edge to the optimizer" << std::endl;
        return false;
      } 

    //add teh edge to the set of inequality set
    _vEqLagrangeMultiplier.push_back(vLagrangian);

    return true;


}








}  // namespace g2o

#endif  // G2O_GRAPH_OPTIMIZER_EQ_H_
