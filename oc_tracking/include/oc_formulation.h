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


#ifndef OC_FORMULATION_H
#define OC_FORMULATION_H

#include <memory>
#include <optimizable_graph.h>
#include <vector>
#include "oc_parameters.h"
#include <fstream>


namespace g2o::oc {

template <typename Optimizer>
class OCFormulation {
public:
    OCFormulation(int N, std::shared_ptr<Optimizer> optimizer, std::shared_ptr<OCParameters> param);

    // Add vertices
    void addVertex(int offset, int k);
    void addVertices(int k);
    void setVertexScalarEstimate(int id, double value);
    double getVertexScalarEstimate(int id);
    std::vector<double> getResults(std::string filename ="output.txt");
    double getForceInput(int k);

    // Add edges
    void addEdgeCost_state(int k);
    void addEdgeCost_input(int k);
    void addEdgeDynamics_v_h(int k);
    void addEdges(int k);

    void setN(int N);
    int getN() const;
    void setVLagrangianInitial(double vInitial);
    double getVLagrangianInitial();

    // OC setup
    void setupOC();

    // Set the initial guess for the optimization
    bool setInitialGuess(bool update_all = false);  

    void computeOffsets();

protected:
    int _N;
    std::shared_ptr<Optimizer> _optimizer;
    std::shared_ptr<OCParameters> _param;
    std::vector<double> _results;//v_d[0,N], f[0,N-1], eqDynamics_v_h[0,N-1]
    double _vLagrangianInitial = 0.0; // Initial value for the Lagrangian vertex



    // Offsets for easy ID calculation
    int _offset_v_h;
    int _offset_f;  
    int _offset_eqDynamics_v_h; // for the Lagrangian variables of the equality constraints


    // Vertices
    std::vector<std::shared_ptr<g2o::OptimizableGraph::Vertex>> _vertices;
    std::vector<std::shared_ptr<g2o::OptimizableGraph::Edge>> _edges;
};


}  // namespace g2o::oc

// Include the implementation
#include "oc_formulation.hpp"

#endif  // OC_FORMULATION_H


