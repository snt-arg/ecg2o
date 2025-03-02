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
    std::vector<double> getResults();
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
    std::vector<double> _results;//v_d[0,N], d_h[0,N], f_t[0,N-1], f_b[0,N-1], slack_1[1,N], slack_2[0,N-1], slack_3[0,N-1]
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


