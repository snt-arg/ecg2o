#ifndef MPC_FORMULATION_HPP
#define MPC_FORMULATION_HPP

#include "oc_vertices.h"
#include "oc_formulation.h"
#include <optimization_algorithm_property.h>
#include <ostream>
#include <stdexcept>
#include <memory>
#include <optimizable_graph.h>
#include <unistd.h>
#include <valarray>
#include <vector>
#include "sparse_optimizer.h"
#include "oc_vertices.h"
#include "oc_edges.h"
#include "g2o/core/base_fixed_sized_edge.h"

#include "ecg2o/sparse_optimizer_eq.h" // for utilizing the factor graph solver
#include "ecg2o/sparse_optimizer_al.h" // using the AL solver

  
namespace g2o::oc {


template <typename Optimizer>
OCFormulation<Optimizer>::OCFormulation(int N, std::shared_ptr<Optimizer> optimizer, std::shared_ptr<OCParameters> param)
    : _N(N), _optimizer(optimizer), _param(param) {
    computeOffsets();
}

template <typename Optimizer>
void OCFormulation<Optimizer>::computeOffsets() {
    _offset_v_h = 0 * (_N + 1);
    _offset_f = 1 * (_N + 1);
    _offset_eqDynamics_v_h = 2 * (_N + 1);
}


template <typename Optimizer>
void OCFormulation<Optimizer>::setN(int N) { _N = N; };   

template <typename Optimizer>
int OCFormulation<Optimizer>::getN() const { return _N; };


template <typename Optimizer>
void OCFormulation<Optimizer>::setVLagrangianInitial(double vInitial) {_vLagrangianInitial = vInitial; }

template <typename Optimizer>
double OCFormulation<Optimizer>::getVLagrangianInitial() { return _vLagrangianInitial; }



template <typename Optimizer>
void OCFormulation<Optimizer>::addVertex(int offset, int k) {
        auto vertex = std::make_shared<Vertex_scalar>();
        vertex->setId(offset + k);
        //vertex->setId(offset + num_var_types * k); another way to set the id by grouping the varaible in based on k
        // if we have k continously increasingly num_var_types =7 and offset = 0,1, ... not 0*N, 1*N, ...
        _optimizer->addVertex(vertex.get());
        _vertices.push_back(vertex);
    }
// add vertices

template <typename Optimizer>
void OCFormulation<Optimizer>::addVertices(int k) {
    addVertex(_offset_v_h, k);
    if (k ==0)  {
       _optimizer->vertex(_offset_v_h + k)->setFixed(true);
    }

 
    if (k <= _N - 1) { // x_{N-1} is the last state
      addVertex(_offset_f, k);
    }
}


template <typename Optimizer>
void OCFormulation<Optimizer>::setVertexScalarEstimate(int id, double value) {
    static_cast<Vertex_scalar*>(_optimizer->vertex(id))->setEstimate(value);
}


template <typename Optimizer>
double OCFormulation<Optimizer>::getVertexScalarEstimate(int id) {
    try {
        auto vertex = _optimizer->vertex(id); // Direct call to optimizer
        if (!vertex) {
            throw std::runtime_error("Vertex with ID " + std::to_string(id) + " not found.");
        }
        return static_cast<Vertex_scalar*>(vertex)->estimate();
    } catch (const std::runtime_error& e) {
        std::cerr << "Error in getVertexScalarEstimate: " << e.what() << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
  
 }


// get the results
template <typename Optimizer>
std::vector<double> OCFormulation<Optimizer>::getResults() {
    _results.clear();
    for (int k = 0; k <= _N; k++) {
        _results.push_back(getVertexScalarEstimate(_offset_v_h + k));
    } 
    _results.push_back(101010101);

    for (int k = 0; k <= _N - 1; k++) {
        _results.push_back(getVertexScalarEstimate(_offset_f + k));
    }
    

    _results.push_back(101010101);
   

    for (int k = 0; k <= _N - 1; k++) {
      _results.push_back(getVertexScalarEstimate(_offset_eqDynamics_v_h + k));
    }


    std::ofstream outFile("output.txt"); // Open file for writing

    if (outFile.is_open()) {
        outFile << "[ "; // Start with an opening bracket
        for (size_t i = 0; i < _results.size(); ++i) {
            outFile << _results[i];
            if (i != _results.size() - 1) outFile << ", "; // Add comma except for the last element
        }
        outFile << " ]\n"; // End with a closing bracket

        outFile.close(); // Close the file
        std::cout << "Vector saved to output.txt with bracket format.\n";
    } else {
        std::cerr << "Error opening file!\n";
    }
        

    return _results;
}
// get the force input
template <typename Optimizer>
double OCFormulation<Optimizer>::getForceInput(int k) {
    return getVertexScalarEstimate(_offset_f + k);
}




template <typename Optimizer>
void OCFormulation<Optimizer>::addEdgeCost_state(int k) {
    auto edge = std::make_shared<EdgeCost_state>(k, _param);
    edge->setVertex(0, _optimizer->vertex(_offset_v_h + k));
    Eigen::Matrix<double, 1, 1> information;
    if (k == _N) {
        information << _param->get_P();
     } else {
        information << _param->get_Q();
    }
    edge->setInformation(information);
    _optimizer->addEdge(edge.get());
    _edges.push_back(edge);
}

template <typename Optimizer>
void OCFormulation<Optimizer>::addEdgeCost_input(int k) {
    auto edge = std::make_shared<EdgeCost_input>(k, _param);
    edge->setVertex(0, _optimizer->vertex(_offset_f + k));
    Eigen::Matrix<double, 1, 1> information;
    information << _param->get_R();
    edge->setInformation(information);
    _optimizer->addEdge(edge.get());
    _edges.push_back(edge);
}

template <typename Optimizer>
void OCFormulation<Optimizer>::addEdgeDynamics_v_h(int k) {
    auto edge = std::make_shared<EdgeDynamics_v_h>(k, _param);
    edge->setVertex(0, _optimizer->vertex(_offset_v_h + k));
    edge->setVertex(1, _optimizer->vertex(_offset_v_h + k + 1));
    edge->setVertex(2, _optimizer->vertex(_offset_f + k));
     edge->setVertexLagrangeMultiplierId(_offset_eqDynamics_v_h + k);
    _optimizer->addEdgeEq(edge.get());
    _edges.push_back(edge);
}

template <typename Optimizer>
void OCFormulation<Optimizer>::addEdges(int k){
    std::cout << "Adding edges for k = " << k << std::endl;

    //addEdgeCost_input(k);
    if (k > 0 ) { // v_h_0 is fixed
      addEdgeCost_state(k);
    }    
    
    if (k <= _N - 1) {
        addEdgeCost_input(k);
        addEdgeDynamics_v_h(k);
    }
 

     
  
}

template <typename Optimizer>
bool OCFormulation<Optimizer>::setInitialGuess(bool update_all) {
    // set the initial guess for the optimization

    // uddate v_h_0 and d_h_0
    setVertexScalarEstimate(_offset_v_h, _param->get_initial_v_h_0());
   
    if (update_all) {

        for (int k = 1; k <= _N; k++) {
            setVertexScalarEstimate(_offset_v_h + k, _param->get_initial_v_h());
          }

        for (int k = 0; k <= _N - 1; k++) {
            setVertexScalarEstimate(_offset_f + k, _param->get_initial_f());        

        }

      _optimizer->resetVLagrangian(); // automatically done by the optimizer but we do it to see the initial values
    } 
     
     // print the resutls 
    std::vector<double> results = this->getResults();
    bool print_results = false;

    if (print_results) {
        std::cout << "Initial guess set.\nResults: f_t_0, f_b_0, ..., slack_1, slack_2, slack_3,v_h, d_h \n"
                  << Eigen::Map<Eigen::VectorXd>(results.data(), results.size()).transpose() << std::endl;
    }
   
    
    return true;
}

template <typename Optimizer>
void OCFormulation<Optimizer>::setupOC() {
    std::cout << "Setting up MPC formulation with Horizon "<< _N << std::endl;
    

     

    // creeate verteices
        for (int k = 0; k <= _N; k++) {
            addVertices(k);
         }



    // create edges
        for (int k = 0; k <= _N; k++) {
            addEdges(k);
        }

    // set the initial guess for the optimization
        setInitialGuess(true);    
    
}

 
}
#endif // MPC_FORMULATION_HPP