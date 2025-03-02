 // the implementation of the vertexes is in the header file

#ifndef VERTEX_OC_HPP_
#define VERTEX_OC_HPP_

//#include "mpc/mpc_vertices.h"

#include <Eigen/Core>
  
namespace g2o::oc {
    
// VertexOC class
template <int D>
void VertexOC<D>::setToOriginImpl() {
    if constexpr (D == 1) {
        this->_estimate = 0.0; // Single value for D = 1
    } else {
        this->_estimate.setZero(); // Vector of zeros for D > 1
    }
}

// Apply an increment (delta) to the state
template <int D>
void VertexOC<D>::oplusImpl(const double* update) {
    if constexpr (D == 1) {
        this->_estimate += update[0]; // Single value update
    } else {
        Eigen::Map<const Eigen::Matrix<double, D, 1>> delta(update);
        this->_estimate += delta; // Vector update
    }
}

// Write the vertex state to an output stream
template <int D>
bool VertexOC<D>::write(std::ostream& os) const {
    if constexpr (D == 1) {
        os << this->_estimate;
    } else {
        os << this->_estimate.transpose();
    }
    return os.good();
}

// Read the vertex state from an input stream
template <int D>
bool VertexOC<D>::read(std::istream& is) {
    if constexpr (D == 1) {
        is >> this->_estimate;
    } else {
        for (int i = 0; i < D; ++i) {
            is >> this->_estimate[i];
        }
    }
    return is.good();
}


}

#endif  // VERTEX_OC_HPP_