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