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


#ifndef EDGE_COST_INPUT_HPP
#define EDGE_COST_INPUT_HPP

 #include <Eigen/Dense>

namespace g2o::oc {

EdgeCost_input::EdgeCost_input(int k, std::shared_ptr<OCParameters> param = nullptr)
    : _k(k), _param(param){}

void EdgeCost_input::computeError() {
    const Vertex_f* v_f = static_cast<const Vertex_f*>(this->_vertices[0]);

    double f = v_f->estimate();
 
    _error[0] = f;

}

void EdgeCost_input::linearizeOplus() {
    // The derivative of _error[0] with respect to the variable f is 1.0
     // Set the Jacobian's value for the scalar relationship
     
     
    auto& J_v0 = std::get<0>(this->_jacobianOplus);  
    J_v0 << 1.0;   // J(e_0,v_0[0])

}


bool EdgeCost_input::write(std::ostream& os) const {
    os << "EdgeCost_brake_change, ";
    for (int i = 0; i < 1; i++) {
        os << _error[i] << " ";
    }

    return os.good();

}

bool EdgeCost_input::read(std::istream& is) {
    for (int i = 0; i < 1; i++) {
        is >> _error[i];
    }

    return is.good();
    
}

 

}  // namespace g2o::oc

#endif  // EDGE_COST_INPUT_HPP

