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

#ifndef OS_SYS_DYNAMICS_SIM_H
#define OS_SYS_DYNAMICS_SIM_H
 

#include "oc_parameters.h"
#include <memory>
#include <cmath>

namespace g2o::oc {

class system_dynamics_simulator {
public:
 
    system_dynamics_simulator(std::shared_ptr<OCParameters> param) : _param(param) {};
  
    double simulate(double v_h_prev, double u_car) {
    
    double alpha = 0.0; // Gradient in radians, assumed flat

    double Ts = _param->get_Ts();
    double g = _param->get_g();
    double m_total = _param->get_m_v() + _param->get_m_p();
    double m_eq = m_total + _param->get_e_i() * _param->get_m_v();
    double c_r = _param->get_x1() * _param->get_c_r();
    double c_a = _param->get_x2() * _param->get_c_a();
    
    double f_roll = (v_h_prev > 0.05) ? (c_r * m_total * g * cos(alpha)) : 0.0;
    
    // Calculate new velocity
    double v_h_1 = (Ts / m_eq) * (u_car - m_total * g * sin(alpha) - f_roll - c_a * v_h_prev * v_h_prev) + v_h_prev;
    
    // Ensure velocity is non-negative
    return std::max(v_h_1, 0.0);
    }

private:   
    std::shared_ptr<OCParameters> _param;
   
};

}  // namespace g2o::oc

#endif  // EDGE_COST_ENERGY_H
