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
