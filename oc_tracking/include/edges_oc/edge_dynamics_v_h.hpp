#ifndef EDGE_DYNAMICS_V_H_HPP
#define EDGE_DYNAMICS_V_H_HPP

#include <cmath>

namespace g2o::oc {

EdgeDynamics_v_h::EdgeDynamics_v_h(int k, std::shared_ptr<OCParameters> param = nullptr)
    : _k(k), _param(param){}

void EdgeDynamics_v_h::computeEq() {
    const Vertex_v_h* v_v_k = static_cast<const Vertex_v_h*>(this->_vertices[0]);
    const Vertex_v_h* v_v_kp1 = static_cast<const Vertex_v_h*>(this->_vertices[1]);
    const Vertex_f* v_f = static_cast<const Vertex_f*>(this->_vertices[2]);
 
    double v_h_k = v_v_k->estimate();
    double v_h_kp1 = v_v_kp1->estimate();
    double f = v_f->estimate();
     

    double delta_t = _param->get_delta_t();
    double m_eq = _param->get_m_eq();
    double g = _param->get_g();
    double alpha_k = _param->get_alpha_prediction(_k);;
    double p1 = _param->get_p1();
    double p2 = _param->get_p2();

    if (_param->isLinearDynamics()) {
        double w_1 = -delta_t / m_eq *
                     (_param->get_m_total() * g * sin(alpha_k * M_PI / 180) +
                      _param->get_c_r() * _param->get_m_total() * g * cos(alpha_k * M_PI / 180) + p1);
        
        _eq[0] = v_h_kp1 - ((1 - delta_t / m_eq * p2) * v_h_k + delta_t / m_eq * (f) + w_1);
    } else {
        double f_roll =   _param->get_c_r() *  _param->get_m_total() * g * cos(alpha_k * M_PI / 180) ;
        _eq[0] = v_h_kp1 - ( v_h_k + (delta_t / m_eq) * (f - _param->get_m_total() * g * sin(alpha_k * M_PI / 180) 
                             - f_roll - _param->get_c_a() * v_h_k * v_h_k));

    }





    
}

void EdgeDynamics_v_h::linearizeOplus() {
    // The derivative of _eq[0]  with respect to the vertex[0] (v_h_k)  is - delta_t / m_eq * (1 - delta_t / m_eq * p2)
    // The derivative of _eq[0]  with respect to the vertex[1] (v_h_kp1)  is 1
    // The derivative of _eq[0]  with respect to the vertex[2] (f_t_k)  is  delta_t / m_eq
    // The derivative of _eq[0]  with respect to the vertex[3] (f_b_k)  is -delta_t / m_eq
    // Set the Jacobian's value for the scalar relationship
    // we need here carefully to define the Jacobian matrix because the solver might add   
    // Langrangian varibales so to be general we define it in the equatily edges as following
    // the dimension of this
    // initialize the
    // Jacobians to zero

    initializeJacobians();
    // The Jacobian for the Lagrangian vertex
    double delta_t = _param->get_delta_t();
    double m_eq = _param->get_m_eq();
    double p2 = _param->get_p2();
    
       const int D = _eq.size();
    auto&& J_v0 = std::get<0>(_jacobianOplus).topRows(D);
    auto&& J_v1 = std::get<1>(_jacobianOplus).topRows(D);
    auto&& J_v2 = std::get<2>(this->_jacobianOplus).topRows(D);
 
   
    
    J_v0 << - (1 - delta_t / m_eq * p2);
    J_v1 << 1;
    J_v2 << -delta_t / m_eq;
 }

bool EdgeDynamics_v_h::write(std::ostream& os) const {
    os << "EdgeDynamics_v_h, ";
    for (int i = 0; i < 1; i++) {
        os << _eq[i] << " ";
    }

    return os.good();
}

bool EdgeDynamics_v_h::read(std::istream& is) {
    for (int i = 0; i < 1; i++) {
        is >> _eq[i];
    }

    return is.good();
}


}  // namespace g2o::oc

#endif  // EDGE_DYNAMICS_V_H_HPP