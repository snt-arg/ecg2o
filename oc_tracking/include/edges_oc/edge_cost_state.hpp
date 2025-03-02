#ifndef EDGE_COST_STATE_HPP
#define EDGE_COST_STATE_HPP

namespace g2o::oc {

EdgeCost_state::EdgeCost_state(int k, std::shared_ptr<OCParameters> param = nullptr)
    : _k(k), _param(param){}

void EdgeCost_state::computeError() {
    const Vertex_v_h* v_v_h_k = static_cast<const Vertex_v_h*>(_vertices[0]);
 
    double v_h_k = v_v_h_k->estimate();   

    auto v_ref_k = _param->get_v_h_ref()[_k];

   
    _error[0]= v_h_k - v_ref_k;  
}

void EdgeCost_state::linearizeOplus() {
    // The derivative of _error[0] (v_h_k)  with respect to the vertex[0] (v_h_k)  is 1
     // Set the Jacobian's value for the scalar relationship
     
     
    auto& J_v0 = std::get<0>(this->_jacobianOplus);  
    J_v0 << 1.0;   // J(e_0,v_0[0])

}


bool EdgeCost_state::write(std::ostream& os) const {
    os << "EdgeCost_brake_change, ";
    for (int i = 0; i < 1; i++) {
        os << _error[i] << " ";
    }

    return os.good();

}

bool EdgeCost_state::read(std::istream& is) {
    for (int i = 0; i < 1; i++) {
        is >> _error[i];
    }

    return is.good();
    
}
}  // namespace g2o::oc

#endif