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

