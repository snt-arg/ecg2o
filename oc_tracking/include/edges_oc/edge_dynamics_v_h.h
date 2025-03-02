#ifndef EDGE_DYNAMICS_V_H_H
#define EDGE_DYNAMICS_V_H_H


#include "oc_parameters.h"
#include "oc_vertices.h"
#include "ecg2o/base_fixed_sized_edge_eq.h"

namespace g2o::oc {


class EdgeDynamics_v_h : public g2o::BaseFixedSizedEdgeEq<1, double, Vertex_v_h, Vertex_v_h, Vertex_f> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    EdgeDynamics_v_h(int k, std::shared_ptr<OCParameters> param);
    void computeEq() override;
    void linearizeOplus() override ; //optional

    bool write(std::ostream& os) const override;
    bool read(std::istream& is) override;


private:
    int _k;
     std::shared_ptr<OCParameters> _param;

};

}  // namespace g2o::oc
#include "edges_oc/edge_dynamics_v_h.hpp"
#endif  // EDGE_DYNAMICS_V_H_H
