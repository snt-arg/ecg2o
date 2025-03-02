#ifndef EDGE_COST_STATE_H
#define EDGE_COST_STATE_H

#include "g2o/core/base_fixed_sized_edge.h"
#include "oc_parameters.h"
#include "oc_vertices.h"


namespace g2o::oc {

class EdgeCost_state : public g2o::BaseFixedSizedEdge<1, double, Vertex_v_h> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    EdgeCost_state(int k, std::shared_ptr<OCParameters> param);

    void computeError() override;
    void linearizeOplus() override; //optional
    
    bool write(std::ostream& os) const override;
    bool read(std::istream& is) override;

private:
    int _k;
    std::shared_ptr<OCParameters> _param;

};

}  // namespace g2o::oc

#include "edges_oc/edge_cost_state.hpp"

#endif  // EDGE_COST_BRAKE_CHANGE_H
