#ifndef EDGE_COST_INPUT_H
#define EDGE_COST_INPUT_H


#include "oc_parameters.h"
#include "oc_vertices.h"
#include "g2o/core/base_fixed_sized_edge.h"


namespace g2o::oc {

class EdgeCost_input : public g2o::BaseFixedSizedEdge<1, double, Vertex_f> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    EdgeCost_input(int k, std::shared_ptr<OCParameters> param);
    void computeError() override;
    void linearizeOplus() override; //optional
    
    bool write(std::ostream& os) const override;
    bool read(std::istream& is) override;

 

private:
    int _k;
    std::shared_ptr<OCParameters> _param;
   
};

}  // namespace g2o::oc

#include "edges_oc/edge_cost_input.hpp" 

#endif  // EDGE_COST_INPUT_H
