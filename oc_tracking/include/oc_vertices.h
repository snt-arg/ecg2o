#ifndef  G2O_GRAPH_VERTEXES_OC_H_
#define  G2O_GRAPH_VERTEXES_OC_H_

#include <g2o/core/base_vertex.h>
#include <Eigen/Core>

namespace g2o::oc {


// Traits for defining the Estimate type
template <int D>
struct EstimateTraits {
    using Type = Eigen::Matrix<double, D, 1>;
};

// Specialization for D = 1
template <>
struct EstimateTraits<1> {
    using Type = double;
};

// VertexOC class template declaration
template <int D>
class VertexOC : public g2o::BaseVertex<D, typename EstimateTraits<D>::Type> {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    using Base = g2o::BaseVertex<D, typename EstimateTraits<D>::Type>;
    using EstimateType = typename EstimateTraits<D>::Type;

    // Set the vertex state to the origin
    virtual void setToOriginImpl() override;

    // Apply an increment (delta) to the state
    virtual void oplusImpl(const double* update) override;

    // Write the vertex state to an output stream
    virtual bool write(std::ostream& os) const override;

    // Read the vertex state from an input stream
    virtual bool read(std::istream& is) override;
};

}  // namespace g2o::oc

# include "oc_vertices.hpp"


// List the need vertices 
class Vertex_scalar : public g2o::oc::VertexOC<1> { };   // 1D vertex for longitudinal velocity

class Vertex_v_h : public g2o::oc::VertexOC<1> { };   // 1D vertex for longitudinal velocity
class Vertex_f : public g2o::oc::VertexOC<1> { };   // 1D vertex for traction force
 
 
#endif  // G2O_GRAPH_VERTEXES_OC_H_

