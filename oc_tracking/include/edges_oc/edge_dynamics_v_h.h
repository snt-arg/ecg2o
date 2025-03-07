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
