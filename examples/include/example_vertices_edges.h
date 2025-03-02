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

#ifndef EXAMPLE_EDGES_H
#define EXAMPLE_EDGES_H

#include <g2o/core/base_fixed_sized_edge.h>
#include <g2o/core/base_unary_edge.h>
#include <g2o/core/base_binary_edge.h>
#include <g2o/core/base_vertex.h>
#include "ecg2o/base_fixed_sized_edge_eq.h"


class VertexXY : public g2o::BaseVertex<2, Eigen::Vector2d> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
        void setToOriginImpl() override {
            _estimate << -40,-20;
        }
    
        void oplusImpl(const double* update) override {
            _estimate += Eigen::Vector2d(update);
         //   std::cout << "------------------oplusImpl: -------- " << _estimate << std::endl;
        }
    
        virtual bool read(std::istream& in) override {
            in >> _estimate[0] >> _estimate[1];
            return true;
        }
    
        virtual bool write(std::ostream& out) const override {
            out << _estimate[0] << " " << _estimate[1];
            return true;
        }
    };
    
    class VertexZ : public g2o::BaseVertex<1, double> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
        VertexZ() {}
    
        void setToOriginImpl() override {
            _estimate = -40.0;
        }
    
        void oplusImpl(const double* update) override {
           // std::cout << "------------------oplusImpl: -------- " << _estimate << std::endl;
            _estimate += update[0];
            }
    
        bool read(std::istream& /*is*/) override {
            return false;
        }
    
        bool write(std::ostream& /*os*/) const override {
            return false;
        }
    };
    
    
    // Define the EdgeZ which represents the cost function (x-2)^2 + (y-9)^2
    class  EdgeXY : public g2o::BaseUnaryEdge<2, Eigen::Vector2d, VertexXY> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
         EdgeXY() {}
    
        void computeError() override {
            const VertexXY* v = static_cast<const VertexXY*>(_vertices[0]);
            const Eigen::Vector2d& xy = v->estimate();
            _error= (xy - _measurement);
        }
    
        virtual bool read(std::istream& in) override {
            in >> _measurement[0] >> _measurement[1];
            return true;
        }
    
        virtual bool write(std::ostream& out) const override {
            out << _measurement[0] << " " << _measurement[1];
            return true;
        }
    };
    
    
    
    // Unary EdgeZ representing a single constraint
    class EdgeZ : public g2o::BaseUnaryEdge<1, double, VertexZ> {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
        EdgeZ() {}
    
        void computeError() override {
            const VertexZ* v = static_cast<const VertexZ*>(_vertices[0]);
            _error[0] = v->estimate() - _measurement;
        }
    
        bool read(std::istream& /*is*/) override {
            return false;
        }
    
        bool write(std::ostream& /*os*/) const override {
            return false;
        }
    
        void linearizeOplus() override {
            _jacobianOplusXi(0, 0) = 1;
        }
    
    
    };
    
    
    // Binary EdgeZ representing a constraint between two variables
    class  EdgeEq : public g2o::BaseFixedSizedEdgeEq<2,double, VertexXY, VertexZ > {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW    
        EdgeEq(){};        
        void computeEq() override {
            const VertexXY* vxy = static_cast<const VertexXY*>(_vertices[0]);
            const VertexZ* vz = static_cast<const VertexZ*>(_vertices[1]);
            
                auto xy = vxy->estimate();
            auto z = vz->estimate();         
            _eq[0] = (xy[1] + z - 3);
            _eq[1] = (xy[0] +xy[1] -2.5);
        } 
    
        void linearizeOplus() override {
                initializeJacobians(); // to handle bottomRows (related to the Lagrange multiplier)
                const int D = _eq.size();
                auto&& J_v0 = std::get<0>(_jacobianOplus).topRows(D); // vertex xy , v_0[0] = xy[0], v_0[1] = xy[1]
                auto&& J_v1 = std::get<1>(_jacobianOplus).topRows(D); // vertex z = v_1
                J_v0 << 0, 1,     //  J(e_0,v_0[0])  J(e_0,v_0[1])
                        1, 1;    //  J(e_1,v_0[0])  J(e_1,v_0[1])                
                J_v1 << 1,       //  J(e_0,v_1)
                        0;       //  J(e_1,v_1)
            
            
         }         
    };

    
#endif