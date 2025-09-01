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

#ifndef G2O_BASE_FIXED_SIZED_EDGE_EQ_GN
#define G2O_BASE_FIXED_SIZED_EDGE_EQ_GN

#include "g2o/core/base_fixed_sized_edge.h"
#include "vertex_lagrange_multiplier.h"
#include <Eigen/Core>
#include <vector>

namespace g2o {

template <int D, typename E, typename... VertexTypes>
class G2O_CORE_API BaseFixedSizedEdgeEq
    : public BaseFixedSizedEdge<2 * D, E, VertexTypes...,
                                VertexLagrangeMultiplier<D>> {
public:
  using VoidEdgeFuncType =
      std::function<void(BaseFixedSizedEdgeEq<D, E, VertexTypes...> &)>;
  using InformationType =
      typename BaseFixedSizedEdge<D, E, VertexTypes...>::InformationType;

  BaseFixedSizedEdgeEq();

  // User-defined equality computation
  virtual void computeEq() = 0;

  void computeError() override;
  void constructQuadraticForm() override;

  void setVertexLagrangeMultiplierId(int id);
  int getVertexLagrangeMultiplierId() const;

  // Setter for the function pointers
  void setConstructQuadraticFormImpl(const VoidEdgeFuncType &func);

  // Setter
  void
  setLagrangeMultiplier(const Eigen::Matrix<double, D, 1> &lagrangeMuliplier);
  void setRho(const Eigen::Matrix<double, D, 1> &rho);
  void setConstraintViolationPrev(
      Eigen::Matrix<double, D, 1> &&input,
      size_t start = 0); // in AL method we use the information matrix to save
                         // the constraint violation and rho_bar

  // Getter functions
  Eigen::Matrix<double, D, 1> lagrangeMultiplier() const;
  Eigen::Matrix<double, D, 1> rho() const;

  bool read(std::istream &is) override;
  bool write(std::ostream &os) const override;
  void initializeInformationMatrix();

protected:
  void setInformation(const InformationType &information);
  Eigen::Matrix<double, D, 1> _eq = Eigen::Matrix<double, D, 1>::Zero();
  int _VertexLagrangeMultiplierId = -1;
  void initializeJacobians();

  VoidEdgeFuncType _constructQuadraticFormImpl = nullptr;

  // In case Rho is used
  Eigen::Matrix<double, D, 1> _rho;
};

} // namespace g2o

namespace g2o {

template <int D, typename E, typename... VertexTypes>
BaseFixedSizedEdgeEq<D, E, VertexTypes...>::BaseFixedSizedEdgeEq() {
  this->_information.setZero();
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeEq<D, E, VertexTypes...>::initializeJacobians() {
  // Helper lambda to set each Jacobian to zero
  auto setZeroHelper = [](auto &jacobian) { jacobian.setZero(); };

  // Apply the lambda to each Jacobian in the tuple
  std::apply(
      [&setZeroHelper](auto &...jacobians) {
        (setZeroHelper(jacobians),
         ...); // Expand and apply to all tuple elements
      },
      this->_jacobianOplus);

  // Set the Jacobian for the Lagrangian vertex
  if constexpr (sizeof...(VertexTypes) > 0) {
    auto &J_vLagrangian =
        std::get<sizeof...(VertexTypes)>(this->_jacobianOplus);
    J_vLagrangian.setZero();
    J_vLagrangian.template block<D, D>(D, 0).setIdentity();
  }
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeEq<D, E, VertexTypes...>::computeError() {
  if (this->_vertices.empty() ||
      this->_vertices.size() < sizeof...(VertexTypes) + 1) {
    throw std::runtime_error("Lagrangian vertex missing or invalid index");
  }

  const auto *vLagrangian =
      static_cast<const VertexLagrangeMultiplier<D> *>(this->_vertices.back());
  Eigen::Matrix<double, D, 1> Lagrangian = vLagrangian->estimate();

  // Compute equality using the user-defined method
  this->computeEq();

  // Update the error vector
  this->_error.segment(0, D) = _eq;
  this->_error.segment(D, D) = Lagrangian.segment(0, D);
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeEq<D, E, VertexTypes...>::initializeInformationMatrix() {
  this->_information.setZero();
  this->_information.template block<D, D>(0, D).setIdentity();
  this->_information.template block<D, D>(D, 0).setIdentity();
}

template <int D, typename E, typename... VertexTypes>
int BaseFixedSizedEdgeEq<D, E, VertexTypes...>::getVertexLagrangeMultiplierId()
    const {
  return _VertexLagrangeMultiplierId;
}

template <int D, typename E, typename... VertexTypes>
bool BaseFixedSizedEdgeEq<D, E, VertexTypes...>::write(std::ostream &os) const {
  for (int i = 0; i < D; ++i) {
    os << this->_error[i] << " ";
  }
  return os.good();
}

template <int D, typename E, typename... VertexTypes>
bool BaseFixedSizedEdgeEq<D, E, VertexTypes...>::read(std::istream &is) {
  for (int i = 0; i < D; ++i) {
    is >> this->_error[i];
  }
  return is.good();
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeEq<D, E, VertexTypes...>::setConstructQuadraticFormImpl(
    const VoidEdgeFuncType &func) {
  _constructQuadraticFormImpl = func;
}

// Override the constructQuadraticForm method
template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeEq<D, E, VertexTypes...>::constructQuadraticForm() {
  // Ensure that a function has been assigned before calling it
  if (_constructQuadraticFormImpl) {
    _constructQuadraticFormImpl(
        *this); // Call the assigned function (or lambda)
  } else {
    // If no function is assigned, call the base class implementation (if
    // needed)
    this->BaseFixedSizedEdge<
        2 * D, E, VertexTypes...,
        VertexLagrangeMultiplier<D>>::constructQuadraticForm();
  }
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeEq<D, E, VertexTypes...>::setInformation(
    const InformationType &information) {
  this->_information = information;
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeEq<D, E, VertexTypes...>::setConstraintViolationPrev(
    Eigen::Matrix<double, D, 1> &&input, size_t start) {
  this->_information.diagonal().segment(start, D) = input;
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeEq<D, E, VertexTypes...>::setVertexLagrangeMultiplierId(
    int id) {
  _VertexLagrangeMultiplierId = id;
}

// Setter functions
template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeEq<D, E, VertexTypes...>::setLagrangeMultiplier(
    const Eigen::Matrix<double, D, 1> &lagrangeMultiplier) {
  // Ensure the last vertex is of type VertexLagrangeMultiplier<D>
  if (this->_vertices.empty() ||
      this->_vertices.size() < sizeof...(VertexTypes) + 1) {
    throw std::runtime_error("Lagrangian vertex missing or invalid index");
  }

  // Get the Lagrangian multiplier vertex
  auto *vLagrangian =
      static_cast<VertexLagrangeMultiplier<D> *>(this->_vertices.back());

  // Set the estimate of the Lagrangian multiplier vertex
  vLagrangian->setEstimate(lagrangeMultiplier);
}

template <int D, typename E, typename... VertexTypes>
void BaseFixedSizedEdgeEq<D, E, VertexTypes...>::setRho(
    const Eigen::Matrix<double, D, 1> &rho) {
  _rho = rho;
}

// Getter functions

template <int D, typename E, typename... VertexTypes>
Eigen::Matrix<double, D, 1>
BaseFixedSizedEdgeEq<D, E, VertexTypes...>::lagrangeMultiplier() const {
  // Ensure the vertex list is valid
  if (this->_vertices.empty() ||
      this->_vertices.size() < sizeof...(VertexTypes) + 1) {
    throw std::runtime_error("Lagrangian vertex missing or invalid index");
  }

  // Get the Lagrangian multiplier vertex
  const auto *vLagrangian =
      static_cast<const VertexLagrangeMultiplier<D> *>(this->_vertices.back());

  // Return the estimate of the Lagrange multiplier vertex
  return vLagrangian->estimate();
}

template <int D, typename E, typename... VertexTypes>
Eigen::Matrix<double, D, 1>
BaseFixedSizedEdgeEq<D, E, VertexTypes...>::rho() const {
  return _rho;
}

} // namespace g2o

#endif // G2O_BASE_FIXED_SIZED_EDGE_EQ_H
