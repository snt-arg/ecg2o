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


#ifndef OC_PARAMETERS_H
#define OC_PARAMETERS_H

#include <vector>
#include <data_base.h>
namespace g2o::oc {

class OCParameters {
private:
    // Parameters as const private members
    const bool linearDynamics = false ;
    const int N = 2;
    const int N_max = 10;
    const double Ts = 0.1;
    const double c_a = 0.421022;
    const double c_r = 0.00708572690577572;
    const double e_i = 0.07632073023674692;
    const double g = 9.81;
    const double m_eq = 1639.0408163265306;
    const double m_p = 200.0;
    const double m_total = 1537.0;
    const double m_v = 1337.0;
    const double p1 = 0.0;
    const double p2 = 20.0;
    const double x1 = 1.1;
    const double x2 = 2.0;
    const double Omega1 = 0.01;
    const double Omega6 = 1000.0;
    const double BrakeChangePenalty = 0.1;
    const double BrakePenalty = 0.08;
    const double TractionChangePenalty = 0.01;
    const double Brems_a1 = -47.01763813449595;
    const double Brems_b1 = 3068.0;
    const double Regen_max = 4000.0;
    const double a1 = 0.0;
    const double a2 = -206.46604377690898;
    const double b1 = 0.0;
    const double b2 = 9669.0;
    const double b3 = 0.0;
    const double b4 = 6903.0;
    const double b_max = 5000.0;
    const double d_min = 4.5;
    //  const double delta_t = 1.0;
    const double h_min = 0.75;
    const double h_soft = 1.25;
    const double delta_f_t_max = 0;
    const double delta_f_b_max = 0;
    const std::vector<double> eff_map = {
        14510.811959872348, -8.984433145528733, 0.10270476623078269,
        -1521.2579176564311, 0.8265216693384867, 61.289232062679915};
        
        const double TracForce_Prev = 111;
        const double BrakeForce_Prev = 222;
        //const std::vector<double> v_h_ref  = data_v_h_ref;
        const double P = 1000.0;
        const double Q = 1000.0;
        const double R = .0007;
        const double delta_t = 1;
        const std::vector<double> v_h_ref  = data_v_h_104_1;
        
        const std::vector<double> alpha_prediction = data_alpha_prediction;
        
        const std::vector<double> v_min_prediction = data_v_min_prediction;
        
        const std::vector<double> v_max_prediction = data_v_max_prediction;    
        
        const double initial_f = 0;
        const double initial_f_t = 99;                                        
        const double initial_f_b = 199;
        const double initial_slack = 1000;
        const double initial_v_h = 0.0;
        const double initial_d_h = 20.0;
        const double initial_d_h_0 = 10.2;
        const double initial_v_h_0 = 0.001;
        
        
        
        
        public:
        // Getters for read-only access
        bool isLinearDynamics() const { return linearDynamics; }
        int get_N() const { return N; }
        double get_Ts() const { return Ts; }
        double get_c_a() const { return c_a; }
        double get_c_r() const { return c_r; }
        double get_e_i() const { return e_i; }
        double get_g() const { return g; }
        double get_m_eq() const { return m_eq; }
        double get_m_p() const { return m_p; }
        double get_m_total() const { return m_total; }
        double get_m_v() const { return m_v; }
        double get_p1() const { return p1; }
        double get_p2() const { return p2; }
        double get_x1() const { return x1; }
        double get_x2() const { return x2; }
        double get_BrakeChangePenalty() const { return BrakeChangePenalty; }
        double get_BrakePenalty() const { return BrakePenalty; }
        double get_Brems_a1() const { return Brems_a1; }
        double get_Brems_b1() const { return Brems_b1; }
        int get_N_max() const { return N_max; }
        double get_P() const { return P; }
        double get_Q() const { return Q; }
        double get_R() const { return R; }
        double get_Omega1() const { return Omega1; }
    double get_Omega6() const { return Omega6; }
    double get_Regen_max() const { return Regen_max; }
    double get_TractionChangePenalty() const { return TractionChangePenalty; }
    double get_a1() const { return a1; }
    double get_a2() const { return a2; }
    double get_b1() const { return b1; }
    double get_b2() const { return b2; }
    double get_b3() const { return b3; }
    double get_b4() const { return b4; }
    double get_b_max() const { return b_max; }
    double get_d_min() const { return d_min; }
    double get_delta_t() const { return delta_t; }
    double get_h_min() const { return h_min; }
    double get_h_soft() const { return h_soft; }
    double get_delta_f_t_max() const { return delta_f_t_max; }
    double get_delta_f_b_max() const { return delta_f_b_max; } 
    std::vector<double> get_eff_map() const { return eff_map; }
    double get_TracForce_Prev() const { return TracForce_Prev; }
    double get_BrakeForce_Prev() const { return BrakeForce_Prev; } 
    double get_v_h_ref(int k) const { return v_h_ref.at(k); }
    double get_alpha_prediction(int k) const { return alpha_prediction.at(k); }
    double get_v_min_prediction(int k) const { return v_min_prediction.at(k); }
    double get_v_max_prediction(int k) const { return v_max_prediction.at(k); }
    std::vector<double> get_v_h_ref() const { return v_h_ref; }
    std::vector<double> get_alpha_prediction() const { return alpha_prediction; }
    std::vector<double> get_v_min_prediction() const { return v_min_prediction; }
    std::vector<double> get_v_max_prediction() const { return v_max_prediction; }
    double get_initial_f() const { return initial_f; }
    double get_initial_f_t() const { return initial_f_t; }
    double get_initial_f_b() const { return initial_f_b; }
    double get_initial_slack() const { return initial_slack; }
    double get_initial_v_h() const { return initial_v_h; }
    double get_initial_d_h() const { return initial_d_h; }
    double get_initial_d_h_0() const { return initial_d_h_0; }
    double get_initial_v_h_0() const { return initial_v_h_0; }



    // Setters for controlled modification (using const_cast)
    void set_linearDynamics(bool value) { const_cast<bool&>(linearDynamics) = value; }  
    void set_N(int value) { const_cast<int&>(N) = value; }
    void set_Ts(double value) { const_cast<double&>(Ts) = value; }
    void set_c_a(double value) { const_cast<double&>(c_a) = value; }
    void set_c_r(double value) { const_cast<double&>(c_r) = value; }
    void set_e_i(double value) { const_cast<double&>(e_i) = value; }
    void set_m_eq(double value) { const_cast<double&>(m_eq) = value; }
    void set_m_p(double value) { const_cast<double&>(m_p) = value; }
    void set_m_total(double value) { const_cast<double&>(m_total) = value; }
    void set_m_v(double value) { const_cast<double&>(m_v) = value; }
    void set_p1(double value) { const_cast<double&>(p1) = value; }
    void set_p2(double value) { const_cast<double&>(p2) = value; }
    void set_x1(double value) { const_cast<double&>(x1) = value; }
    void set_x2(double value) { const_cast<double&>(x2) = value; }
    void set_BrakeChangePenalty(double value) { const_cast<double&>(BrakeChangePenalty) = value; }
    void set_BrakePenalty(double value) { const_cast<double&>(BrakePenalty) = value; }
    void set_Brems_a1(double value) { const_cast<double&>(Brems_a1) = value; }
    void set_Brems_b1(double value) { const_cast<double&>(Brems_b1) = value; }
    void set_P(double value) { const_cast<double&>(P) = value; }
    void set_Q(double value) { const_cast<double&>(Q) = value; }
    void set_R(double value) { const_cast<double&>(R) = value; }
    void set_Omega1(double value) { const_cast<double&>(Omega1) = value; }
    void set_Omega6(double value) { const_cast<double&>(Omega6) = value; }
    void set_Regen_max(double value) { const_cast<double&>(Regen_max) = value; }
    void set_TractionChangePenalty(double value) { const_cast<double&>(TractionChangePenalty) = value; }
    void set_a1(double value) { const_cast<double&>(a1) = value; }
    void set_a2(double value) { const_cast<double&>(a2) = value; }
    void set_b1(double value) { const_cast<double&>(b1) = value; }
    void set_b2(double value) { const_cast<double&>(b2) = value; }
    void set_b3(double value) { const_cast<double&>(b3) = value; }
    void set_b4(double value) { const_cast<double&>(b4) = value; }
    void set_b_max(double value) { const_cast<double&>(b_max) = value; }
    void set_d_min(double value) { const_cast<double&>(d_min) = value; }
    void set_delta_t(double value) { const_cast<double&>(delta_t) = value; }
    void set_h_min(double value) { const_cast<double&>(h_min) = value; }
    void set_h_soft(double value) { const_cast<double&>(h_soft) = value; }
    void set_delta_f_t_max(double value) { const_cast<double&>(delta_f_t_max) = value; }
    void set_delta_f_b_max(double value) { const_cast<double&>(delta_f_b_max) = value; }
    void set_eff_map(std::vector<double> value) { const_cast<std::vector<double>&>(eff_map) = value; }
    void set_TracForce_Prev(double value) { const_cast<double&>(TracForce_Prev) = value; }
    void set_BrakeForce_Prev(double value) { const_cast<double&>(BrakeForce_Prev) = value; }
    void set_v_h_ref(std::vector<double> value) { const_cast<std::vector<double>&>(v_h_ref) = value; }
    void set_alpha_prediction(std::vector<double> value) { const_cast<std::vector<double>&>(alpha_prediction) = value; }
    void set_v_min_prediction(std::vector<double> value) { const_cast<std::vector<double>&>(v_min_prediction) = value; }
    void set_v_max_prediction(std::vector<double> value) { const_cast<std::vector<double>&>(v_max_prediction) = value; }
    void set_initial_f(double value) { const_cast<double&>(initial_f) = value; }
    void set_initial_f_t(double value) { const_cast<double&>(initial_f_t) = value; }
    void set_initial_f_b(double value) { const_cast<double&>(initial_f_b) = value; }
    void set_initial_slack(double value) { const_cast<double&>(initial_slack) = value; }
    void set_initial_v_h(double value) { const_cast<double&>(initial_v_h) = value; }
    void set_initial_d_h(double value) { const_cast<double&>(initial_d_h) = value; }
    void set_initial_d_h_0(double value) { const_cast<double&>(initial_d_h_0) = value; }
    void set_initial_v_h_0(double value) { const_cast<double&>(initial_v_h_0) = value; }
    

};

}  // namespace g2o::oc

#endif  // OC_PARAMETERS_H
