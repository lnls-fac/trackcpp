// TRACKCPP - Particle tracking code
// Copyright (C) 2015  LNLS Accelerator Physics Group
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef _FIELD3D_H
#define _FIELD3D_H

#include "auxiliary.h"
#include <string>
#include <vector>


template <typename T>
T ay(const double& brho, const double& kx, const double& ks, const std::vector<std::vector<double>>& coefs, const T& x, const T& y, const double& s) 


template <typename T>
T ax(const double& brho, const double& kx, const double& ks, const std::vector<std::vector<double>>& coefs, const T& x, const T& y, const double& s)


template <typename T>
T inty_day_dx(const double& brho, const double& kx, const double& ks, const std::vector<std::vector<double>>& coefs, const T& x, const T& y, const double& s)


template <typename T>
T intx_dax_dy(const double& brho, const double& kx, const double& ks, const std::vector<std::vector<double>>& coefs, const T& x, const T& y, const double& s)


template <typename T>
T calc_D(const double& beta0, const T& delta)

template <typename T>
void exp_h1_z(const double& beta0, Pos<T>& map, double step) 


void exp_h1_s(double& s, double step)

template <typename T>
void exp_iy_px(const double& brho, const double& kx, const double& ks, const std::vector<std::vector<double>>& coefs, Pos<T>& map, double s, int sign = 1, double step)


template <typename T>
void exp_iy_py(const double& brho, const double& kx, const double& ks, const std::vector<std::vector<double>>& coefs, Pos<T>& map, double s, int sign = 1, double step)

template <typename T>
void exp_h2_y(const double& beta0, Pos<T>& map, double step)


template <typename T>
void exp_h2_z(const double& beta0, Pos<T>& map, double step)


template <typename T>
void exp_ix_px(const double& brho, const double& kx, const double& ks, const std::vector<std::vector<double>>& coefs, Pos<T>& map, double s, int sign = 1, double step)


template <typename T>
void exp_ix_py(const double& brho, const double& kx, const double& ks, const std::vector<std::vector<double>>& coefs, Pos<T>& map, double s, int sign = 1, double step)


template <typename T>
void exp_h3_x(const double& beta0, Pos<T>& map, double step)

template <typename T>
void exp_h3_z(const double& beta0, Pos<T>& map, double step)

template <typename T>
void prop_h1(const double& beta0, Pos<T>& map, double& s, double step)

template <typename T>
void prop_h2(const double& beta0, Pos<T>& map, double step)

template <typename T>
void prop_h3(const double& beta0, Pos<T>& map, double step)

template <typename T>
void prop_ix(const double& brho, const double& kx, const double& ks, const std::vector<std::vector<double>>& coefs, Pos<T>& map, double s, int sign = 1, double step)

template <typename T>
void prop_iy(const double& brho, const double& kx, const double& ks, const std::vector<std::vector<double>>& coefs, Pos<T>& map, double s, int sign = 1, double step)

template <typename T>
void prop_step(const double& beta0, const double& brho, const double& kx, const double& ks, const std::vector<std::vector<double>>& coefs, Pos<T>& map, double& s, double step)


#endif
