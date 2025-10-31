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

#include <trackcpp/auxiliary.h>
#include <string>
#include <fstream>
#include <cmath>


template <typename T>
T ay(const double& brho, const double& kx, const double& ks, const std::vector<std::vector<double>>& coefs, const T& x, const T& y, const double& s) 
{
    T ay_ = 0.0;
    int M = coefs.size();
    int N = coefs[0].size();

    for (int m = 1; m <= M; ++m) {
        for (int n = 1; n <= N; ++n) {
            double ky = std::sqrt(std::pow(m * kx, 2) + std::pow(n * ks, 2));
            double fac = coefs[m - 1][n - 1] * (m * kx) / (n * ks * ky);
            ay_ += fac * sin(m * kx * x) * sinh(ky * y) * std::cos(n * ks * s);
        }
    }

    return -1 * ay_/brho;
}

template <typename T>
T ax(const double& brho, const double& kx, const double& ks, const std::vector<std::vector<double>>& coefs, const T& x, const T& y, const double& s)
{
    T ax_ = 0.0;
    const int M = coefs.size();
    const int N = coefs[0].size();

    for (int m = 1; m <= M; ++m) {
        for (int n = 1; n <= N; ++n) {
            double ky = std::sqrt(std::pow(m * kx, 2) + std::pow(n * ks, 2));
            double fac = coefs[m - 1][n - 1] / (n * ks);
            ax_ += fac * cos(m * kx * x) * cosh(ky * y) * std::cos(n * ks * s);
        }
    }

    return -1 * ax_/brho;
}

template <typename T>
T inty_day_dx(const double& brho, const double& kx, const double& ks, const std::vector<std::vector<double>>& coefs, const T& x, const T& y, const double& s)
{
    T day_dx = 0.0;
    int M = coefs.size();
    int N = coefs[0].size();

    for (int m = 1; m <= M; ++m) {
        for (int n = 1; n <= N; ++n) {
            double ky = std::sqrt(std::pow(m * kx, 2) + std::pow(n * ks, 2));
            double fac = coefs[m - 1][n - 1] * std::pow(m * kx, 2) / (n * ks * std::pow(ky, 2));
            day_dx += fac * cos(m * kx * x)* std::cos(n * ks * s)* (cosh(ky * y) - 1.0);
        }
    }

    return -1 * day_dx/brho;
}

template <typename T>
T intx_dax_dy(const double& brho, const double& kx, const double& ks, const std::vector<std::vector<double>>& coefs, const T& x, const T& y, const double& s)
{
    T dax_dy = 0.0;
    int M = coefs.size();
    int N = coefs[0].size();

    for (int m = 1; m <= M; ++m) {
        for (int n = 1; n <= N; ++n) {
            double ky = std::sqrt(std::pow(m * kx, 2) + std::pow(n * ks, 2));
            double fac = coefs[m - 1][n - 1] * ky / (n * ks * m * kx);
            dax_dy += fac * sin(m * kx * x) * std::cos(n * ks * s) * sinh(ky * y);
        }
    }

    return -1* dax_dy/brho;
}

template <typename T>
T calc_D(const double& beta0, const T& delta) {
    return sqrt(1.0 + 2.0 * delta / beta0 + delta * delta);
}

template <typename T>
void exp_h1_z(const double& beta0, Pos<T>& map, double step) {
    T d = calc_D(beta0, map.de);
    T factor = (1.0 / beta0 - (1.0 / beta0 + map.de) / d);
    map.dl += factor * step / 2.0;
}

template <typename T>
void exp_h1_s(T& s, T step) {
    s += step / 2.0;
}

template <typename T>
void exp_iy_px(const double& brho, const double& kx, const double& ks, const std::vector<std::vector<double>>& coefs, Pos<T>& map, double s, int sign, double step) {
    T factor = inty_day_dx(brho, kx, ks, coefs, map.rx, map.ry, s) * sign * -1;
    map.px += factor;
}


template <typename T>
void exp_iy_py(const double& brho, const double& kx, const double& ks, const std::vector<std::vector<double>>& coefs, Pos<T>& map, double s, int sign, double step) {
    T factor = ay(brho, kx, ks, coefs, map.rx, map.ry, s) * sign * -1;
    map.py += factor;
}

template <typename T>
void exp_h2_y(const double& beta0, Pos<T>& map, double step) {
    T d = calc_D(beta0, map.de);
    map.ry += map.py/d*step/2.0;
}


template <typename T>
void exp_h2_z(const double& beta0, Pos<T>& map, double step) {
    T d = calc_D(beta0, map.de);
    T factor = (pow(map.py, 2)*(1/beta0+map.de)/(2*pow(d, 3)));
    map.dl -= factor*step/2;
}


template <typename T>
void exp_ix_px(const double& brho, const double& kx, const double& ks, const std::vector<std::vector<double>>& coefs, Pos<T>& map, double s, int sign, double step) {
    T factor = ax(brho, kx, ks, coefs, map.rx, map.ry, s) * sign * -1;
    map.px += factor;
}


template <typename T>
void exp_ix_py(const double& brho, const double& kx, const double& ks, const std::vector<std::vector<double>>& coefs, Pos<T>& map, double s, int sign, double step) {
    T factor = intx_dax_dy(brho, kx, ks, coefs, map.rx, map.ry, s) * sign * -1;
    map.py += factor;
}


template <typename T>
void exp_h3_x(const double& beta0, Pos<T>& map, double step) {
    T d = calc_D(beta0, map.de);
    map.rx += map.px/d*step;
}

template <typename T>
void exp_h3_z(const double& beta0, Pos<T>& map, double step) {
    T d = calc_D(beta0, map.de);
    T factor = (pow(map.px, 2)*(1/beta0+map.de)/(2*pow(d, 3)));
    map.dl -= factor*step;
}

template <typename T>
void prop_h1(const double& beta0, Pos<T>& map, double& s, double step) {
    exp_h1_z(beta0, map, step);
    exp_h1_s(s, step);
}

template <typename T>
void prop_h2(const double& beta0, Pos<T>& map, double step) {
    exp_h2_y(beta0, map, step);
    exp_h2_z(beta0, map, step);
}

template <typename T>
void prop_h3(const double& beta0, Pos<T>& map, double step) {
    exp_h3_x(beta0, map, step);
    exp_h3_z(beta0, map, step);
}

template <typename T>
void prop_ix(const double& brho, const double& kx, const double& ks, const std::vector<std::vector<double>>& coefs, Pos<T>& map, double s, int sign, double step) {
    exp_ix_px(brho, kx, ks, coefs, map, s, sign, step);
    exp_ix_py(brho, kx, ks, coefs, map, s, sign, step);

}

template <typename T>
void prop_iy(const double& brho, const double& kx, const double& ks, const std::vector<std::vector<double>>& coefs, Pos<T>& map, double s, int sign, double step) {
    exp_iy_px(brho, kx, ks, coefs, map, s, sign, step);
    exp_iy_py(brho, kx, ks, coefs, map, s, sign, step);

}

template <typename T>
void prop_step(const double& beta0, const double& brho, const double& kx, const double& ks, const std::vector<std::vector<double>>& coefs, Pos<T>& map, double& s, double step) {
    prop_h1(beta0, map, s, step);
    prop_iy(brho, kx, ks, coefs, map, s, +1, step);
    prop_h2(beta0, map, step);
    prop_iy(brho, kx, ks, coefs, map, s, -1, step);
    prop_ix(brho, kx, ks, coefs, map, s, +1, step);
    prop_h3(beta0, map, step);
    prop_ix(brho, kx, ks, coefs, map, s, -1, step);
    prop_iy(brho, kx, ks, coefs, map, s, +1, step);
    prop_h2(beta0, map, step);
    prop_iy(brho, kx, ks, coefs, map, s, -1, step);
    prop_h1(beta0, map, s, step);
}
