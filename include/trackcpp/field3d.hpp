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
T ay(const double& brho, const double& kx, const double& ks, 
              const std::vector<std::vector<double>>& coefs, const std::vector<std::vector<double>>& coefs2,
              const std::vector<std::vector<double>>& coefs3, const std::vector<std::vector<double>>& coefs4,
              const T& x, const T& y, const double& s) 
{
    T ay_ = 0.0;
    int M = coefs.size();
    int N = coefs[0].size();

    for (int m = 1; m <= M; ++m) {
        for (int n = 1; n <= N; ++n) {
            double ky = std::sqrt(std::pow(m * kx, 2) + std::pow(n * ks, 2));
            double fac1 = coefs[m - 1][n - 1] * (m * kx) / (n * ks * ky);
            double fac2 = coefs2[m - 1][n - 1] * (m * kx) / (n * ks * ky);
            double fac3 = -coefs3[m - 1][n - 1] * (m * kx) / (n * ks * ky);
            double fac4 = -coefs4[m - 1][n - 1] * (m * kx) / (n * ks * ky);
            ay_ += fac1 * sin(m * kx * x) * sinh(ky * y) * std::sin(n * ks * s);
            ay_ += fac2 * cos(m * kx * x) * sinh(ky * y) * std::sin(n * ks * s);
            ay_ += fac3 * sin(m * kx * x) * sinh(ky * y) * std::cos(n * ks * s);
            ay_ += fac4 * cos(m * kx * x) * sinh(ky * y) * std::cos(n * ks * s);
        }
    }

    return 1 * ay_/brho;
}

template <typename T>
T ax(const double& brho, const double& kx, const double& ks, 
              const std::vector<std::vector<double>>& coefs, const std::vector<std::vector<double>>& coefs2,
              const std::vector<std::vector<double>>& coefs3, const std::vector<std::vector<double>>& coefs4,
              const T& x, const T& y, const double& s)
{
    T ax_ = 0.0;
    const int M = coefs.size();
    const int N = coefs[0].size();

    for (int m = 1; m <= M; ++m) {
        for (int n = 1; n <= N; ++n) {
            double ky = std::sqrt(std::pow(m * kx, 2) + std::pow(n * ks, 2));
            double fac1 = coefs[m - 1][n - 1] / (n * ks);
            double fac2 = coefs2[m - 1][n - 1] / (n * ks);
            double fac3 = -coefs3[m - 1][n - 1] / (n * ks);
            double fac4 = -coefs4[m - 1][n - 1] / (n * ks);
            ax_ += fac1 * cos(m * kx * x) * cosh(ky * y) * std::sin(n * ks * s);
            ax_ += fac2 * sin(m * kx * x) * cosh(ky * y) * std::sin(n * ks * s);
            ax_ += fac3 * cos(m * kx * x) * cosh(ky * y) * std::cos(n * ks * s);
            ax_ += fac4 * sin(m * kx * x) * cosh(ky * y) * std::cos(n * ks * s);
        }
    }

    return 1 * ax_/brho;
}

template <typename T>
T inty_day_dx(const double& brho, const double& kx, const double& ks, 
              const std::vector<std::vector<double>>& coefs, const std::vector<std::vector<double>>& coefs2,
              const std::vector<std::vector<double>>& coefs3, const std::vector<std::vector<double>>& coefs4,
              const T& x, const T& y, const double& s)
{
    T day_dx = 0.0;
    int M = coefs.size();
    int N = coefs[0].size();

    for (int m = 1; m <= M; ++m) {
        for (int n = 1; n <= N; ++n) {
            double ky = std::sqrt(std::pow(m * kx, 2) + std::pow(n * ks, 2));
            double fac1 = coefs[m - 1][n - 1] * std::pow(m * kx, 2) / (n * ks * std::pow(ky, 2));
            double fac2 = -coefs2[m - 1][n - 1] * std::pow(m * kx, 2) / (n * ks * std::pow(ky, 2));
            double fac3 = -coefs3[m - 1][n - 1] * std::pow(m * kx, 2) / (n * ks * std::pow(ky, 2));
            double fac4 = coefs4[m - 1][n - 1] * std::pow(m * kx, 2) / (n * ks * std::pow(ky, 2));
            day_dx += fac1 * cos(m * kx * x) * std::sin(n * ks * s) * (cosh(ky * y) - 1.0);
            day_dx += fac2 * sin(m * kx * x) * std::sin(n * ks * s) * (cosh(ky * y) - 1.0);
            day_dx += fac3 * cos(m * kx * x) * std::cos(n * ks * s) * (cosh(ky * y) - 1.0);
            day_dx += fac4 * sin(m * kx * x) * std::cos(n * ks * s) * (cosh(ky * y) - 1.0);
        }
    }

    return 1 * day_dx/brho;
}

template <typename T>
T intx_dax_dy(const double& brho, const double& kx, const double& ks, 
              const std::vector<std::vector<double>>& coefs, const std::vector<std::vector<double>>& coefs2,
              const std::vector<std::vector<double>>& coefs3, const std::vector<std::vector<double>>& coefs4,
              const T& x, const T& y, const double& s)
{
    T dax_dy = 0.0;
    int M = coefs.size();
    int N = coefs[0].size();

    for (int m = 1; m <= M; ++m) {
        for (int n = 1; n <= N; ++n) {
            double ky = std::sqrt(std::pow(m * kx, 2) + std::pow(n * ks, 2));
            double fac1 = coefs[m - 1][n - 1] * ky / (n * ks * m * kx);
            double fac2 = -coefs2[m - 1][n - 1] * ky / (n * ks * m * kx);
            double fac3 = -coefs3[m - 1][n - 1] * ky / (n * ks * m * kx);
            double fac4 = coefs4[m - 1][n - 1] * ky / (n * ks * m * kx);
            dax_dy += fac1 * sin(m * kx * x) * std::sin(n * ks * s) * sinh(ky * y);
            dax_dy += fac2 * cos(m * kx * x) * std::sin(n * ks * s) * sinh(ky * y);
            dax_dy += fac3 * sin(m * kx * x) * std::cos(n * ks * s) * sinh(ky * y);
            dax_dy += fac4 * cos(m * kx * x) * std::cos(n * ks * s) * sinh(ky * y);
        }
    }

    return 1* dax_dy/brho;
}

template <typename T>
T calc_D(const double& beta0, const T& delta) {
    // return sqrt(1.0 + 2.0 * delta / beta0 + delta * delta);
    return 1.0 + delta;
}

template <typename T>
void exp_h1_z(const double& beta0, Pos<T>& map, const T& d, double step) {
    T factor = (1.0 / beta0 - (1.0 / beta0 + map.de) / d);
    map.dl += factor * step / 2.0;
}

template <typename T>
void exp_h1_s(T& s, T step) {
    s += step / 2.0;
}

template <typename T>
void exp_iy_px(const double& brho, const double& kx, const double& ks, 
               const std::vector<std::vector<double>>& coefs, const std::vector<std::vector<double>>& coefs2,
               const std::vector<std::vector<double>>& coefs3, const std::vector<std::vector<double>>& coefs4,
               Pos<T>& map, double s, int sign, double step) {
    T factor = inty_day_dx(brho, kx, ks, coefs, coefs2, coefs3, coefs4, map.rx, map.ry, s) * sign * -1.0;
    map.px += factor;
}


template <typename T>
void exp_iy_py(const double& brho, const double& kx, const double& ks, 
               const std::vector<std::vector<double>>& coefs, const std::vector<std::vector<double>>& coefs2,
               const std::vector<std::vector<double>>& coefs3, const std::vector<std::vector<double>>& coefs4,
               Pos<T>& map, double s, int sign, double step) {
    T factor = ay(brho, kx, ks, coefs, coefs2, coefs3, coefs4, map.rx, map.ry, s) * sign * -1.0;
    map.py += factor;
}

template <typename T>
void exp_h2_y(const double& beta0, Pos<T>& map, const T& pnorm, double step) {
    map.ry += 0.5*step*pnorm*map.py;
}


template <typename T>
void exp_h2_z(const double& beta0, Pos<T>& map, const T& pnorm, double step) {
    // T factor = (map.py * map.py)*(1/beta0+map.de)/(2*d*d*d);
    // T factor = (map.py * map.py)/d/d*0.5;
    map.dl += 0.25*step*pnorm*pnorm*(map.py * map.py);
}


template <typename T>
void exp_ix_px(const double& brho, const double& kx, const double& ks, 
               const std::vector<std::vector<double>>& coefs, const std::vector<std::vector<double>>& coefs2,
               const std::vector<std::vector<double>>& coefs3, const std::vector<std::vector<double>>& coefs4,
               Pos<T>& map, double s, int sign, double step) {
    T factor = ax(brho, kx, ks, coefs, coefs2, coefs3, coefs4, map.rx, map.ry, s) * sign * -1.0;
    map.px += factor;
}


template <typename T>
void exp_ix_py(const double& brho, const double& kx, const double& ks, 
               const std::vector<std::vector<double>>& coefs, const std::vector<std::vector<double>>& coefs2,
               const std::vector<std::vector<double>>& coefs3, const std::vector<std::vector<double>>& coefs4,
               Pos<T>& map, double s, int sign, double step) {
    T factor = intx_dax_dy(brho, kx, ks, coefs, coefs2, coefs3, coefs4, map.rx, map.ry, s) * sign * -1.0;
    map.py += factor;
}


template <typename T>
void exp_h3_x(const double& beta0, Pos<T>& map, const T& pnorm, double step) {
    map.rx += step*pnorm*map.px;
}

template <typename T>
void exp_h3_z(const double& beta0, Pos<T>& map, const T& pnorm, double step) {
    // T factor = (map.px*map.px)*(1/beta0+map.de)/(2*d*d*d);
    // T factor = (map.px*map.px)/d/d*0.5;
    map.dl += 0.5*step*pnorm*pnorm*(map.px * map.px);
}

template <typename T>
void prop_h1(const double& beta0, Pos<T>& map, const T& d, double& s, double step) {
    // exp_h1_z(beta0, map, d, step);
    exp_h1_s(s, step);
}

template <typename T>
void prop_h2(const double& beta0, Pos<T>& map, const T& pnorm, double step) {
    exp_h2_y(beta0, map, pnorm, step);
    exp_h2_z(beta0, map, pnorm, step);
}

template <typename T>
void prop_h3(const double& beta0, Pos<T>& map, const T& pnorm, double step) {
    exp_h3_x(beta0, map, pnorm, step);
    exp_h3_z(beta0, map, pnorm, step);
}

template <typename T>
void prop_ix(const double& brho, const double& kx, const double& ks, 
             const std::vector<std::vector<double>>& coefs, const std::vector<std::vector<double>>& coefs2,
             const std::vector<std::vector<double>>& coefs3, const std::vector<std::vector<double>>& coefs4,
             Pos<T>& map, double s, int sign, double step) {
    exp_ix_px(brho, kx, ks, coefs, coefs2, coefs3, coefs4, map, s, sign, step);
    exp_ix_py(brho, kx, ks, coefs, coefs2, coefs3, coefs4, map, s, sign, step);

}

template <typename T>
void prop_iy(const double& brho, const double& kx, const double& ks, 
             const std::vector<std::vector<double>>& coefs, const std::vector<std::vector<double>>& coefs2,
             const std::vector<std::vector<double>>& coefs3, const std::vector<std::vector<double>>& coefs4,
             Pos<T>& map, double s, int sign, double step) {
    exp_iy_px(brho, kx, ks, coefs, coefs2, coefs3, coefs4, map, s, sign, step);
    exp_iy_py(brho, kx, ks, coefs, coefs2, coefs3, coefs4, map, s, sign, step);

}

template <typename T>
void prop_step(const double& beta0, const double& brho, const double& kx, const double& ks,
               const std::vector<std::vector<double>>& coefs, const std::vector<std::vector<double>>& coefs2,
               const std::vector<std::vector<double>>& coefs3, const std::vector<std::vector<double>>& coefs4,
               Pos<T>& map, const T& pnorm, double& s, double step) {
    prop_h1(beta0, map, pnorm, s, step);
    prop_iy(brho, kx, ks, coefs, coefs2, coefs3, coefs4, map, s, +1, step);
    prop_h2(beta0, map, pnorm, step);
    prop_iy(brho, kx, ks, coefs, coefs2, coefs3, coefs4, map, s, -1, step);
    prop_ix(brho, kx, ks, coefs, coefs2, coefs3, coefs4, map, s, +1, step);
    prop_h3(beta0, map, pnorm, step);
    prop_ix(brho, kx, ks, coefs, coefs2, coefs3, coefs4, map, s, -1, step);
    prop_iy(brho, kx, ks, coefs, coefs2, coefs3, coefs4, map, s, +1, step);
    prop_h2(beta0, map, pnorm, step);
    prop_iy(brho, kx, ks, coefs, coefs2, coefs3, coefs4, map, s, -1, step);
    prop_h1(beta0, map, pnorm, s, step);
}
