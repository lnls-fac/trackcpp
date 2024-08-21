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

#ifndef _KICKTABLE_H
#define _KICKTABLE_H

#include "auxiliary.h"
#include <string>
#include <vector>


class Kicktable {

  // Kicktable: assumes x_kick and y_kick tables sampled on the same regular (x,y)-point grid.

public:

  std::string         filename;
  double              length;
  unsigned int        x_nrpts, y_nrpts;
  double              x_min, x_max;
  double              y_min, y_max;
  std::vector<double> x_kick,  y_kick;

  Kicktable(const std::string& filename_ = "");
  Kicktable(const Kicktable &) = default;

  Status::type load_from_file(const std::string& filename_);

  unsigned int get_idx(unsigned int ix, unsigned int iy) const { return iy*x_nrpts+ix; }
  double       get_x(unsigned int ix) const { return x_min + ix * (x_max - x_min) / (x_nrpts - 1.0); }
  double       get_y(unsigned int iy) const { return y_min + iy * (y_max - y_min) / (y_nrpts - 1.0); }
  template <typename T> unsigned int get_ix(const T& x) const { return (int) ((x - x_min) / ((x_max - x_min) / (x_nrpts - 1))); }
  template <typename T> unsigned int get_iy(const T& y) const { return (int) ((y - y_min) / ((y_max - y_min) / (y_nrpts - 1))); }

  bool operator==(const Kicktable& o) const;
  bool operator!=(const Kicktable& o) const { return !(*this == o); }

};


extern std::vector<Kicktable> kicktable_list;


int add_kicktable(
  const double x_min,
  const double x_max,
  const unsigned int x_nrpts,
  const std::vector<double>& x_kick,
  const double y_min,
  const double y_max,
  const unsigned int y_nrpts,
  const std::vector<double>& y_kick,
  const double length=1
);
int add_kicktable(const std::string& filename);
void clear_kicktables(std::vector<Kicktable>& kicktable_list);

template <typename T>
Status::type kicktable_getkicks_bilinear(const int& kicktable_idx, const T& rx, const T& ry, T& hkick, T& vkick) {

  const Kicktable* kicktable = &kicktable_list[kicktable_idx];

  // checks x limits
  const double& xmin = kicktable->x_min;
  const double& xmax = kicktable->x_max;
  if ((rx < xmin) or (rx > xmax)) {
    //std::cout << "rx: " << double(rx) << std::endl;
    hkick = nan("");
    return Status::kicktable_out_of_range;
  }

  // checks y limits
  const double& ymin = kicktable->y_min;
  const double& ymax = kicktable->y_max;
  if ((ry < ymin) or (ry > ymax)) {
    //std::cout << "ry: " << double(ry) << std::endl;
    vkick = nan("");
    return Status::kicktable_out_of_range;
  }

  // gets indices
  const unsigned int ix = kicktable->get_ix(rx);
  const unsigned int iy = kicktable->get_iy(ry);

  /* coordinates */
  const double x1  = kicktable->get_x(ix);
  const double x2  = kicktable->get_x(ix+1);
  const double y1  = kicktable->get_y(iy);
  const double y2  = kicktable->get_y(iy+1);

  // Bilinear interpolation - https://en.wikipedia.org/wiki/Bilinear_interpolation

  {  /* hkick */
    const double fq11 = kicktable->x_kick[kicktable->get_idx(ix+0, iy+0)];
    const double fq12 = kicktable->x_kick[kicktable->get_idx(ix+0, iy+1)];
    const double fq21 = kicktable->x_kick[kicktable->get_idx(ix+1, iy+0)];
    const double fq22 = kicktable->x_kick[kicktable->get_idx(ix+1, iy+1)];
    const T fR1 = ((x2 - rx) * fq11 + (rx - x1) * fq21) / (x2 - x1);
    const T fR2 = ((x2 - rx) * fq12 + (rx - x1) * fq22) / (x2 - x1);
    const T fP  = ((y2 - ry) * fR1  + (ry - y1) * fR2 ) / (y2 - y1);
    hkick = fP;
  }
  {  /* vkick */
    const double fq11 = kicktable->y_kick[kicktable->get_idx(ix+0, iy+0)];
    const double fq12 = kicktable->y_kick[kicktable->get_idx(ix+0, iy+1)];
    const double fq21 = kicktable->y_kick[kicktable->get_idx(ix+1, iy+0)];
    const double fq22 = kicktable->y_kick[kicktable->get_idx(ix+1, iy+1)];
    const T fR1 = ((x2 - rx) * fq11 + (rx - x1) * fq21) / (x2 - x1);
    const T fR2 = ((x2 - rx) * fq12 + (rx - x1) * fq22) / (x2 - x1);
    const T fP  = ((y2 - ry) * fR1  + (ry - y1) * fR2 ) / (y2 - y1);
    vkick = fP;
  }

  return Status::success;
}

template <typename T>
Status::type kicktable_getkicks(const int& kicktable_idx, const T& rx, const T& ry, T& hkick, T& vkick) {
  return kicktable_getkicks_bilinear(kicktable_idx, rx, ry, hkick, vkick);
}

#endif
