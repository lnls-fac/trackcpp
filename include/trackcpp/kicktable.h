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
  static std::vector<Kicktable> kicktable_list;

  Kicktable(
    const double x_min,
    const double x_max,
    const unsigned int x_nrpts,
    const std::vector<double>& x_kick,
    const double y_min,
    const double y_max,
    const unsigned int y_nrpts,
    const std::vector<double>& y_kick,
    const double length = 1
  );
  Kicktable(const std::string& filename_ = "");
  Kicktable(const Kicktable &) = default;

  Status::type load_from_file(const std::string& filename_);

  bool is_valid_kicktable() const;
  unsigned int get_idx(unsigned int ix, unsigned int iy) const
  {
    return iy*x_nrpts+ix;
  }
  double get_x(unsigned int ix) const
  {
    return x_min + ix * (x_max - x_min) / (x_nrpts - 1.0);
  }
  double get_y(unsigned int iy) const
  {
    return y_min + iy * (y_max - y_min) / (y_nrpts - 1.0);
  }

  template <typename T>
  unsigned int get_ix(const T& x) const
  {
    return (int) ((x - x_min) / ((x_max - x_min) / (x_nrpts - 1)));
  }
  template <typename T>
  unsigned int get_iy(const T& y) const
  {
    return (int) ((y - y_min) / ((y_max - y_min) / (y_nrpts - 1)));
  }

  template <typename T>
  Status::type getkicks_bilinear(
    const T& rx, const T& ry, T& hkick, T& vkick
  ) const
  {
    // checks x limits
    if ((rx < x_min) or (rx > x_max)) {
      //std::cout << "rx: " << double(rx) << std::endl;
      hkick = nan("");
      return Status::kicktable_out_of_range;
    }
    // checks y limits
    if ((ry < y_min) or (ry > y_max)) {
      //std::cout << "ry: " << double(ry) << std::endl;
      vkick = nan("");
      return Status::kicktable_out_of_range;
    }

    // gets indices
    const unsigned int ix = get_ix(rx);
    const unsigned int iy = get_iy(ry);

    /* coordinates */
    const double x1 = get_x(ix);
    const double x2 = get_x(ix+1);
    const double y1 = get_y(iy);
    const double y2 = get_y(iy+1);

    // Bilinear interpolation
    // https://en.wikipedia.org/wiki/Bilinear_interpolation
    {  /* hkick */
      const double fq11 = x_kick[get_idx(ix+0, iy+0)];
      const double fq12 = x_kick[get_idx(ix+0, iy+1)];
      const double fq21 = x_kick[get_idx(ix+1, iy+0)];
      const double fq22 = x_kick[get_idx(ix+1, iy+1)];
      const T fR1 = ((x2 - rx) * fq11 + (rx - x1) * fq21) / (x2 - x1);
      const T fR2 = ((x2 - rx) * fq12 + (rx - x1) * fq22) / (x2 - x1);
      const T fP  = ((y2 - ry) * fR1  + (ry - y1) * fR2 ) / (y2 - y1);
      hkick = fP;
    }
    {  /* vkick */
      const double fq11 = y_kick[get_idx(ix+0, iy+0)];
      const double fq12 = y_kick[get_idx(ix+0, iy+1)];
      const double fq21 = y_kick[get_idx(ix+1, iy+0)];
      const double fq22 = y_kick[get_idx(ix+1, iy+1)];
      const T fR1 = ((x2 - rx) * fq11 + (rx - x1) * fq21) / (x2 - x1);
      const T fR2 = ((x2 - rx) * fq12 + (rx - x1) * fq22) / (x2 - x1);
      const T fP  = ((y2 - ry) * fR1  + (ry - y1) * fR2 ) / (y2 - y1);
      vkick = fP;
    }

    return Status::success;
  }

  template <typename T>
  Status::type getkicks(const T& rx, const T& ry, T& hkick, T& vkick) const
  {
    return getkicks_bilinear(rx, ry, hkick, vkick);
  }

  static int add_kicktable(
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
  static int add_kicktable(const std::string& filename);
  static int add_kicktable(const Kicktable &new_kicktable);
  static Kicktable& get_kicktable(const int& kicktable_idx)
  {
    if (not is_valid_kicktable_index(kicktable_idx))
      throw std::out_of_range("kicktable_idx is out of range.");
    return kicktable_list[kicktable_idx];
  }
  static bool is_valid_kicktable_index(const int idx)
  {
    return ((idx>0) & (idx < kicktable_list.size()));
  }
  static void clear_kicktables();

  bool operator==(const Kicktable& o) const;
  bool operator!=(const Kicktable& o) const { return !(*this == o); }

};

#endif
