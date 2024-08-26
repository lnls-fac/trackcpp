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
// #include "../alglib/interpolation.h"
#include <string>
#include <vector>
#include <algorithm>


class Kicktable {

  // Kicktable: assumes x_kick and y_kick tables sampled on the same regular (x,y)-point grid.

public:

  std::string         filename;
  double              length;
  std::vector<double> x_pos, y_pos;
  std::vector<double> x_kick,  y_kick;
  static std::vector<Kicktable> kicktable_list;

  Kicktable(
    const std::vector<double>& x_pos,
    const std::vector<double>& x_kick,
    const std::vector<double>& y_pos,
    const std::vector<double>& y_kick,
    const double length = 1
  );
  Kicktable(const std::string& filename_ = "");
  Kicktable(const Kicktable &) = default;

  Status::type load_from_file(const std::string& filename_);

  bool is_valid_kicktable() const;
  unsigned int get_idx(unsigned int ix, unsigned int iy) const
  {
    return iy * x_pos.size() + ix;
  }
  double get_x(unsigned int ix) const {return x_pos[ix];}
  double get_y(unsigned int iy) const {return y_pos[iy];}
  template <typename T> unsigned int get_ix(const T& x) const
  {
    return std::lower_bound(x_pos.begin(), x_pos.end(), x) - x_pos.begin();
  }
  template <typename T> unsigned int get_iy(const T& y) const
  {
    return std::lower_bound(y_pos.begin(), y_pos.end(), y) - y_pos.begin();
  }

  template <typename T> Status::type getkicks_bilinear(
    const T& rx, const T& ry, T& hkick, T& vkick
  ) const
  {
    // gets indices
    unsigned int ix = get_ix(rx);
    unsigned int ixp1, iyp1;
    if (ix >= x_pos.size()-1){ixp1 = ix; --ix;}
    else ixp1 = ix + 1;

    unsigned int iy = get_iy(ry);
    if (iy >= y_pos.size()-1){iyp1 = iy; --iy;}
    else iyp1 = iy + 1;

    /* coordinates */
    const double x1 = get_x(ix);
    const double x2 = get_x(ixp1);
    const double y1 = get_y(iy);
    const double y2 = get_y(iyp1);
    double denom = (x2 - x1) * (y2 - y1);

    // Bilinear interpolation
    // https://en.wikipedia.org/wiki/Bilinear_interpolation
    {  /* hkick */
      const double f11 = x_kick[get_idx(ix, iy)];
      const double f12 = x_kick[get_idx(ix, iyp1)];
      const double f21 = x_kick[get_idx(ixp1, iy)];
      const double f22 = x_kick[get_idx(ixp1, iyp1)];
      hkick =
        f11 * (x2 - rx) * (y2 - ry) + f21 * (rx - x1) * (y2 - ry) +
        f12 * (x2 - rx) * (ry - y1) + f22 * (rx - x1) * (ry - y1);
      hkick /= denom;
    }
    {  /* vkick */
      const double f11 = y_kick[get_idx(ix, iy)];
      const double f12 = y_kick[get_idx(ix, iyp1)];
      const double f21 = y_kick[get_idx(ixp1, iy)];
      const double f22 = y_kick[get_idx(ixp1, iyp1)];
      vkick =
        f11 * (x2 - rx) * (y2 - ry) + f21 * (rx - x1) * (y2 - ry) +
        f12 * (x2 - rx) * (ry - y1) + f22 * (rx - x1) * (ry - y1);
      vkick /= denom;
    }

    return Status::success;
  }

  template <typename T> Status::type getkicks(
    const T& rx, const T& ry, T& hkick, T& vkick
  ) const
  {
    return getkicks_bilinear(rx, ry, hkick, vkick);
  }

  static int add_kicktable(
    const std::vector<double>& x_pos,
    const std::vector<double>& x_kick,
    const std::vector<double>& y_pos,
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
