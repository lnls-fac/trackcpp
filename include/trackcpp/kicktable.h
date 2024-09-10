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

public:

  std::string         filename;
  double              length;
  std::vector<double> x_pos, y_pos;
  std::vector<double> x_kick,  y_kick;
  static std::vector<Kicktable> kicktable_list;

  Kicktable(
    const std::vector<double>& x_pos,
    const std::vector<double>& y_pos,
    const std::vector<double>& x_kick,
    const std::vector<double>& y_kick,
    const double length = 1
  );
  Kicktable(const std::string& filename_ = "");
  Kicktable(const Kicktable &) = default;

  Status::type load_from_file(
    const std::string& filename_, bool file_flag=true
  );
  Status::type _dump_to_stream(
    std::ostream& fp, const std::string author_name = " "
  );
  Status::type save_to_file(
    const std::string filename_, const std::string author_name = " "
  );
  std::string save_to_string(const std::string author_name = " ");

  bool is_valid_kicktable() const;
  unsigned int get_idx(unsigned int ix, unsigned int iy) const
  {
    return iy * x_pos.size() + ix;
  }
  double get_x(unsigned int ix) const {return x_pos[ix];}
  double get_y(unsigned int iy) const {return y_pos[iy];}
  template <typename T> unsigned int get_ix(const T& x) const
  {
    int idx;
    idx = std::lower_bound(x_pos.begin(), x_pos.end(), x) - x_pos.begin() - 1;
    return idx >= 0 ? (unsigned int) idx : 0u;
  }
  template <typename T> unsigned int get_iy(const T& y) const
  {
    int idx;
    idx = std::lower_bound(y_pos.begin(), y_pos.end(), y) - y_pos.begin() - 1;
    return idx >= 0 ? (unsigned int) idx : 0u;
  }
  unsigned int get_ix(const double& x) const {return get_ix<double>(x);}
  unsigned int get_iy(const double& y) const {return get_iy<double>(y);}

  template <typename T> Status::type getkicks_bilinear(
    const T& rx, const T& ry, T& hkick, T& vkick
  ) const
  {

    // gets indices
    unsigned int ix = get_ix(rx);
    unsigned int iy = get_iy(ry);
    unsigned int ixp1, iyp1;
    if (ix >= x_pos.size()-1) ixp1 = ix--;
    else ixp1 = ix + 1;
    if (iy >= y_pos.size()-1) iyp1 = iy--;
    else iyp1 = iy + 1;
    const unsigned int i00 = get_idx(ix, iy);
    const unsigned int i01 = get_idx(ix, iyp1);
    const unsigned int i10 = get_idx(ixp1, iy);
    const unsigned int i11 = get_idx(ixp1, iyp1);

    /* coordinates */
    const double x1 = get_x(ix);
    const double x2 = get_x(ixp1);
    const double y1 = get_y(iy);
    const double y2 = get_y(iyp1);
    const double denom = (x2 - x1) * (y2 - y1);
    const T dx1 = rx - x1;
    const T dx2 = x2 - rx;
    const T dy1 = ry - y1;
    const T dy2 = y2 - ry;

    // Bilinear interpolation
    // https://en.wikipedia.org/wiki/Bilinear_interpolation
    {  /* hkick */
      const double& f11 = x_kick[i00];
      const double& f12 = x_kick[i01];
      const double& f21 = x_kick[i10];
      const double& f22 = x_kick[i11];
      hkick = (f11 * dx2 + f21 * dx1) * dy2 + (f12 * dx2 + f22 * dx1) * dy1;
      hkick /= denom;
    }
    {  /* vkick */
      const double& f11 = y_kick[i00];
      const double& f12 = y_kick[i01];
      const double& f21 = y_kick[i10];
      const double& f22 = y_kick[i11];
      vkick = (f11 * dx2 + f21 * dx1) * dy2 + (f12 * dx2 + f22 * dx1) * dy1;
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

  Status::type getkicks(
    const double& rx, const double& ry, double& hkick__, double& vkick__
  ) const
  {
    return getkicks<double>(rx, ry, hkick__, vkick__);
  }

  static int add_kicktable(
    const std::vector<double>& x_pos,
    const std::vector<double>& y_pos,
    const std::vector<double>& x_kick,
    const std::vector<double>& y_kick,
    const double length=1
  );
  static int add_kicktable(const std::string filename);
  static int add_kicktable(const Kicktable &new_kicktable);
  static Kicktable& get_kicktable(const int& kicktable_idx)
  {
    if (not is_valid_kicktable_index(kicktable_idx))
      throw std::out_of_range("kicktable_idx is out of range.");
    return kicktable_list[kicktable_idx];
  }
  static bool is_valid_kicktable_index(const int idx)
  {
    return ((idx>=0) & (idx < kicktable_list.size()));
  }
  static void clear_kicktables();
  static size_t get_kicktable_list_size(){return kicktable_list.size();}

  bool operator==(const Kicktable& o) const;
  bool operator!=(const Kicktable& o) const { return !(*this == o); }

};

#endif
