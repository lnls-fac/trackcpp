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

#include <trackcpp/kicktable.h>
#include <trackcpp/auxiliary.h>
#include <string>
#include <fstream>
#include <cmath>


std::vector<Kicktable> Kicktable::kicktable_list;

Kicktable::Kicktable(const std::string& filename_) :
  filename("")
{
  if (filename_ != "")
    this->load_from_file(filename_);
}


Kicktable::Kicktable(
  const std::vector<double>& x_pos,
  const std::vector<double>& y_pos,
  const std::vector<double>& x_kick,
  const std::vector<double>& y_kick,
  const double length
) :
  filename(""),
  x_pos(x_pos), y_pos(y_pos),
  x_kick(x_kick), y_kick(y_kick),
  length(length)
{}

bool Kicktable::is_valid_kicktable() const
{

  if (x_pos.size() <= 1)
    return false;
  if (y_pos.size() <= 1)
    return false;
  if (x_kick.size() != x_pos.size() * y_pos.size())
    return false;
  if (y_kick.size() != x_kick.size())
    return false;
  if (not std::is_sorted(x_pos.cbegin(), x_pos.cend()))
    return false;
  if (not std::is_sorted(y_pos.cbegin(), y_pos.cend()))
    return false;
  return true;
}

Status::type Kicktable::load_from_file(const std::string& filename_) {

  std::ifstream fp(filename_);
  if (fp.fail()) {
    std::cout << "Could not find kicktable file " << filename_ << "!" << std::endl;
    return Status::file_not_found;
  }
  filename = filename_;

  std::string str;
  unsigned int x_nrpts, y_nrpts;

  // HEADER
  getline(fp, str);   // name of kicktable line
  getline(fp, str);   // author line
  getline(fp, str);   // label 'ID length[m]'
  fp >> length; // length of element
  getline(fp, str);   // advances to new line
  getline(fp, str);   // label 'number of horizontal points'
  fp >> x_nrpts;      // number of horizontal points
  getline(fp, str);   // advances to new line
  getline(fp, str);   // label 'number of vertical points'
  fp >> y_nrpts;      // number of vertical points
  getline(fp, str);   // advances to new line

  x_kick.resize(x_nrpts * y_nrpts, 0);
  y_kick.resize(x_nrpts * y_nrpts, 0);
  x_pos.resize(x_nrpts, 0);
  // used to invert tables ordering, if necessary
  y_pos.resize(y_nrpts, 0);

  // HORIZONTAL KICK TABLE
  getline(fp, str);   // label 'Horizontal KickTable in T^2.m^2'
  getline(fp, str);   // label 'START'
  for(unsigned int i=0; i<x_nrpts; ++i) {fp >> x_pos[i];}
  for(int j=y_nrpts-1; j>=0; --j) {
    fp >> y_pos[j];
    for(unsigned int i=0; i<x_nrpts; ++i)
      fp >> x_kick[get_idx(i, j)];
  }
  getline(fp, str);   // advances to new line

  // VERTICAL KICK TABLE
  getline(fp, str);   // label 'Vertical KickTable in T^2.m^2'
  getline(fp, str);   // label 'START'
  getline(fp, str);   // horizontal position. Already set from x_kick
  for(int j=y_nrpts-1; j>=0; --j) {
    double posy; fp >> posy;
    for(unsigned int i=0; i<x_nrpts; ++i)
      fp >> y_kick[get_idx(i, j)];
  }

  // invert tables, if necessary
  if (y_pos.size() > 1 and y_pos[1] > y_pos[0]) {
    for(unsigned int i=0; i<x_nrpts; ++i) {
      for(unsigned int j=0; j<y_nrpts/2; ++j) {
        std::swap(x_kick[get_idx(i, j)], x_kick[get_idx(i, y_nrpts-j-1)]);
        std::swap(y_kick[get_idx(i, j)], y_kick[get_idx(i, y_nrpts-j-1)]);
      }
    }
  }
  return Status::success;
}


static const int hw = 18; // header field width
static const int pw = 16; // parameter field width
static const int np = 17; // number precision
Status::type Kicktable::save_to_file(
  std::string& filename_,
  const std::string author_name,
  const bool file_flag
)
{
  // done with the help of chatgpt:
  std::unique_ptr<std::ostream> fp;
  if (file_flag)
  {
    fp = std::make_unique<std::ofstream>(filename.c_str());
    if (!fp->good())
      return Status::file_not_found;
  }
  else
    fp = std::make_unique<std::stringstream>();

  fp->setf(
    std::ios_base::left |
    std::ios_base::scientific |
    std::ios_base::uppercase
  );
  fp->precision(np);

  // HEADER
  *fp << "# Author: " << author_name << '\n';
  *fp << "#" << '\n';
  *fp << "# Total Length of Longitudinal Interval [m]" << '\n';
  *fp << length << '\n';
  *fp << "# Number of Horizontal Points" << '\n';
  *fp << x_pos.size() << '\n';
  *fp << "# Number of Vertical Points" << '\n';
  *fp << y_pos.size() << '\n';

  // x_kick:
  *fp << "# Total Horizontal 2nd Order Kick [T2m2]" << '\n';
  *fp << "START" << '\n';
  *fp << std::setw(pw) << ' ';
  for (auto xi: x_pos)
    *fp << std::setw(pw) << xi << ' ';
  *fp << '\n';
  for (auto j=y_pos.size()-1; j>=0; ++j)
  {
    *fp << std::setw(pw) << y_pos[j] << ' ';
    for (auto i=0; i<x_pos.size(); ++i)
    {
      *fp << std::setw(pw) << x_kick[get_idx(i, j)] << ' ';
    }
    *fp << '\n';
  }

  // y_kick:
  *fp << "# Total Vertical 2nd Order Kick [T2m2]" << '\n';
  *fp << "START" << '\n';
  *fp << std::setw(pw) << ' ';
  for (auto xi: x_pos)
    *fp << std::setw(pw) << xi << ' ';
  *fp << '\n';
  for (auto j=y_pos.size()-1; j>=0; ++j)
  {
    *fp << std::setw(pw) << y_pos[j] << ' ';
    for (auto i=0; i<x_pos.size(); ++i)
    {
      *fp << std::setw(pw) << y_kick[get_idx(i, j)] << ' ';
    }
    *fp << '\n';
  }

  if (!file_flag)
    // done with the help of chatgpt:
    filename = std::move(static_cast<std::stringstream*>(fp.get())->str());
  return Status::success;

}

int Kicktable::add_kicktable(
  const std::vector<double>& x_pos,
  const std::vector<double>& y_pos,
  const std::vector<double>& x_kick,
  const std::vector<double>& y_kick,
  const double length
)
{
  Kicktable new_kicktable = Kicktable(x_pos, y_pos, x_kick, y_kick, length);
  return Kicktable::add_kicktable(new_kicktable);
}


int Kicktable::add_kicktable(const std::string& filename)
{
  // loads a new kicktable from file and inserts it into vector of kicktables
  Kicktable new_kicktable(filename);
  return Kicktable::add_kicktable(new_kicktable);
}


int Kicktable::add_kicktable(const Kicktable &new_kicktable)
{
  if (not new_kicktable.is_valid_kicktable())
    return -1;

  // looks through vector of kicktables...
  for(int i=0; i<kicktable_list.size(); ++i)
    if (kicktable_list[i] == new_kicktable)
      return i;

  kicktable_list.push_back(new_kicktable);
  return kicktable_list.size() - 1;
}

void Kicktable::clear_kicktables() {
  kicktable_list.clear();
}

bool Kicktable::operator==(const Kicktable& o) const {
  if (this == &o) return true;
  if (this->length != o.length) return false;
  if (this->x_pos != o.x_pos) return false;
  if (this->y_pos != o.y_pos) return false;
  if (this->x_kick != o.x_kick) return false;
  if (this->y_kick != o.y_kick) return false;
  return true;
}
