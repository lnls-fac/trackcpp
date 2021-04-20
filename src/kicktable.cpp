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

Kicktable::Kicktable(const std::string& filename_) :
  filename(""),
  x_nrpts(0), y_nrpts(0),
  x_min(nan("")), x_max(nan("")),
  y_min(nan("")), y_max(nan("")) {

  if (filename_ != "") {
    this->load_from_file(filename_);
  }

}

Status::type Kicktable::load_from_file(const std::string& filename_) {

  std::ifstream fp(filename_);
  if (fp.fail()) {
    std::cout << "Could not find kicktable file " << filename_ << "!" << std::endl;
    return Status::file_not_found;
  }
  this->filename = filename_;

  std::string str;

  // HEADER
  getline(fp, str);    // name of kicktable line
  getline(fp, str);    // author line
  getline(fp, str);    // label 'ID length[m]'
  fp >> this->length;  // length of element
  getline(fp, str);    // advances to new line
  getline(fp, str);    // label 'number of horizontal points'
  fp >> this->x_nrpts; // number of horizontal points
  getline(fp, str);    // advances to new line
  getline(fp, str);    // label 'number of vertical points'
  fp >> this->y_nrpts; // number of vertical points
  getline(fp, str);    // advances to new line

  this->x_kick.resize(x_nrpts * y_nrpts, 0);
  this->y_kick.resize(x_nrpts * y_nrpts, 0);

  std::vector<double> yvec(y_nrpts); // used to invert tables ordering, if necessary

  // HORIZONTAL KICK TABLE
  getline(fp, str);   // label 'Horizontal KickTable in T^2.m^2'
  getline(fp, str);   // label 'START'
  for(unsigned int i=0; i<this->x_nrpts; ++i) {
    double posx; fp >> posx;
    if (std::isnan(x_min) or posx < x_min) x_min = posx;
    if (std::isnan(x_max) or posx > x_max) x_max = posx;
  }
  for(int j=y_nrpts-1; j>=0; --j) {
    double posy; fp >> posy; yvec.push_back(posy);
    if (std::isnan(y_min) or posy < y_min) y_min = posy;
    if (std::isnan(y_max) or posy > y_max) y_max = posy;
    for(unsigned int i=0; i<x_nrpts; ++i)
      fp >> x_kick[this->get_idx(i,j)];
  }
  getline(fp, str);   // advances to new line

  // VERTICAL KICK TABLE
  getline(fp, str);   // label 'Vertical KickTable in T^2.m^2'
  getline(fp, str);   // label 'START'
  for(unsigned int i=0; i<this->x_nrpts; ++i) { double posx; fp >> posx; }
  for(int j=y_nrpts-1; j>=0; --j) {
    double posy; fp >> posy;
    for(unsigned int i=0; i<x_nrpts; ++i)
      fp >> y_kick[this->get_idx(i,j)];
  }

  // invert tables, if necessary
  if (yvec.size() > 1 and yvec[1] > yvec[0]) {
    for(unsigned int i=0; i<this->x_nrpts; ++i) {
      for(unsigned int j=0; j<this->y_nrpts/2; ++j) {
        std::swap(x_kick[this->get_idx(i,j)], x_kick[this->get_idx(i,this->y_nrpts-j-1)]);
        std::swap(y_kick[this->get_idx(i,j)], y_kick[this->get_idx(i,this->y_nrpts-j-1)]);
      }
    }
  }
  return Status::success;

}

Status::type add_kicktable(const std::string& filename, std::vector<Kicktable*>& kicktable_list, const Kicktable*& kicktable_pointer) {

  // looks through vector of kickmaps...
  for(unsigned int i=0; i<kicktable_list.size(); ++i) {
    if (kicktable_list[i]->filename == filename) {
      kicktable_pointer = kicktable_list[i];
      return Status::success;
    }
  }

  // loads a new kicktable from file and inserts it into vector of kicktables
  Kicktable* new_kicktable = new Kicktable();
  Status::type status = new_kicktable->load_from_file(filename);
  if (status == Status::success) {
    kicktable_list.push_back(new_kicktable);
    kicktable_pointer = new_kicktable;
  } else {
    kicktable_pointer = nullptr;
  }
  return status;

}


void del_kicktables(std::vector<Kicktable*>& kicktable_list) {
  for(unsigned int i=0; i<kicktable_list.size(); ++i) {
    delete [] kicktable_list[i];
  }
}

bool Kicktable::operator==(const Kicktable& o) const {
  if (this == &o) return true;
  if (this->length != o.length) return false;
  if (this->x_min != o.x_min) return false;
  if (this->x_max != o.x_max) return false;
  if (this->y_min != o.y_min) return false;
  if (this->y_max != o.y_max) return false;
  if (this->x_kick != o.x_kick) return false;
  if (this->y_kick != o.y_kick) return false;
  return true;
}
