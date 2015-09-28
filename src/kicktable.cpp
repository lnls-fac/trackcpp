// The MIT License (MIT)
//
// Copyright (c) 2015 LNLS Accelerator Division
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

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
  if (fp.fail()) return Status::file_not_found;
  this->filename = filename_;

  std::string str;

  // HEADER
  getline(fp, str);   // name of kicktable line
  getline(fp, str);   // author line
  getline(fp, str);   // label 'ID length[m]'
  fp >> this->length; // length of element
  getline(fp, str);   // advances to new line
  getline(fp, str);   // label 'number of horizontal points'
  fp >> this->x_nrpts;      // number of horizontal points
  getline(fp, str);   // advances to new line
  getline(fp, str);   // label 'number of vertical points'
  fp >> this->y_nrpts;      // number of vertical points
  getline(fp, str);   // advances to new line

  this->x_kick.resize(x_nrpts * y_nrpts, 0);
  this->y_kick.resize(x_nrpts * y_nrpts, 0);

  // HORIZONTAL KICK TABLE
  getline(fp, str);   // label 'Horizontal KickTable in T^2.m^2'
  getline(fp, str);   // label 'START'
  for(unsigned int i=0; i<this->x_nrpts; ++i) {
    double posx; fp >> posx;
    if (isnan(x_min) or posx < x_min) x_min = posx;
    if (isnan(x_max) or posx > x_max) x_max = posx;
  }
  for(int j=y_nrpts-1; j>=0; --j) {
    double posy; fp >> posy;
    if (isnan(y_min) or posy < y_min) y_min = posy;
    if (isnan(y_max) or posy > y_max) y_max = posy;
    for(unsigned int i=0; i<x_nrpts; ++i) fp >> x_kick[this->get_idx(i,j)];
  }
  getline(fp, str);   // advances to new line

  // VERTICAL KICK TABLE
  getline(fp, str);   // label 'Vertical KickTable in T^2.m^2'
  getline(fp, str);   // label 'START'
  for(unsigned int i=0; i<this->x_nrpts; ++i) { double posx; fp >> posx; }
  for(int j=y_nrpts-1; j>=0; --j) {
    double posy; fp >> posy;
    for(unsigned int i=0; i<x_nrpts; ++i) fp >> y_kick[this->get_idx(i,j)];
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

Kick Kicktable::interpolate_kicks(const double& rx, const double& ry) const {
    Kick kick;
    kick.status = kicktable_getkicks(this, rx, ry, kick.hkick, kick.vkick);
    return kick;
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
