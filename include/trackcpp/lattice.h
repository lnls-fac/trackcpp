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

#ifndef _LATTICE_H
#define _LATTICE_H

#include "accelerator.h"
#include "elements.h"
#include "auxiliary.h"
#include <vector>

std::vector<Element> latt_join(const std::vector<std::vector<Element> >& v_);
std::vector<Element> latt_reverse(const std::vector<Element>& v_);
void                 latt_print(const std::vector<Element>& lattice);
std::vector<int>     latt_range(const std::vector<Element>& lattice);
std::vector<double>  latt_findspos(const std::vector<Element>& lattice, const std::vector<int>& idx);
double               latt_findspos(const std::vector<Element>& lattice, const int idx);
void                 latt_setcavity(std::vector<Element>& lattice, const std::string& state);
std::vector<Element> latt_set_num_integ_steps(const std::vector<Element>& orig_lattice);
std::vector<Element> latt_read_flat_file(const std::string& filename);
Status::type         latt_read_flat_file(const std::string& filename, Accelerator& accelerator);
std::vector<int>     latt_findcells_fam_name    (const std::vector<Element>& lattice, const std::string& value, bool reverse = false);
std::vector<int>     latt_findcells_angle       (const std::vector<Element>& lattice, const double& value,      bool reverse = false);
std::vector<int>     latt_findcells_frequency   (const std::vector<Element>& lattice, const double& value,      bool reverse = false);
std::vector<int>     latt_findcells_polynom_b   (const std::vector<Element>& lattice, unsigned int n, const double& value, bool reverse = false);
std::vector<int>     latt_findcells_polynom_a   (const std::vector<Element>& lattice, unsigned int n, const double& value, bool reverse = false);
std::vector<int>     latt_findcells_pass_method (const std::vector<Element>& lattice, const std::string& value, bool reverse = false);

template <typename T>
std::vector<T> latt_getcellstruct(const std::vector<Element>& lattice, const std::string& field, std::vector<int> idx) {
	std::vector<T> r;
	if (field == "fam_name") {
		for(unsigned int i=0; i<idx.size(); ++i) {
			r.push_back(*((T*)(&(lattice[idx[i]].fam_name))));
		}
	} else if (field == "length") {
		for(unsigned int i=0; i<idx.size(); ++i) {
			r.push_back(*((T*)(&(lattice[idx[i]].length))));
		}
	}
	return r;
}
#endif
