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

#ifndef _POS_H
#define _POS_H

#include <vector>
#include <iostream>

template <typename T = double>
class Pos {
public:
	T rx;
	T px;
	T ry;
	T py;
	T de;
	T dl;
	typedef T type;
	//Pos(const Pos<T>& v_);
	Pos(const T& rx_, const T& px_, const T& ry_, const T& py_, const T& de_, const T& dl_);
	Pos(const T& v_ = 0);
	Pos& operator+=(const Pos<T>& v) {
		this->rx += v.rx; this->px += v.px;
		this->ry += v.ry; this->py += v.py;
		this->de += v.de; this->de += v.de;
		return *this;
	}
};


template <typename T>
std::ostream& operator<< (std::ostream &out, const Pos<T>& el) {
	char buffer[255];
	sprintf(buffer, "%+.16e\n", el.rx); out << buffer;
	sprintf(buffer, "%+.16e\n", el.px); out << buffer;
	sprintf(buffer, "%+.16e\n", el.ry); out << buffer;
	sprintf(buffer, "%+.16e\n", el.py); out << buffer;
	sprintf(buffer, "%+.16e\n", el.de); out << buffer;
	sprintf(buffer, "%+.16e\n", el.dl); out << buffer;
	return out;
}

template <typename T>
Pos<T>::Pos(const T& rx_, const T& px_, const T& ry_, const T& py_, const T& de_, const T& dl_):
rx(rx_), px(px_),
ry(ry_), py(py_),
de(de_), dl(dl_)
{

}

// template <typename T>
// Pos<T>::Pos(const Pos<T>& v_):
// rx(v_.rx), px(v_.px),
// ry(v_.ry), py(v_.py),
// de(v_.de), dl(v_.dl)
// {
//
// }

template <typename T>
Pos<T>::Pos(const T& v_):
rx(v_), px(v_),
ry(v_), py(v_),
de(v_), dl(v_)
{

}

template <typename T>
inline
Pos<T> operator+ (const Pos<T>& v1, const Pos<T>& v2) {
	return Pos<T>(v1.rx + v2.rx, v1.px + v2.px, v1.ry + v2.ry, v1.py + v2.py, v1.de + v2.de, v1.dl + v2.dl);
}


template <typename T>
inline
std::vector<Pos<T> > operator+ (const std::vector<Pos<T> >& m1, const std::vector<Pos<T> >& m2) {
	std::vector<Pos<T> > r;
	for(unsigned int i=0; i<m1.size(); ++i) {
		r.push_back(m1[i] + m2[i]);
	}
	return r;
}

template <typename T>
inline
std::vector<Pos<T> > operator- (const std::vector<Pos<T> >& m1, const std::vector<Pos<T> >& m2) {
	std::vector<Pos<T> > r;
	for(unsigned int i=0; i<m1.size(); ++i) {
		r.push_back(m1[i] - m2[i]);
	}
	return r;
}

template <typename T>
inline
Pos<T> operator- (const Pos<T>& v1, const Pos<T>& v2) {
	return Pos<T>(v1.rx - v2.rx, v1.px - v2.px, v1.ry - v2.ry, v1.py - v2.py, v1.de - v2.de, v1.dl - v2.dl);
}

template <typename T>
inline
Pos<T> operator* (const Pos<T>& v1, const T& v) {
	return Pos<T>(v*v1.rx, v*v1.px, v*v1.ry, v*v1.py, v*v1.de, v*v1.dl);
}

template <typename T>
inline
Pos<T> operator* (const T& v, const Pos<T>& v1) {
	return Pos<T>(v*v1.rx, v*v1.px, v*v1.ry, v*v1.py, v*v1.de, v*v1.dl);
}

template <typename T>
inline
std::vector<Pos<T> > operator* (const T& v, const std::vector<Pos<T> >& m1) {
	std::vector<Pos<T> > r;
	for(unsigned int i=0; i<m1.size(); ++i) {
		r.push_back(v * m1[i]);
	}
	return r;
}

template <typename T>
inline
Pos<T> operator/ (const Pos<T>& v1, const T& v) {
	return (1/v)*v1;
}


template <typename T>
T abs(const T& x) {
	if (x>0) return x; else return -x;
}

template <typename T>
T get_max(const Pos<T>& v) {
	T max = abs(v.rx);
	if (abs(v.px) > max) max = abs(v.px);
	if (abs(v.ry) > max) max = abs(v.ry);
	if (abs(v.py) > max) max = abs(v.py);
	if (abs(v.de) > max) max = abs(v.de);
	if (abs(v.dl) > max) max = abs(v.dl);
	//std::cout << max << std::endl;
	return max;
}


#endif
