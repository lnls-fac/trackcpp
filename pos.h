#ifndef _POS_H
#define _POS_H

// TRACKC++
// ========
// Author: 		Ximenes R. Resende
// email:  		xresende@gmail.com, ximenes.resende@lnls.br
// affiliation:	LNLS - Laboratorio Nacional de Luz Sincrotron
// Date: 		Tue Dec 10 17:57:20 BRST 2013

#include <vector>

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
Pos<T>::Pos(const T& rx_, const T& px_, const T& ry_, const T& py_, const T& de_, const T& dl_):
rx(rx_), px(px_),
ry(ry_), py(py_),
de(de_), dl(dl_)
{

}

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
inline
void vector6_sum(Pos<T>& v, const T& a1, const Pos<T>& v1, const T& a2, const Pos<T>& v2) {
	unsigned int i = 0;
	v.rx = a1 * v1.rx + a2 * v2.rx; v.px = a1 * v1.px + a2 * v2.px;
	v.ry = a1 * v1.ry + a2 * v2.ry; v.py = a1 * v1.py + a2 * v2.py;
	v.de = a1 * v1.de + a2 * v2.de; v.de = a1 * v1.de + a2 * v2.de;
}

template <typename T>
inline
void matrix6_set_identity(std::vector<Pos<T> >& m, const T& a = 1) {
	for(unsigned int i=0; i<m.size(); ++i) m[i] = Pos<T>(0,0,0,0,0,0);
	m[0].rx = m[1].px = m[2].ry = m[3].py = m[4].de = m[5].dl = a;
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
