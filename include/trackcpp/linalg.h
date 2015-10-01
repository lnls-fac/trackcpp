#ifndef _LINALG_H
#define _LINALG_H

#include <trackcpp/auxiliary.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

#include "pos.h"

typedef double Vector2[2];
typedef double Matrix2[2][2];

void linalg_solve2(Vector2&, const Matrix2& m, const Vector2& b);

Pos<double> linalg_solve2_posvec(const std::vector<Pos<double> >& M, const Pos<double>& B);
Pos<double> linalg_solve4_posvec(const std::vector<Pos<double> >& M, const Pos<double>& B);
Pos<double> linalg_solve6_posvec(const std::vector<Pos<double> >& M, const Pos<double>& B);

void getmx(const Matrix& m, Matrix2& mx);
void getmy(const Matrix& m, Matrix2& my);

void matrix2_eye(Matrix2& m);
void matrix2_lc(Matrix2& m, const double& a1, const Matrix2& m1, const double& a2, const Matrix2& m2);

template <typename T>
inline
void vector6_lc_posvec(Pos<T>& v, const T& a1, const Pos<T>& v1, const T& a2, const Pos<T>& v2) {
	unsigned int i = 0;
	v.rx = a1 * v1.rx + a2 * v2.rx; v.px = a1 * v1.px + a2 * v2.px;
	v.ry = a1 * v1.ry + a2 * v2.ry; v.py = a1 * v1.py + a2 * v2.py;
	v.de = a1 * v1.de + a2 * v2.de; v.de = a1 * v1.de + a2 * v2.de;
}

template <typename T>
inline
void matrix6_set_identity_posvec(std::vector<Pos<T> >& m, const T& a = 1) {
	for(unsigned int i=0; i<m.size(); ++i) m[i] = Pos<T>(0,0,0,0,0,0);
	m[0].rx = m[1].px = m[2].ry = m[3].py = m[4].de = m[5].dl = a;
}

template <typename T>
inline
void matrix6_lc_posvec(std::vector<Pos<T> >& m, const T& a1, const std::vector<Pos<T>>& m1, const T& a2, const std::vector<Pos<T>>& m2) {
	m = m1;
	vector6_sum(m[0], a1, m1[0], a2, m2[0]);
	vector6_sum(m[1], a1, m1[1], a2, m2[1]);
	vector6_sum(m[2], a1, m1[2], a2, m2[2]);
	vector6_sum(m[3], a1, m1[3], a2, m2[3]);
	vector6_sum(m[4], a1, m1[4], a2, m2[4]);
	vector6_sum(m[5], a1, m1[5], a2, m2[5]);
}

#endif
