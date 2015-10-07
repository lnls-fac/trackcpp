#include <limits>
#include <trackcpp/linalg.h>

Vector& Vector::multiplication(const Matrix& m, const Vector& b) {
  Vector& v = *this;
  v.clear();
  for(unsigned int i=0; i<m.size(); i++) {
    v.push_back(0);
    for(unsigned int j=0; j<b.size(); ++j) v[i] += m[i][j] * b[j];
  }
  return v;
}

void    Vector::print() const {
  const Vector& m = *this;
  for(unsigned int i=0; i<m.size(); ++i) printf("%+.4e ", m[i]);
  printf("\n");
}

double Matrix::norm() const {
  double max = 0;
  const Matrix& m = *this;
  for(unsigned int i=0; i<m.size(); ++i)
    for(unsigned int j=0; j<m[i].size(); ++j)
      max = std::max(max, abs(m[i][j]));
  return max;
}

Matrix& Matrix::eye(const double& v) {
  Matrix& m = *this;
  for(unsigned int i=0; i<m.size(); ++i)
    for(unsigned int j=0; j<m[i].size(); ++j)
      if (j==i) m[i][j] = v; else m[i][j] = 0.0;
  return m;
}

Matrix& Matrix::scalar(const double& v) {
  Matrix& m = *this;
  for(unsigned int i=0; i<m.size(); ++i)
    for(unsigned int j=0; j<m[i].size(); ++j)
      m[i][j] *= v;
  return m;
}

Matrix& Matrix::transpose(int size, unsigned int r, unsigned int c) {
  Matrix& m = *this;
  if (size < 0) size = m.size();
  for(unsigned int i=0; i<size; ++i)
    for(unsigned int j=i+1; j<size; ++j)
      std::swap(m[r+i][c+j], m[r+j][c+i]);
  return m;
}

void    Matrix::print() const {
  const Matrix& m = *this;
  for(unsigned int i=0; i<m.size(); ++i) {
    for(unsigned int j=0; j<m[i].size(); ++j) {
      printf("%+.4e ", m[i][j]);
    }
    printf("\n");
  }
}

Matrix& Matrix::linear_combination(const double& a1, const Matrix& m1, const double& a2, const Matrix& m2) {

  Matrix& m = *this;
  m.clear();
  for(unsigned int i=0; i<m1.size(); ++i) {
    Vector v(m2[i].size());
    for(unsigned int j=0; j<v.size(); ++j)
      v[j] = a1 * m1[i][j] + a2 * m2[i][j];
    m.push_back(v);
  }
  return m;
}

Matrix& Matrix::multiplication(const Matrix& m1, const Matrix& m2) {
  Matrix& m = *this;
  unsigned int nr = m1.size();
  unsigned int nc = m2[0].size();
  m.clear();
  for(unsigned int i=0; i<nr; ++i) {
    Vector v(std::vector<double>(nc,0));
    m.push_back(v);
    for(unsigned int j=0; j<nc; ++j)
      for(unsigned int k=0; k<m2.size(); ++k) m[i][j] += m1[i][k] * m2[k][j];
  }
  return m;
}

Matrix& Matrix::getM(Matrix& s, int nr, int nc, unsigned int r, unsigned int c) const {
  const Matrix& m = *this;
  s.clear();
  for(unsigned int i=0; i<nr; ++i) {
    s.push_back(std::vector<double>(nc,0.0));
    for(unsigned int j=0; j<nc; ++j) s[i][j] = m[r+i][c+j];
  }
  return s;
}

Matrix& Matrix::setM(Matrix& s, int nr, int nc, unsigned int r, unsigned int c) {
  Matrix& m = *this;
  for(unsigned int i=0; i<nr; ++i) {
    s.push_back(std::vector<double>(nc,0.0));
    for(unsigned int j=0; j<nc; ++j) m[r+i][c+j] = s[i][j];
  }
  return m;
}

Matrix& Matrix::getMx(Matrix& s) const {
  this->getM(s, 2, 2, 0, 0);
  return s;
}

Matrix& Matrix::getMy(Matrix& s) const {
  this->getM(s, 2, 2, 2, 2);
  return s;
}

Matrix& Matrix::inverse_symplectic(int size, unsigned int r, unsigned int c) {
  Matrix& m = *this;
  if (size < 0) size = m.size();
  if (size == 2) {
    std::swap(m[r+0][c+0],m[r+1][c+1]);
    m[r+0][c+1] = - m[r+0][c+1];
    m[r+1][c+0] = - m[r+1][c+0];
  } else {
    Matrix J({{+0, +1, +0,+0, +0, +0},
              {-1, +0, +0,+0, +0, +0},
              {+0, +0, +0,+1, +0, +0},
              {+0, +0, -1,+0, +0, +0},
              {+0, +0, +0,+0, +0, +1},
              {+0, +0, +0,+0, -1, +0}}
              );
    Matrix m1(size);
    m.transpose();
    m1.multiplication(m,J);
    m.multiplication(J,m1);
    m.scalar(-1);
  }
  return m;
}

static void matrix_inverse6_lu(Matrix& M) {

  gsl_matrix* m = gsl_matrix_alloc(6,6);
  gsl_matrix* inv_m = gsl_matrix_alloc(6,6);
  gsl_permutation* p = gsl_permutation_alloc(6);
  for(unsigned int i=0; i<6; ++i)
    for(unsigned int j=0; j<6; ++j) gsl_matrix_set(m,i,j,M[i][j]);
  int s;
  gsl_linalg_LU_decomp(m, p, &s);
  gsl_linalg_LU_invert(m, p, inv_m);
  for(unsigned int i=0; i<6; ++i)
    for(unsigned int j=0; j<6; ++j) M[i][j] = gsl_matrix_get(inv_m,i,j);
  gsl_matrix_free(m);
  gsl_matrix_free(inv_m);
  gsl_permutation_free(p);

}

Matrix& Matrix::inverse(int size, unsigned int r, unsigned int c) {
  Status::type status = Status::success;
  Matrix& m = *this;
  if (size < 0) size = m.size();
  if (size == 2) {
    double det = m[r+0][c+0] * m[r+1][c+1] - m[r+0][c+1] * m[r+1][c+0];
    std::swap(m[r+0][c+0],m[r+1][c+1]);
    m[r+0][c+1] = - m[r+0][c+1] / det;
    m[r+1][c+0] = - m[r+1][c+0] / det;
    m[r+0][c+0] /= det; m[r+1][c+1] /= det;
  } else if (size == 6) {
    matrix_inverse6_lu(m);
  } else {
    status = Status::not_implemented;
  }
  return m;
}

Status::type matrix_inverse4_blockwise(Matrix& m) {

  Matrix A({{m[0][0],m[0][1]},{m[1][0],m[1][1]}});
  Matrix B({{m[0][2],m[0][3]},{m[1][2],m[1][3]}});
  Matrix C({{m[2][0],m[2][1]},{m[3][0],m[3][1]}});
  Matrix D({{m[2][2],m[2][3]},{m[3][2],m[3][3]}});

  Matrix t1(2);
  Matrix t2(2);
  Matrix t3(2);
  Matrix t4(2);
  Matrix t5(2);

  A.inverse();
  t1.multiplication(C,A);            // t1 = C * A^-1
  t2.multiplication(A,B);            // t2 = A^-1 * B

  t4.multiplication(C,t2);           // t4 = C * A^-1 * B
  t3.linear_combination(1,D,-1,t4);  // t3 = D - C * A^-1 * B
  t3.inverse();               // t3 = (D - C * A^-1 * B)^-1

  m[2][2] = t3[0][0]; m[2][3] = t3[0][1];
  m[3][2] = t3[1][0]; m[3][3] = t3[1][1];

  t5.multiplication(t3,t1);          // t1 = (D - C * A^-1 * B)^-1 * C * A^-1

  m[2][0] = -t5[0][0]; m[2][1] = -t5[0][1];
  m[3][0] = -t5[1][0]; m[3][1] = -t5[1][1];

  std::cerr << "matrix_inverse4_blockwise not finished" << std::endl;


}

Vector operator+(const Vector& v1, const Vector& v2) {
  Vector v(v1);
  for(unsigned int i=0; i<v.size(); ++i) v[i] += v2[i];
  return v;
}

Pos<double> linalg_solve4_posvec(const std::vector<Pos<double> >& M, const Pos<double>& B) {

  gsl_matrix* m = gsl_matrix_alloc(4,4);
  gsl_vector* b = gsl_vector_alloc(4);
  gsl_vector* x = gsl_vector_alloc(4);
  gsl_permutation* p = gsl_permutation_alloc(4);

  gsl_vector_set(b,0,B.rx); gsl_vector_set(b,1,B.px);
  gsl_vector_set(b,2,B.ry); gsl_vector_set(b,3,B.py);
  for(unsigned int i=0; i<4; ++i) {
    gsl_matrix_set(m,0,i,M[i].rx); gsl_matrix_set(m,1,i,M[i].px);
    gsl_matrix_set(m,2,i,M[i].ry); gsl_matrix_set(m,3,i,M[i].py);
  }

  int s; gsl_linalg_LU_decomp(m, p, &s);
  gsl_linalg_LU_solve(m, p, b, x);
  Pos<double> X(gsl_vector_get(x,0),gsl_vector_get(x,1),gsl_vector_get(x,2),gsl_vector_get(x,3),0,0);

  gsl_matrix_free(m);
  gsl_vector_free(b);
  gsl_vector_free(x);
  gsl_permutation_free(p);
  return X;
}

Pos<double> linalg_solve6_posvec(const std::vector<Pos<double> >& M, const Pos<double>& B) {

  gsl_matrix* m = gsl_matrix_alloc(6,6);
  gsl_vector* b = gsl_vector_alloc(6);
  gsl_vector* x = gsl_vector_alloc(6);
  gsl_permutation* p = gsl_permutation_alloc(6);

  gsl_vector_set(b,0,B.rx); gsl_vector_set(b,1,B.px);
  gsl_vector_set(b,2,B.ry); gsl_vector_set(b,3,B.py);
  gsl_vector_set(b,4,B.de); gsl_vector_set(b,5,B.dl);
  for(unsigned int i=0; i<6; ++i) {
    gsl_matrix_set(m,0,i,M[i].rx); gsl_matrix_set(m,1,i,M[i].px);
    gsl_matrix_set(m,2,i,M[i].ry); gsl_matrix_set(m,3,i,M[i].py);
    gsl_matrix_set(m,4,i,M[i].de); gsl_matrix_set(m,5,i,M[i].dl);
  }

  int s; gsl_linalg_LU_decomp(m, p, &s);
  gsl_linalg_LU_solve(m, p, b, x);
  Pos<double> X(gsl_vector_get(x,0),gsl_vector_get(x,1),gsl_vector_get(x,2),gsl_vector_get(x,3),gsl_vector_get(x,4),gsl_vector_get(x,5));

  gsl_matrix_free(m);
  gsl_vector_free(b);
  gsl_vector_free(x);
  gsl_permutation_free(p);
  return X;
}
