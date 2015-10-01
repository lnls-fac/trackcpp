#include <trackcpp/linalg.h>

void getmx(const Matrix& m, Matrix2& mx) {
  mx[0][0] = m[0][0]; mx[0][1] = m[0][1];
  mx[1][0] = m[1][0]; mx[1][1] = m[1][1];
}

void getmy(const Matrix& m, Matrix2& my) {
  Matrix2 r;
  my[0][0] = m[2][2]; my[0][1] = m[2][3];
  my[1][0] = m[3][2]; my[1][1] = m[3][3];
}

void matrix2_eye(Matrix2& m) {
  m[0][0] = 1; m[0][1] = 0;
  m[1][0] = 0; m[1][1] = 1;
}
void matrix2_lc(Matrix2& m, const double& a1, const Matrix2& m1, const double& a2, const Matrix2& m2) {
  m[0][0] = a1 * m1[0][0] + a2 * m2[0][0]; m[0][1] = a1 * m1[0][1] + a2 * m2[0][1];
  m[1][0] = a1 * m1[1][0] + a2 * m2[1][0]; m[1][1] = a1 * m1[1][1] + a2 * m2[1][1];
}

void linalg_solve2(Vector2& X, const Matrix2& M, const Vector2& B) {
  gsl_matrix* m = gsl_matrix_alloc(2,2);
  gsl_vector* b = gsl_vector_alloc(2);
  gsl_vector* x = gsl_vector_alloc(2);
  gsl_permutation* p = gsl_permutation_alloc(2);
  gsl_vector_set(b,0,B[0]); gsl_vector_set(b,1,B[1]);
  for(unsigned int i=0; i<2; ++i) {
    gsl_matrix_set(m,0,i,M[i][0]); gsl_matrix_set(m,1,i,M[i][1]);
  }
  int s; gsl_linalg_LU_decomp(m, p, &s);
  gsl_linalg_LU_solve(m, p, b, x);
  X[0] = gsl_vector_get(x,0);
  X[1] = gsl_vector_get(x,1);
  gsl_matrix_free(m);
  gsl_vector_free(b);
  gsl_vector_free(x);
  gsl_permutation_free(p);
}

Pos<double> linalg_solve2_posvec(const std::vector<Pos<double> >& M, const Pos<double>& B) {

  gsl_matrix* m = gsl_matrix_alloc(2,2);
  gsl_vector* b = gsl_vector_alloc(2);
  gsl_vector* x = gsl_vector_alloc(2);
  gsl_permutation* p = gsl_permutation_alloc(2);

  gsl_vector_set(b,0,B.rx); gsl_vector_set(b,1,B.px);
  for(unsigned int i=0; i<2; ++i) {
    gsl_matrix_set(m,0,i,M[i].rx); gsl_matrix_set(m,1,i,M[i].px);
  }

  int s; gsl_linalg_LU_decomp(m, p, &s);
  gsl_linalg_LU_solve(m, p, b, x);
  Pos<double> X(gsl_vector_get(x,0),gsl_vector_get(x,1),0,0,0,0);

  gsl_matrix_free(m);
  gsl_vector_free(b);
  gsl_vector_free(x);
  gsl_permutation_free(p);
  return X;
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
