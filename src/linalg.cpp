#include <limits>
#include <trackcpp/linalg.h>


double matrix_norm(const Matrix& m) {

  double max = 0;
  { const int i=0, j=0; max = std::max(max,abs(m[i][j])); }
  { const int i=0, j=1; max = std::max(max,abs(m[i][j])); }
  { const int i=0, j=2; max = std::max(max,abs(m[i][j])); }
  { const int i=0, j=3; max = std::max(max,abs(m[i][j])); }
  { const int i=0, j=4; max = std::max(max,abs(m[i][j])); }
  { const int i=0, j=5; max = std::max(max,abs(m[i][j])); }

  { const int i=1, j=0; max = std::max(max,abs(m[i][j])); }
  { const int i=1, j=1; max = std::max(max,abs(m[i][j])); }
  { const int i=1, j=2; max = std::max(max,abs(m[i][j])); }
  { const int i=1, j=3; max = std::max(max,abs(m[i][j])); }
  { const int i=1, j=4; max = std::max(max,abs(m[i][j])); }
  { const int i=1, j=5; max = std::max(max,abs(m[i][j])); }

  { const int i=2, j=0; max = std::max(max,abs(m[i][j])); }
  { const int i=2, j=1; max = std::max(max,abs(m[i][j])); }
  { const int i=2, j=2; max = std::max(max,abs(m[i][j])); }
  { const int i=2, j=3; max = std::max(max,abs(m[i][j])); }
  { const int i=2, j=4; max = std::max(max,abs(m[i][j])); }
  { const int i=2, j=5; max = std::max(max,abs(m[i][j])); }

  { const int i=3, j=0; max = std::max(max,abs(m[i][j])); }
  { const int i=3, j=1; max = std::max(max,abs(m[i][j])); }
  { const int i=3, j=2; max = std::max(max,abs(m[i][j])); }
  { const int i=3, j=3; max = std::max(max,abs(m[i][j])); }
  { const int i=3, j=4; max = std::max(max,abs(m[i][j])); }
  { const int i=3, j=5; max = std::max(max,abs(m[i][j])); }

  { const int i=4, j=0; max = std::max(max,abs(m[i][j])); }
  { const int i=4, j=1; max = std::max(max,abs(m[i][j])); }
  { const int i=4, j=2; max = std::max(max,abs(m[i][j])); }
  { const int i=4, j=3; max = std::max(max,abs(m[i][j])); }
  { const int i=4, j=4; max = std::max(max,abs(m[i][j])); }
  { const int i=4, j=5; max = std::max(max,abs(m[i][j])); }

  { const int i=5, j=0; max = std::max(max,abs(m[i][j])); }
  { const int i=5, j=1; max = std::max(max,abs(m[i][j])); }
  { const int i=5, j=2; max = std::max(max,abs(m[i][j])); }
  { const int i=5, j=3; max = std::max(max,abs(m[i][j])); }
  { const int i=5, j=4; max = std::max(max,abs(m[i][j])); }
  { const int i=5, j=5; max = std::max(max,abs(m[i][j])); }

  return max;
}

void matrix_eye(Matrix& m, const double& v) {
  for(unsigned int i=0; i<6; ++i) {
      for(unsigned int j=0; j<6; ++j) m[i][j] = 0.0;
      m[i][i] = v;
    }
}

void matrix_print(const Matrix& m) {
  for(unsigned int i=0; i<6; ++i) {
    for(unsigned int j=0; j<6; ++j) {
      printf("%+.4e ", m[i][j]);
    }
    printf("\n");
  }
}

void matrix_transpose(Matrix& m) {
  for(unsigned int i=0; i<6; ++i) {
    for(unsigned int j=i+1; j<6; ++j) {
        double t = m[i][j]; m[i][j] = m[j][i]; m[j][i] = t;
    }
  }
}

void matrix_scalar(Matrix& m, double scalar) {

  { const int i=0, j=0; m[i][j] *= scalar; }
  { const int i=0, j=1; m[i][j] *= scalar; }
  { const int i=0, j=2; m[i][j] *= scalar; }
  { const int i=0, j=3; m[i][j] *= scalar; }
  { const int i=0, j=4; m[i][j] *= scalar; }
  { const int i=0, j=5; m[i][j] *= scalar; }

  { const int i=1, j=0; m[i][j] *= scalar; }
  { const int i=1, j=1; m[i][j] *= scalar; }
  { const int i=1, j=2; m[i][j] *= scalar; }
  { const int i=1, j=3; m[i][j] *= scalar; }
  { const int i=1, j=4; m[i][j] *= scalar; }
  { const int i=1, j=5; m[i][j] *= scalar; }

  { const int i=2, j=0; m[i][j] *= scalar; }
  { const int i=2, j=1; m[i][j] *= scalar; }
  { const int i=2, j=2; m[i][j] *= scalar; }
  { const int i=2, j=3; m[i][j] *= scalar; }
  { const int i=2, j=4; m[i][j] *= scalar; }
  { const int i=2, j=5; m[i][j] *= scalar; }

  { const int i=3, j=0; m[i][j] *= scalar; }
  { const int i=3, j=1; m[i][j] *= scalar; }
  { const int i=3, j=2; m[i][j] *= scalar; }
  { const int i=3, j=3; m[i][j] *= scalar; }
  { const int i=3, j=4; m[i][j] *= scalar; }
  { const int i=3, j=5; m[i][j] *= scalar; }

  { const int i=4, j=0; m[i][j] *= scalar; }
  { const int i=4, j=1; m[i][j] *= scalar; }
  { const int i=4, j=2; m[i][j] *= scalar; }
  { const int i=4, j=3; m[i][j] *= scalar; }
  { const int i=4, j=4; m[i][j] *= scalar; }
  { const int i=4, j=5; m[i][j] *= scalar; }

  { const int i=5, j=0; m[i][j] *= scalar; }
  { const int i=5, j=1; m[i][j] *= scalar; }
  { const int i=5, j=2; m[i][j] *= scalar; }
  { const int i=5, j=3; m[i][j] *= scalar; }
  { const int i=5, j=4; m[i][j] *= scalar; }
  { const int i=5, j=5; m[i][j] *= scalar; }

}

void matrix_linear_combination2(Matrix& m, const double& a1, const Matrix& m1, const double& a2, const Matrix& m2) {

  { const int i=0, j=0; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=0, j=1; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=1, j=0; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=1, j=1; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }

}

void matrix_linear_combination(Matrix& m, const double& a1, const Matrix& m1, const double& a2, const Matrix& m2) {

  { const int i=0, j=0; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=0, j=1; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=0, j=2; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=0, j=3; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=0, j=4; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=0, j=5; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }

  { const int i=1, j=0; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=1, j=1; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=1, j=2; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=1, j=3; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=1, j=4; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=1, j=5; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }

  { const int i=2, j=0; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=2, j=1; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=2, j=2; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=2, j=3; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=2, j=4; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=2, j=5; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }

  { const int i=3, j=0; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=3, j=1; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=3, j=2; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=3, j=3; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=3, j=4; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=3, j=5; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }

  { const int i=4, j=0; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=4, j=1; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=4, j=2; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=4, j=3; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=4, j=4; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=4, j=5; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }

  { const int i=5, j=0; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=5, j=1; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=5, j=2; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=5, j=3; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=5, j=4; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }
  { const int i=5, j=5; m[i][j] = a1 * m1[i][j] + a2 * m2[i][j]; }


}

void matrix_multiplication2(Matrix& m, const Matrix& m1, const Matrix& m2) {
  m[0][0] = m1[0][0] * m2[0][0] + m1[0][1] * m2[1][0];
  m[0][1] = m1[0][0] * m2[0][1] + m1[0][1] * m2[1][1];
  m[1][0] = m1[1][0] * m2[0][0] + m1[1][1] * m2[1][0];
  m[1][1] = m1[1][0] * m2[0][1] + m1[1][1] * m2[1][1];
}

void matrix_multiplication(Matrix& m, const Matrix& m1, const Matrix& m2) {

  { const int i=0, j=0; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=0, j=1; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=0, j=2; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=0, j=3; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=0, j=4; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=0, j=5; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }

  { const int i=1, j=0; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=1, j=1; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=1, j=2; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=1, j=3; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=1, j=4; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=1, j=5; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }

  { const int i=2, j=0; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=2, j=1; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=2, j=2; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=2, j=3; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=2, j=4; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=2, j=5; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }

  { const int i=3, j=0; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=3, j=1; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=3, j=2; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=3, j=3; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=3, j=4; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=3, j=5; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }

  { const int i=4, j=0; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=4, j=1; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=4, j=2; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=4, j=3; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=4, j=4; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=4, j=5; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }

  { const int i=5, j=0; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=5, j=1; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=5, j=2; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=5, j=3; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=5, j=4; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }
  { const int i=5, j=5; m[i][j] = m1[i][0]*m2[0][j]+m1[i][1]*m2[1][j]+m1[i][2]*m2[2][j]+m1[i][3]*m2[3][j]+m1[i][4]*m2[4][j]+m1[i][5]*m2[5][j]; }

}

void matrix_symplectic_inverse2(Matrix& m, unsigned int i, unsigned int j) {
  //std::swap(m[i+0][j+0],m[i+1][j+1]);
  double t = m[i+0][j+0]; m[i+0][j+0] = m[i+1][j+1]; m[i+1][j+1] = t;
  m[i+0][j+1] = - m[i+0][j+1];
  m[i+1][j+0] = - m[i+1][j+0];
}

void matrix_inverse2(Matrix& m, unsigned int i, unsigned int j) {
  double det = m[i+0][j+0]*m[i+1][j+1] - m[i+0][j+1]*m[i+1][j+0];
  std::swap(m[i+0][j+0],m[i+1][j+1]);
  m[i+0][j+0] /= det;
  m[i+1][j+1] /= det;
  m[i+0][j+1] = - m[i+0][j+1]/det;
  m[i+1][j+0] = - m[i+1][j+0]/det;
}

void matrix_symplectic_inverse6(Matrix& m) {
  Matrix J = {{+0,+1, +0,+0, +0,+0},
              {-1,+0, +0,+0, +0,+0},
              {+0,+0, +0,+1, +0,+0},
              {+0,+0, -1,+0, +0,+0},
              {+0,+0, +0,+0, +0,+1},
              {+0,+0, +0,+0, -1,+0},
            };

  Matrix m1 = m;
  matrix_transpose(m);
  matrix_multiplication(m1,m,J);
  matrix_multiplication(m,J,m1);
  matrix_scalar(m,-1);

  //std::cout << "mt:" << std::endl;
  //matrix_print(m);
  // matrix_symplectic_inverse2(m, 0, 0); // -J M_xx J
  // matrix_symplectic_inverse2(m, 0, 2); // -J M_xy J
  // matrix_symplectic_inverse2(m, 0, 4); // -J M_xz J
  // matrix_symplectic_inverse2(m, 2, 0); // -J M_yx J
  // matrix_symplectic_inverse2(m, 2, 2); // -J M_yy J
  // matrix_symplectic_inverse2(m, 2, 4); // -J M_yz J
  // matrix_symplectic_inverse2(m, 4, 0); // -J M_zx J
  // matrix_symplectic_inverse2(m, 4, 2); // -J M_zy J
  // matrix_symplectic_inverse2(m, 4, 4); // -J M_zz J
  // matrix_transpose(m);
}

#include <chrono>
Status::type matrix_inverse6_newton(Matrix& m) {

  auto start = std::chrono::steady_clock::now();
  auto end = std::chrono::steady_clock::now();
  auto diff = end - start;

  start = std::chrono::steady_clock::now();
  // Newton algorithm
  Matrix x0(6,std::vector<double>(6,0));
  Matrix x1(6,std::vector<double>(6,0));
  Matrix t1(6,std::vector<double>(6,0));
  Matrix t2(6,std::vector<double>(6,0));
  Matrix dm(6,std::vector<double>(6,0));
  x0 = m;
  matrix_symplectic_inverse6(x0);      //  x0 = -j m j
  unsigned int iter_left = 20;
  while (true) {
    //std::cout << std::endl << "new iteration" << std::endl;
    //std::cout << "x0:" << std::endl;
    //matrix_print(x0);
    matrix_multiplication(t1,m,x0);            //  t1 = m x0
    //std::cout << "t1:" << std::endl;
    //matrix_print(t1);
    matrix_multiplication(t2,x0,t1);           //  t2 = x0 t1 = x0 m x0
    //std::cout << "t2:" << std::endl;
    //matrix_print(t2);
    matrix_linear_combination(x1, 2,x0,-1,t2);  //  x1 = 2 x0 - x0 m x0
    //std::cout << "x1:" << std::endl;
    //matrix_print(m);
    matrix_linear_combination(dm,1,x1,-1,x0);  //  dm = x1 - x0
    double norm = matrix_norm(dm);
    //std::cout << "dm:" << std::endl;
    //matrix_print(dm);
    //std::cout << "norm: " << norm << std::endl;
    //std::cout << "epsilon: " <<  std::numeric_limits<double>::epsilon() << std::endl;
    if (norm < std::numeric_limits<double>::epsilon()) break;
    x0 = x1;
    iter_left--;
    if (iter_left == 0) return Status::newton_not_converged;
  }
  m = x1;

  end = std::chrono::steady_clock::now();
  diff = end - start;
  std::cout << "newton: " << std::chrono::duration <double, std::milli> (diff).count() << " ms" << " in " << 20-iter_left << " iterations" << std::endl;

  return Status::success;

}


// Status::type matrix_inverse4_blockwise(Matrix& m) {
//   Matrix t1;
//   Matrix A = {{m[0][0],m[0][1]},{m[1][0],m[1][1]}};
//   Matrix B = {{m[0][2],m[0][3]},{m[1][2],m[1][3]}};
//   Matrix C = {{m[2][0],m[2][1]},{m[3][0],m[3][1]}};
//   Matrix D = {{m[2][2],m[2][3]},{m[3][2],m[3][3]}};
//   matrix_inverse2(A);
//   //matrix_multiplication2(t1,C,A);
//   //matrix_multiplication2(t2,t1,B);
//   //matrix_linear_combination2(t1,1,D,-1,t2);
// }



Status::type matrix_inverse(Matrix& m, const unsigned int size, const unsigned int i, const unsigned int j) {


  switch (size) {
    case 1: m[i+0][j+0] = 1.0 / m[i+0][j+0]; break;
    case 2: matrix_inverse2(m, i, j); break;
    case 6: matrix_inverse6_newton(m); break;

  }

}

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
