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

class Vector : public std::vector<double> {
public:
    Vector(const unsigned int size = 0) : std::vector<double>(size, 0) {}
    Vector(const std::vector<double>& v) : std::vector<double>(v) {}
    Vector& multiplication(const Matrix& m, const Vector& b);
};

class Matrix : public std::vector<std::vector<double> > {
public:
  Matrix(const unsigned int size = 0) : std::vector<std::vector<double> >(size, std::vector<double>(size,0)) {}
  Matrix(const std::vector<std::vector<double> >& v) : std::vector<std::vector<double> >(v) {}
  double norm() const;
  Matrix& eye(const double& v = 1);
  Matrix& scalar(const double& v);
  Matrix& transpose(int size = -1, unsigned int r=0, unsigned int c=0);
  Matrix& linear_combination(const double& a1, const Matrix& m1, const double& a2, const Matrix& m2);
  Matrix& multiplication(const Matrix& m1, const Matrix& m2);
  Matrix& getM(Matrix& s, int nr, int nc, unsigned int r=0, unsigned int c=0) const;
  Matrix& setM(Matrix& s, int nr, int nc, unsigned int r=0, unsigned int c=0);
  Matrix& getMx(Matrix& s) const;
  Matrix& getMy(Matrix& s) const;
  Matrix& inverse_symplectic(int size=-1, unsigned int r=0, unsigned int c=0);
  Matrix& inverse(int size=-1, unsigned int r=0, unsigned int c=0);
};

namespace std {
    %template(CppMatrixVector) vector<Matrix>;
}
