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

#include <trackcpp/auxiliary.h>

static std::random_device rand_dev;
static std::mt19937 generator(rand_dev());
static std::normal_distribution<double>  distr_gauss(0., 1.);
static std::uniform_real_distribution<double> distr_uniform(-sqrt(3.0), sqrt(3.0));
int choosen_distribution = Distributions::normal;

void set_random_distribution(unsigned value){
  if (value == Distributions::normal){
    choosen_distribution = value;
  }else if (value == Distributions::uniform){
    choosen_distribution = value;
  }else{
    std::cerr << "Not valid input!" << std::endl;
  }
}

void set_random_seed(unsigned rnd_seed) {
    distr_gauss.reset();
    generator.seed(rnd_seed);
}

void set_random_seed_with_random_device() {
    distr_gauss.reset();
    generator.seed(rand_dev());
}

double gen_random_number() {
  if (choosen_distribution==Distributions::uniform){
    return distr_uniform(generator);
  }else{
    return distr_gauss(generator);
  }
  // Box-Muller Transform
  //
  // z1 = sqrt(-2*ln(u1)) * cos(2*pi*u2)
  // z2 = sqrt(-2*ln(u1)) * sin(2*pi*u2)
  //
  // u1, u2: uniform [0, 1]
  // z1, z2: normal dist.
}
