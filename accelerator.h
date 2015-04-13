#ifndef _ACCELERATOR_H
#define _ACCELERATOR_H

#include "kicktable.h"
#include "elements.h"
#include <vector>

struct Accelerator {

  double                  energy;              // [eV]
  bool                    cavity_on;
  bool                    radiation_on;
  bool                    vchamber_on;
  int                     harmonic_number;
  std::vector<Element>    lattice;
  std::vector<Kicktable*> kicktables;

};

#endif
