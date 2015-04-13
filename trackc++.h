#ifndef _TRACKCPP_H
#define _TRACKCPP_H

// TRACKC++
// ========
// Author: 		Ximenes R. Resende
// email:  		xresende@gmail.com, ximenes.resende@lnls.br
// affiliation:	LNLS - Laboratorio Nacional de Luz Sincrotron
// Date: 		Tue Dec 10 17:57:20 BRST 2013


#include "output.h"
#include "commands.h"
#include "optics.h"
#include "dynap.h"
#include "tracking.h"
#include "lattice.h"
#include "flat_file.h"
#include "kicktable.h"
#include "passmethods.h"
#include "accelerator.h"
#include "elements.h"
#include "pos.h"
#include "tpsa.h"
#include "auxiliary.h"

extern bool verbose_on;

void sirius_v500(std::vector<Element>& the_ring);

#endif
