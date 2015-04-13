#ifndef FMAP_MP_H
#define FMAP_MP_H

// TRACKING_MP
// ===========
// Author: 		Ximenes R. Resende
// email:  		xresende@gmail.com, ximenes.resende@lnls.br
// affiliation:	LNLS - Laboratorio Nacional de Luz Sincrotron
// Date: 		2012-08-20


#include <stdlib.h>

//#include "../../tracy_lib.h"
//#include "../../soleilcommon.h"
//#include "../../naffutils.h"
//extern globvalrec globval;

#include "mp_task_mgr.h"

void fmap_mp(int nr_cpus, long Nbx, long Nbz, long Nbtour, double x0, double xmax,
		  double z0, double zmax, double energy, bool diffusion);


#endif

