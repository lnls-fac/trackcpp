#ifndef FMAPDP_MP_H
#define FMAPDP_MP_H

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

void fmapdp_mp(int nr_cpus, long Nbx, long Nbe, long Nbtour, double x0, double xmax,
		double emin, double emax, double z, bool diffusion);


#endif

