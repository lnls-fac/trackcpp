#ifndef MP_TASK_MGR_H
#define MP_TASK_MGR_H

// TRACKING_MP
// ===========
// Author: 		Ximenes R. Resende
// email:  		xresende@gmail.com, ximenes.resende@lnls.br
// affiliation:	LNLS - Laboratorio Nacional de Luz Sincrotron
// Date: 		2012-08-20

#include <string.h>
#include <cmath>
#include <semaphore.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/shm.h>
#include <sys/stat.h>
#include <sys/wait.h>

void mp_task_mgr(int nr_cpus, int nr_pnts, void *data_in, void *data_out, int data_out_size, void (*parallelized_function)(int cpu_id, int data_id, void *data_in, void *shm_data_out));

#endif

