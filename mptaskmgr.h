#ifndef _MP_TASK_MGR_H
#define _MP_TASK_MGR_H

// MP_TASK_MGR
// ===========
// Author: 		Ximenes R. Resende
// email:  		xresende@gmail.com, ximenes.resende@lnls.br
// affiliation:	LNLS - Laboratorio Nacional de Luz Sincrotron
// Date: 		2012-08-20

#define MP_DATA_RAW 0
#define MP_DATA_PRO 1
#define MP_DATA_RDY 2
#define SHMKEY1 0
#define SHMKEY2 1

void mp_task_mgr(int nr_cpus, int nr_pnts, void *data_in, void *data_out, int data_out_size, void (*parallelized_function)(int cpu_id, int data_id, void *data_in, void *shm_data_out));

#endif
