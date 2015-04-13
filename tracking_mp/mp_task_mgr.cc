// TRACKING_MP
// ===========
// Author: 		Ximenes R. Resende
// email:  		xresende@gmail.com, ximenes.resende@lnls.br
// affiliation:	LNLS - Laboratorio Nacional de Luz Sincrotron
// Date: 		2012-08-20

#include "mp_task_mgr.h"


#define FMAP_MP_DATA_RAW 0
#define FMAP_MP_DATA_PRO 1
#define FMAP_MP_DATA_RDY 2
#define SHMKEY1 0
#define SHMKEY2 1

static int  get_next_data_idx(int *data_flag, int nr_pnts);


void mp_task_mgr(int nr_cpus, int nr_pnts, void *data_in, void *data_out, int data_out_size, void (*parallelized_function)(int cpu_id, int data_id, void *data_in, void *shm_data_out)) {

	int i,data_idx,child_status;

	int 	cpu_id;
	sem_t 	mutex;
	//pid_t 	pid;
	int     shmid1, shmid2;
	int    	*shm_data_flag; // 0 : dado ñ processado, 1 : em processamento, 2 : já processado
	void    *shm_data_out;
	pid_t  	*children;

	printf("<mp_task_mgr> START\n");

	children = (pid_t*) malloc ((nr_cpus-1) * sizeof(int));


	// cria segmentos de memória compartilhada (apontado por 'data_flag') para bookkeeping dos dados já processados
	shmid1 = shmget(getpid() + SHMKEY1, sizeof(int) * nr_pnts, S_IRWXU | S_IRWXG | S_IRWXO | IPC_CREAT);		
	if (shmid1 == -1) { perror("shmget"); exit(1); }
	shm_data_flag = (int*) shmat(shmid1, (void*)0, 0);
	for(i=0; i<nr_pnts; ++i) {
		shm_data_flag[i] = FMAP_MP_DATA_RAW;
	}
	shmid2 = shmget(getpid() + SHMKEY2, data_out_size, S_IRWXU | S_IRWXG | S_IRWXO | IPC_CREAT);		
	if (shmid2 == -1) { perror("shmget"); exit(1); }
	shm_data_out = (void*) shmat(shmid2, (void*)0, 0);
	bzero(shm_data_out, data_out_size);


	// cria mutex para acesso à memória compartilhada
	sem_init(&mutex, 1, 1);

	// gerar processos
	cpu_id = 0;
	for (i=0; i<nr_cpus-1; ++i) {
		children[i] = fork();
		if (children[i] == 0) { 
			cpu_id = i+1;
			break; 
		};
	}

	printf("<mp_task_mgr> cpu %i running...\n", cpu_id+1);

	// loop principal de execução (todos processos)
	while(1) {
		sem_wait(&mutex);
		data_idx = get_next_data_idx(shm_data_flag, nr_pnts);
		if (data_idx == nr_pnts) {
			sem_post(&mutex); // libera acesso a outros processos
			free(children);
			if (cpu_id == 0) {
				break;
			} else {
				printf("<mp_task_mgr> cpu %i exiting...\n", cpu_id+1);
				_exit(1);
			}
		}
		shm_data_flag[data_idx] = 1; // em calculo...
		sem_post(&mutex);

		// calculo
		//printf("<mp_task_mgr> cpu %i processing data %i/%i...\n", cpu_id+1, data_idx+1, nr_pnts);

		// run_tracking(cpu_id, data_idx, data_in, shm_data_out);
		parallelized_function(cpu_id, data_idx, data_in, shm_data_out);

		//double x1 = ((double*)shm_data_out)[data_idx * 2 + 0];
		//double z1 = ((double*)shm_data_out)[data_idx * 2 + 1]; 
		//printf("%f %f\n", x1, z1);
		sem_wait(&mutex);
		shm_data_flag[data_idx] = 2; // calculado.
		sem_post(&mutex);
	}

	
	for(i=0; i<nr_cpus-1; ++i) {
		waitpid(children[i], &child_status, 0);
	}

	
	memcpy(data_out, shm_data_out, data_out_size);

	
	// deleta mutex
	sem_destroy(&mutex);
	// apaga segmentos de memoria compartilhada criados
	if (shmdt(shm_data_flag) == -1) { perror("shmdt"); exit(1); }
	shmctl(shmid1, IPC_RMID, NULL);
	if (shmdt(shm_data_out) == -1) { perror("shmdt"); exit(1); }
	shmctl(shmid2, IPC_RMID, NULL);

	printf("<mp_task_mgr> FINISH\n");

}

#undef FMAP_MP_DATA_RAW
#undef FMAP_MP_DATA_PRO
#undef FMAP_MP_DATA_RDY
#undef SHMKEY1
#undef SHMKEY2

int get_next_data_idx(int *data_flag, int nr_pnts) {

	int i;
	for(i=0; (i<nr_pnts) && (data_flag[i] != 0); ++i);
	return i;

}

