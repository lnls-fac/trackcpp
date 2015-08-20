#ifndef _MULTITHREADS_H
#define _MULTITHREADS_H

#include <pthread.h>

struct ThreadSharedData {
  long task_id;
	long nr_tasks;
	void   (*func)(ThreadSharedData*, int, long);
	pthread_mutex_t *mutex;
};

void start_all_threads(ThreadSharedData& thread_data, unsigned int nr_threads);

#endif
