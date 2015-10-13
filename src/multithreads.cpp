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

#include <trackcpp/trackcpp.h>

static int current_thread_id = 0;
static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

void* start_thread(void* args) {

  // gets thread_id from global variable
  pthread_mutex_lock(&mutex);
  int this_thread_id = current_thread_id++;
  pthread_mutex_unlock(&mutex);

  // gets pointer to shared input data
  ThreadSharedData* data = (ThreadSharedData*) args;

  while (true) {

    pthread_mutex_lock(&mutex);
    long this_task_id = data->task_id++;
    pthread_mutex_unlock(&mutex);

    // breaks if there is no more task to be done.
    if (this_task_id >= data->nr_tasks) break;

    // run main function
    data->func(data, this_thread_id, this_task_id);

  }

  return NULL;
}

void start_all_threads(ThreadSharedData& thread_data, unsigned int nr_threads) {

  thread_data.task_id = 0;
  thread_data.mutex = &mutex;

  pthread_t threads[nr_threads];
  int thread_id[nr_threads];

  for(int i=0; i<nr_threads-1; i++) {
  	thread_id[i] = i;
  	int iret = pthread_create(&(threads[i]), NULL, start_thread, (void*) &thread_data);
  }

  thread_id[nr_threads-1] = current_thread_id;
  start_thread((void*) &thread_data);

  for(int i=0; i<nr_threads-1; ++i) pthread_join(threads[i], NULL);

}
