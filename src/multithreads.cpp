// The MIT License (MIT)
//
// Copyright (c) 2015 LNLS Accelerator Division
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

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
