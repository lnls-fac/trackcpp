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
static HANDLE mutex;

DWORD WINAPI start_thread(void* args) {

  // gets thread_id from global variable
  int this_thread_id;
  WaitForSingleObject(mutex, INFINITE);
  this_thread_id = current_thread_id++;
  ReleaseMutex(mutex);

  // gets pointer to shared input data
  ThreadSharedData* data = (ThreadSharedData*) args;

  while (true) {

    WaitForSingleObject(mutex, INFINITE);
    long this_task_id = data->task_id++;
    ReleaseMutex(mutex);

    // breaks if there is no more task to be done.
    if (this_task_id >= data->nr_tasks) break;

    // run main function
    data->func(data, this_thread_id, this_task_id);

  }

  return 0;
}

void start_all_threads(ThreadSharedData& thread_data, unsigned int nr_threads) {
  // Inicializamos o mutex
  mutex = CreateMutex(NULL, FALSE, NULL);

  // Inicializamos a ID da tarefa
  thread_data.task_id = 0;

  // Criamos as threads
  HANDLE* threads = new HANDLE[nr_threads];
  for (unsigned int i = 0; i < nr_threads; ++i) {
      threads[i] = CreateThread(NULL, 0, start_thread, &thread_data, 0, NULL);
  }

  // Esperamos pelo término das threads
  WaitForMultipleObjects(nr_threads, threads, TRUE, INFINITE);

  // Liberamos os recursos
  CloseHandle(mutex);
  delete[] threads;
}
