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
