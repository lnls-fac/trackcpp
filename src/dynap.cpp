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
#include <trackcpp/output.h>
#include <trackcpp/dynap.h>
#include <trackcpp/tracking.h>
#include <trackcpp/lattice.h>
#include <trackcpp/pos.h>
#include <trackcpp/auxiliary.h>
#include <algorithm>
#include <vector>
#include <cfloat>

extern void naff_traj(const std::vector<Pos<double>>& data, double& tunex, double& tuney);
static const double tiny_y_amp = 1e-7; // [m]


// declaration of auxiliary functions
static Status::type   calc_closed_orbit(const Accelerator& accelerator, std::vector<Pos<double> >& cod, const char* function_name);
//static DynApGridPoint find_momentum_acceptance(const Accelerator& accelerator, const std::vector<Pos<double> >& cod, unsigned int nr_turns, const Pos<double>& p0, double e0, double e_tol, unsigned int element_idx);
//static DynApGridPoint find_fine_momentum_acceptance(const Accelerator& accelerator, const std::vector<Pos<double> >& cod, unsigned int nr_turns, const Pos<double>& p0, double e_init, double e_tol, unsigned int element_idx);
//static DynApGridPoint find_px_acceptance(const Accelerator& accelerator, const std::vector<Pos<double> >& cod, unsigned int nr_turns, const Pos<double>& p0, double px0, double px_tol, unsigned int element_idx);
//static DynApGridPoint find_fine_px_acceptance(const Accelerator& accelerator, const std::vector<Pos<double> >& cod, unsigned int nr_turns, const Pos<double>& p0, double px_init, double px_tol, unsigned int element_idx);
//static DynApGridPoint find_py_acceptance(const Accelerator& accelerator, const std::vector<Pos<double> >& cod, unsigned int nr_turns, const Pos<double>& p0, double py0, double py_tol, unsigned int element_idx);
//static DynApGridPoint find_fine_py_acceptance(const Accelerator& accelerator, const std::vector<Pos<double> >& cod, unsigned int nr_turns, const Pos<double>& p0, double py_init, double py_tol, unsigned int element_idx);

// global thread variables and decl. of aux. thread functions
static unsigned int                     thread_nr_turns = 0;
static std::string                      thread_type = "";
static const Accelerator*               thread_accelerator = NULL;
static const std::vector<Pos<double>>*  thread_cod = NULL;
static std::vector<DynApGridPoint>*     thread_grid = NULL;
static const std::vector<unsigned int>* thread_elements = NULL;
static const double*                    thread_ma_e0    = NULL;
static const double*                    thread_ma_e_tol = NULL;
static const double*                    thread_ma_e_init  = NULL;
static const double*                    thread_ma_e_delta = NULL;
static const unsigned int*              thread_ma_nr_steps_back = NULL;
static const double*                    thread_ma_rescale = NULL;
static const unsigned int*              thread_ma_nr_iterations = NULL;
static const Pos<double>*               thread_ma_p0 = NULL;

static void           thread_dynap_naff(ThreadSharedData* thread_data, int thread_id, long task_id);
static void           thread_dynap_acceptance(ThreadSharedData* thread_data, int thread_id, long task_id);
static void           thread_dynap(ThreadSharedData* thread_data, int thread_id, long task_id);
//static void           thread_dynap_ma(ThreadSharedData* thread_data, int thread_id, long task_id);
//static void           thread_dynap_pxa(ThreadSharedData* thread_data, int thread_id, long task_id);
//static void           thread_dynap_pya(ThreadSharedData* thread_data, int thread_id, long task_id);

// main functions

Status::type dynap_xy(
    const Accelerator& accelerator,
    std::vector<Pos<double> >& cod,
    unsigned int nr_turns,
    const Pos<double>& p0,
    unsigned int nrpts_x, double x_min, double x_max,
    unsigned int nrpts_y, double y_min, double y_max,
    bool calculate_closed_orbit,
    std::vector<DynApGridPoint>& grid,
    unsigned int nr_threads
  ) {

  Status::type status = Status::success;

  // finds 6D closed-orbit
  if (calculate_closed_orbit) {
    status = calc_closed_orbit(accelerator, cod, __FUNCTION__);
    if (status != Status::success) {
      cod.clear();
      for(unsigned int i=0; i<1+accelerator.lattice.size(); ++i) cod.push_back(Pos<double>(nan("")));
    }
  }


  // creates grid with tracking points
  grid.resize(nrpts_x * nrpts_y);
  int idx = 0;
  for(unsigned int i=0; i<nrpts_x; ++i) {
    double x = x_min + i * (x_max - x_min) / (nrpts_x - 1.0);
    for(unsigned int j=0; j<nrpts_y; ++j) {
      double y = y_max + j * (y_min - y_max) / (nrpts_y - 1.0);
      // prepares initial position
      grid[idx].p = p0;     // offset
      grid[idx].p.rx += x;  // dynapt around closed-orbit
      grid[idx].p.ry += y;  // dynapt around closed-orbit
      grid[idx].start_element = 0; grid[idx].lost_turn = 0; grid[idx].lost_element = 0; grid[idx].lost_plane = Plane::no_plane;
      idx++;
    }
  }

  if (status == Status::success) {
    //std::vector<double> output;
    ThreadSharedData thread_data;
    thread_type = "xy";
    thread_data.nr_tasks = grid.size();
    thread_data.func =  thread_dynap;
    thread_nr_turns = nr_turns;
    thread_accelerator = &accelerator;
    thread_cod = &cod;
    thread_grid = &grid;
    start_all_threads(thread_data, nr_threads);
  }

  return Status::success;

}

Status::type dynap_ex(
    const Accelerator& accelerator,
    std::vector<Pos<double> >& cod,
    unsigned int nr_turns,
    const Pos<double>& p0,
    unsigned int nrpts_e, double e_min, double e_max,
    unsigned int nrpts_x, double x_min, double x_max,
    bool calculate_closed_orbit,
    std::vector<DynApGridPoint>& grid,
    unsigned int nr_threads
  ) {

  Status::type status = Status::success;

  // finds 6D closed-orbit
  if (calculate_closed_orbit) {
    status = calc_closed_orbit(accelerator, cod, __FUNCTION__);
    if (status != Status::success) {
      cod.clear();
      for(unsigned int i=0; i<1+accelerator.lattice.size(); ++i) cod.push_back(Pos<double>(nan("")));
    }
  }


  // creates grid with tracking points
  grid.resize(nrpts_e * nrpts_x);
  int idx = 0;
  for(unsigned int i=0; i<nrpts_e; ++i) {
    double e = e_min + i * (e_max - e_min) / (nrpts_e - 1.0);
    for(unsigned int j=0; j<nrpts_x; ++j) {
      double x = x_max + j * (x_min - x_max) / (nrpts_x - 1.0);
      // prepares initial position
      grid[idx].p = p0;     // offset
      grid[idx].p.de += e;  // dynapt around closed-orbit
      grid[idx].p.rx += x;  // dynapt around closed-orbit
      grid[idx].start_element = 0; grid[idx].lost_turn = 0; grid[idx].lost_element = 0; grid[idx].lost_plane = Plane::no_plane;
      idx++;
    }
  }

  if (status == Status::success) {
    //std::vector<double> output;
    ThreadSharedData thread_data;
    thread_type = "ex";
    thread_data.nr_tasks = grid.size();
    thread_data.func =  thread_dynap;
    thread_nr_turns = nr_turns;
    thread_accelerator = &accelerator;
    thread_cod = &cod;
    thread_grid = &grid;
    start_all_threads(thread_data, nr_threads);
  }

  return Status::success;

}


Status::type dynap_acceptance(
    const std::string calc_type,
    const Accelerator& accelerator,
    std::vector<Pos<double> >& cod,
    unsigned int nr_turns,
    const Pos<double>& p0,
    const double& e_init,
    const double& e_delta,
    unsigned int nr_steps_back,
    double rescale,
    unsigned int nr_iterations,
    const double& s_min, const double& s_max,
    const std::vector<std::string>& fam_names,
    bool calculate_closed_orbit,
    std::vector<DynApGridPoint>& grid,
    unsigned int nr_threads
  ) {

  Status::type status = Status::success;

  if (verbose_on) std::cout << get_timestamp() << " calc_type is '" << calc_type << "'" << std::endl;

  // finds 6D closed-orbit
  if (calculate_closed_orbit) {
    status = calc_closed_orbit(accelerator, cod, __FUNCTION__);
    if (status != Status::success) {
      cod.clear();
      for(unsigned int i=0; i<1+accelerator.lattice.size(); ++i) cod.push_back(Pos<double>(nan("")));
    }
  }

  // finds out which in which elements tracking is to be performed
  grid.clear();
  std::vector<unsigned int> elements;
  double s = 0.0;
  for(unsigned int i=0; i<accelerator.lattice.size(); ++i) {
    if ((s >= s_min) && (s <= s_max)) {
      if (std::find(fam_names.begin(), fam_names.end(), accelerator.lattice[i].fam_name) != fam_names.end()) {
        elements.push_back(i);            // calcs at start of element
        DynApGridPoint p;
        p.start_element = 0; p.lost_turn = 0; p.lost_element = 0; p.lost_plane = Plane::no_plane;
        grid.push_back(p);        // for positive pxa && pya acceptance (negative ma)
        if (calc_type == "dynap_ma") {
          grid.push_back(p); // for positive acceptance
        }
      }
    }
    s += accelerator.lattice[i].length;
  }
  if (verbose_on) std::cout << get_timestamp() << " number of elements within range is " << elements.size() << std::endl;

  if (status == Status::success) {
    //std::vector<double> output;
    ThreadSharedData thread_data;
    thread_type = calc_type;
    thread_data.nr_tasks = grid.size();
    thread_data.func =  thread_dynap_acceptance;
    thread_nr_turns = nr_turns;
    thread_accelerator = &accelerator;
    thread_cod = &cod;
    thread_grid = &grid;
    thread_ma_e_init = &e_init;
    thread_ma_e_delta = &e_delta;
    thread_ma_nr_steps_back = &nr_steps_back;
    thread_ma_rescale = &rescale;
    thread_ma_nr_iterations = &nr_iterations;
    thread_ma_p0 = &p0;
    thread_elements = &elements;
    start_all_threads(thread_data, nr_threads);
  }

  return Status::success;

}

Status::type dynap_xyfmap(
    const Accelerator& accelerator,
    std::vector<Pos<double> >& cod,
    unsigned int nr_turns,
    const Pos<double>& p0,
    unsigned int nrpts_x, double x_min, double x_max,
    unsigned int nrpts_y, double y_min, double y_max,
    bool calculate_closed_orbit,
    std::vector<DynApGridPoint>& grid,
    unsigned int nr_threads
  ) {

  Status::type status = Status::success;

  // finds 6D closed-orbit
  if (calculate_closed_orbit) {
    status = calc_closed_orbit(accelerator, cod, __FUNCTION__);
    if (status != Status::success) {
      cod.clear();
      for(unsigned int i=0; i<1+accelerator.lattice.size(); ++i) cod.push_back(Pos<double>(nan("")));
    }
  }


  // creates grid with tracking points
  grid.resize(nrpts_x * nrpts_y);
  //std::cout << "teste: " << nrpts_x << " " << nrpts_y << " " << grid.size() << std::endl;
  int idx = 0;
  for(unsigned int i=0; i<nrpts_x; ++i) {
    double x = x_min + i * (x_max - x_min) / (nrpts_x - 1.0);
    for(unsigned int j=0; j<nrpts_y; ++j) {
      double y = y_max + j * (y_min - y_max) / (nrpts_y - 1.0);
      // prepares initial position
      grid[idx].p = p0;     // offset
      grid[idx].p.rx += x;  // dynapt around closed-orbit
      grid[idx].p.ry += y;  // dynapt around closed-orbit
      grid[idx].start_element = 0; grid[idx].lost_turn = 0; grid[idx].lost_element = 0; grid[idx].lost_plane = Plane::no_plane;
      grid[idx].nux1 = grid[idx].nuy1 = 0.0;
      grid[idx].nux2 = grid[idx].nuy2 = 0.0;
      idx++;
    }
  }

  if (status == Status::success) {
    std::vector<double> output;
    thread_type = "xyfmap";
    ThreadSharedData thread_data;
    thread_data.nr_tasks = grid.size();
    thread_data.func =  thread_dynap_naff;
    thread_nr_turns = nr_turns;
    thread_accelerator = &accelerator;
    thread_cod = &cod;
    thread_grid = &grid;
    start_all_threads(thread_data, nr_threads);
  }

  return Status::success;

}

Status::type dynap_exfmap(
    const Accelerator& accelerator,
    std::vector<Pos<double> >& cod,
    unsigned int nr_turns,
    const Pos<double>& p0,
    unsigned int nrpts_e, double e_min, double e_max,
    unsigned int nrpts_x, double x_min, double x_max,
    bool calculate_closed_orbit,
    std::vector<DynApGridPoint>& grid,
    unsigned int nr_threads
  ) {

  Status::type status = Status::success;

  // finds 6D closed-orbit
  if (calculate_closed_orbit) {
    status = calc_closed_orbit(accelerator, cod, __FUNCTION__);
    if (status != Status::success) {
      cod.clear();
      for(unsigned int i=0; i<1+accelerator.lattice.size(); ++i) cod.push_back(Pos<double>(nan("")));
    }
  }

  // creates grid with tracking points
  grid.resize(nrpts_e * nrpts_x);
  int idx = 0;
  for(unsigned int i=0; i<nrpts_e; ++i) {
    double e = e_min + i * (e_max - e_min) / (nrpts_e - 1.0);
    for(unsigned int j=0; j<nrpts_x; ++j) {
      double x = x_max + j * (x_min - x_max) / (nrpts_x - 1.0);
      // prepares initial position
      grid[idx].p = p0;     // offset
      grid[idx].p.de += e;  // dynapt around closed-orbit
      grid[idx].p.rx += x;  // dynapt around closed-orbit
      grid[idx].start_element = 0; grid[idx].lost_turn = 0; grid[idx].lost_element = 0; grid[idx].lost_plane = Plane::no_plane;
      grid[idx].nux1 = grid[idx].nuy1 = 0.0;
      grid[idx].nux2 = grid[idx].nuy2 = 0.0;
      idx++;
    }
  }

  if (status == Status::success) {
    std::vector<double> output;
    thread_type = "exfmap";
    ThreadSharedData thread_data;
    thread_data.nr_tasks = grid.size();
    thread_data.func = thread_dynap_naff;
    thread_nr_turns = nr_turns;
    thread_accelerator = &accelerator;
    thread_cod = &cod;
    thread_grid = &grid;
    start_all_threads(thread_data, nr_threads);
  }


  return Status::success;

}


// implementation of auxiliary functions

static Status::type calc_closed_orbit(const Accelerator& accelerator, std::vector<Pos<double> >& cod, const char* function_name) {
  if (verbose_on) {
    std::cout << get_timestamp() << " <" << function_name << ">: calculating closed-orbit...";
    std::cout.flush();
  }

  // creates a new accelerator with vchamber off (just for finding closed orbit)
  Accelerator the_ring = accelerator;
  bool vchamber_state = the_ring.vchamber_on;
  the_ring.vchamber_on = false;
  Status::type status = track_findorbit6(the_ring, cod);
  if (status == Status::success) {
    // turns vchamber to original state && checks if found closed_orbit survives
    the_ring.vchamber_on = vchamber_state;
    status = track_findorbit6(the_ring, cod, cod[0]);
  }
  if (verbose_on) std::cout << string_error_messages[status] <<  std::endl;
  return status;
}

static void thread_dynap_naff(ThreadSharedData* thread_data, int thread_id, long task_id) {

  std::vector<DynApGridPoint>& grid = *thread_grid;

  Pos<double> p = grid[task_id].p + (*thread_cod)[0]; // adds closed-orbit
  if (fabs(p.ry) < tiny_y_amp) p.ry = sgn(p.ry) * tiny_y_amp;

  std::vector<Pos<double>> new_pos;
  Status::type lstatus = Status::success;
  lstatus = track_ringpass (*thread_accelerator,
                            p,
                            new_pos,
                            thread_nr_turns/2,
                            grid[task_id].lost_turn,
                            grid[task_id].lost_element,
                            grid[task_id].lost_plane,
                            true);
  //pthread_mutex_lock(thread_data->mutex);
  if (lstatus == Status::success) naff_traj(new_pos, grid[task_id].nux1, grid[task_id].nuy1);
  //pthread_mutex_unlock(thread_data->mutex);
  if (lstatus == Status::success) {
    p = new_pos.back();

    //pthread_mutex_lock(thread_data->mutex);
    //printf("nan thread:%02i|task:%06lu/%06lu  rx:%+.4e|ry:%+.4e\n", thread_id, (1+task_id), thread_data->nr_tasks, p.rx, p.ry);
    //pthread_mutex_unlock(thread_data->mutex);

    std::vector<Pos<double>> new_pos;
    lstatus = track_ringpass (*thread_accelerator,
                              p,
                              new_pos,
                              thread_nr_turns/2,
                              grid[task_id].lost_turn,
                              grid[task_id].lost_element,
                              grid[task_id].lost_plane,
                              true);

    //p = new_pos.back();
    //pthread_mutex_lock(thread_data->mutex);
    //printf("nan2 thread:%02i|task:%06lu/%06lu  rx:%+.4e|ry:%+.4e\n", thread_id, (1+task_id), thread_data->nr_tasks, p.rx, p.ry);
    //pthread_mutex_unlock(thread_data->mutex);

    //pthread_mutex_lock(thread_data->mutex);
    if (lstatus == Status::success) naff_traj(new_pos, grid[task_id].nux2, grid[task_id].nuy2);
    //pthread_mutex_unlock(thread_data->mutex);
  }

  if (thread_type.compare("xyfmap") == 0) {
    pthread_mutex_lock(thread_data->mutex);
    printf("thread:%02i|task:%06lu/%06lu  rx:%+.4e|ry:%+.4e  nu1:%.4e|%.4e  nu2:%.4e|%.4e  dnu:%.4e|%.4e\n", thread_id, (1+task_id), thread_data->nr_tasks, grid[task_id].p.rx, grid[task_id].p.ry, grid[task_id].nux1, grid[task_id].nuy1, grid[task_id].nux2, grid[task_id].nuy2, fabs(grid[task_id].nux2-grid[task_id].nux1), fabs(grid[task_id].nuy2-grid[task_id].nuy1));
    pthread_mutex_unlock(thread_data->mutex);
  } else if (thread_type.compare("exfmap") == 0) {
    pthread_mutex_lock(thread_data->mutex);
    printf("thread:%02i|task:%06lu/%06lu  de:%+.4e|rx:%+.4e  nu1:%.4e|%.4e  nu2:%.4e|%.4e  dnu:%.4e|%.4e\n", thread_id, (1+task_id), thread_data->nr_tasks, grid[task_id].p.de, grid[task_id].p.rx, grid[task_id].nux1, grid[task_id].nuy1, grid[task_id].nux2, grid[task_id].nuy2, fabs(grid[task_id].nux2-grid[task_id].nux1), fabs(grid[task_id].nuy2-grid[task_id].nuy1));
    pthread_mutex_unlock(thread_data->mutex);
  } else {
    pthread_mutex_lock(thread_data->mutex);
    std::cerr << "undefined multithread dynap calculation type" << std::endl;
    pthread_mutex_unlock(thread_data->mutex);
  }

}

static void thread_dynap(ThreadSharedData* thread_data, int thread_id, long task_id) {

  std::vector<DynApGridPoint>& grid = *thread_grid;

  // pthread_mutex_lock(thread_data->mutex);
  // printf("thread:%02i|task:%06lu/%06lu -> rx:%+.4e|ry:%+.4e\n", thread_id, (1+task_id), thread_data->nr_tasks, grid[task_id].p.rx, grid[task_id].p.ry);
  // pthread_mutex_unlock(thread_data->mutex);

  Pos<double> p = grid[task_id].p + (*thread_cod)[0]; // adds closed-orbit
  if (fabs(p.ry) < tiny_y_amp) p.ry = sgn(p.ry) * tiny_y_amp;

  std::vector<Pos<double>> new_pos;
  Status::type lstatus = Status::success;
  lstatus = track_ringpass (*thread_accelerator,
                            p,
                            new_pos,
                            thread_nr_turns,
                            grid[task_id].lost_turn,
                            grid[task_id].lost_element,
                            grid[task_id].lost_plane,
                            true);
  if (thread_type.compare("xy") == 0) {
    pthread_mutex_lock(thread_data->mutex);
    printf("thread:%02i|task:%06lu/%06lu  rx:%+.4e|ry:%+.4e  turn:%05i|element:%05i  status:%s\n", thread_id, (1+task_id), thread_data->nr_tasks, grid[task_id].p.rx, grid[task_id].p.ry, grid[task_id].lost_turn, grid[task_id].lost_element, string_error_messages[lstatus].c_str());
    pthread_mutex_unlock(thread_data->mutex);
  } else if (thread_type.compare("ex") == 0) {
    pthread_mutex_lock(thread_data->mutex);
    printf("thread:%02i|task:%06lu/%06lu  de:%+.4e|dx:%+.4e  turn:%05i|element:%05i  status:%s\n", thread_id, (1+task_id), thread_data->nr_tasks, grid[task_id].p.de, grid[task_id].p.rx, grid[task_id].lost_turn, grid[task_id].lost_element, string_error_messages[lstatus].c_str());
    pthread_mutex_unlock(thread_data->mutex);
  } else {
    std::cerr << "undefined multithread dynap calculation type" << std::endl;
  }

}

static void thread_dynap_acceptance(ThreadSharedData* thread_data, int thread_id, long task_id) {

  std::vector<DynApGridPoint>& grid = *thread_grid;
  const std::vector<unsigned int>& elements = *thread_elements;

  DynApGridPoint p = (*thread_grid)[task_id];
  DynApGridPoint point;
  unsigned int element_nr;
  double p_init, p_delta;


  // checks the rype of calculation
  const int ma = 0; const int pxa = 1; const int pya = 2;
  int calc_type;
  if (thread_type == "dynap_ma") {
    element_nr = task_id / 2;
    calc_type = ma;
    p_init = (*thread_ma_e_init) * ((task_id % 2) ? 1.0 : -1.0);
    p_delta = (*thread_ma_e_delta) * ((task_id % 2) ? 1.0 : -1.0);
  } else if (thread_type == "dynap_pxa") {
    element_nr = task_id;
    calc_type = pxa;
    p_init = (*thread_ma_e_init);
    p_delta = (*thread_ma_e_delta);
  } else if (thread_type == "dynap_pya") {
    element_nr = task_id;
    calc_type = pya;
    p_init = (*thread_ma_e_init);
    p_delta = (*thread_ma_e_delta);
  }

  unsigned int start_element = elements[element_nr];
  double nr_iterations = (*thread_ma_nr_iterations);
  double nr_steps_back = (*thread_ma_nr_steps_back);
  double rescale = (*thread_ma_rescale);

  double pa = p_init;
  while (true) {
    while (true) {
      //std::cout << pa << std::endl;
      point.p = *thread_ma_p0;     // offset
      switch (calc_type) {     // sets trial parameter
        case ma:  point.p.de += pa; break;
        case pxa: point.p.px += pa; break;
        case pya: point.p.py += pa; break;
      }
      point.start_element = start_element; point.lost_turn = 0; point.lost_element = start_element; point.lost_plane = Plane::no_plane;
      std::vector<Pos<double> > new_pos;
      Pos<double> p = point.p + (*thread_cod)[start_element];  // p initial condition for tracking
      if (fabs(p.ry) < tiny_y_amp) p.ry = sgn(p.ry) * tiny_y_amp;
      Status::type status = track_ringpass (*thread_accelerator, p, new_pos, thread_nr_turns, point.lost_turn, point.lost_element, point.lost_plane, false);
      if (status != Status::success) {
        pa -= p_delta;
        if (calc_type == ma)  { point.p.de = pa; break; };
        if (calc_type == pxa) { point.p.px = pa; break; };
        if (calc_type == pya) { point.p.py = pa; break; };
      }
      pa += p_delta;
    }
    if (nr_iterations < 1) break; else nr_iterations--;
    pa -= p_delta * nr_steps_back;
    p_delta *= rescale;
    pa += p_delta;
  }


  grid[task_id] = point;

  if (calc_type == ma) {
    pthread_mutex_lock(thread_data->mutex);
    printf("thread:%02i|task:%06lu/%06lu  element:%04i|de:%+.4e  %s\n", thread_id, (1+task_id), thread_data->nr_tasks, element_nr, grid[task_id].p.de, thread_accelerator->lattice[(*thread_elements)[element_nr]].fam_name.c_str());
    pthread_mutex_unlock(thread_data->mutex);
  } else if (calc_type == pxa) {
    pthread_mutex_lock(thread_data->mutex);
    printf("thread:%02i|task:%06lu/%06lu  element:%04i|px:%+.4e  %s\n", thread_id, (1+task_id), thread_data->nr_tasks, element_nr, grid[task_id].p.px, thread_accelerator->lattice[(*thread_elements)[element_nr]].fam_name.c_str());
    pthread_mutex_unlock(thread_data->mutex);
  } else if (calc_type == pya) {
    pthread_mutex_lock(thread_data->mutex);
    printf("thread:%02i|task:%06lu/%06lu  element:%04i|de:%+.4e  %s\n", thread_id, (1+task_id), thread_data->nr_tasks, element_nr, grid[task_id].p.py, thread_accelerator->lattice[(*thread_elements)[element_nr]].fam_name.c_str());
    pthread_mutex_unlock(thread_data->mutex);
  }


}
