#include "trackcpp.h"
#include "output.h"
#include "dynap.h"
#include "tracking.h"
#include "lattice.h"
#include "pos.h"
#include "auxiliary.h"
#include <algorithm>
#include <vector>
#include <cfloat>

static const double tiny_y_amp = 1e-7; // [m]

// global thread variables
static int                   thread_nr_turns    = 0;
static const Accelerator*    thread_accelerator = NULL;
std::vector<Pos<double>>*    thread_cod = NULL;
std::vector<DynApGridPoint>* thread_grid = NULL;

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
    // turns vchamber to original state and checks if found closed_orbit survives
    the_ring.vchamber_on = vchamber_state;
    status = track_findorbit6(the_ring, cod, cod[0]);
  }
  if (verbose_on) std::cout << string_error_messages[status] <<  std::endl;
  return status;
}

static DynApGridPoint find_momentum_acceptance(const Accelerator& accelerator, std::vector<Pos<double> >& cod, unsigned int nr_turns, const Pos<double>& p0, const double& e0, const double& e_tol, unsigned int element_idx);

static DynApGridPoint find_fine_momentum_acceptance(const Accelerator& accelerator, std::vector<Pos<double> >& cod, unsigned int nr_turns, const Pos<double>& p0, const double& e_init, const double& e_tol, unsigned int element_idx);

Status::type dynap_xy(
    const Accelerator& accelerator,
    std::vector<Pos<double> >& cod,
    unsigned int nr_turns,
    const Pos<double>& p0,
    unsigned int nrpts_x, double x_min, double x_max,
    unsigned int nrpts_y, double y_min, double y_max,
    bool calculate_closed_orbit,
    std::vector<DynApGridPoint>& grid
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

  char buffer[1024];

  // main loop
  if (verbose_on) std::cout << get_timestamp() << " looping over dynap_xy points" << std::endl;
  grid.resize(nrpts_x * nrpts_y);
  int idx = 0;
  for(unsigned int i=0; i<nrpts_x; ++i) {
    double x = x_min + i * (x_max - x_min) / (nrpts_x - 1.0);
    for(unsigned int j=0; j<nrpts_y; ++j) {
      double y = y_max + j * (y_min - y_max) / (nrpts_y - 1.0);

      if (verbose_on) {
        sprintf(buffer, "(%03i,%03i): x=%+11.4E [m], y=%+11.4E [m] -> ", i+1, j+1, x, y);
        std::cout << buffer; std::cout.flush();
      }

      // prepares initial position
      grid[idx].p = p0;     // offset
      grid[idx].p.rx += x;  // dynapt around closed-orbit
      grid[idx].p.ry += y;  // dynapt around closed-orbit
      grid[idx].start_element = 0; grid[idx].lost_turn = 0; grid[idx].lost_element = 0; grid[idx].lost_plane = Plane::no_plane;

      std::vector<Pos<double> > new_pos;
      Pos<double> p = grid[idx].p + cod[0]; // adds closed-orbit
      if (fabs(p.ry) < tiny_y_amp) p.ry = sgn(p.ry) * tiny_y_amp;
      Status::type lstatus = status;
      if (status == Status::success) lstatus = track_ringpass (accelerator, p, new_pos, nr_turns, grid[idx].lost_turn, grid[idx].lost_element, grid[idx].lost_plane, false);

      if (verbose_on) {
        sprintf(buffer, "turn=%05i element=%05i  %20s", grid[idx].lost_turn, grid[idx].lost_element, string_error_messages[lstatus].c_str());
        std::cout << buffer << std::endl;
      }
      idx++;
    }
  }

  if (verbose_on) std::cout << get_timestamp() << " end of dynap_xy loop" << std::endl;

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
    std::vector<DynApGridPoint>& grid
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

  char buffer[1024];

  // main loop
  if (verbose_on) std::cout << get_timestamp() << " looping over dynap_ex points" << std::endl;
  grid.resize(nrpts_e * nrpts_x);
  int idx = 0;
  for(unsigned int i=0; i<nrpts_e; ++i) {
    double e = e_min + i * (e_max - e_min) / (nrpts_e - 1.0);
    for(unsigned int j=0; j<nrpts_x; ++j) {
      double x = x_max + j * (x_min - x_max) / (nrpts_x - 1.0);

      if (verbose_on) {
        sprintf(buffer, "(%03i,%03i): e=%+11.4E, x=%+11.4E [m] -> ", i+1, j+1, e, x);
        std::cout << buffer; std::cout.flush();
      }

      // prepares initial position
      grid[idx].p = p0;     // offset
      grid[idx].p.de += e;  // dynapt around closed-orbit
      grid[idx].p.rx += x;  // dynapt around closed-orbit
      grid[idx].start_element = 0; grid[idx].lost_turn = 0; grid[idx].lost_element = 0; grid[idx].lost_plane = Plane::no_plane;

      std::vector<Pos<double> > new_pos;
      Pos<double> p = grid[idx].p + cod[0]; // adds closed-orbit
      if (fabs(p.ry) < tiny_y_amp) p.ry = sgn(p.ry) * tiny_y_amp;
      Status::type lstatus = status;
      if (status == Status::success) lstatus = track_ringpass (accelerator, p, new_pos, nr_turns, grid[idx].lost_turn, grid[idx].lost_element, grid[idx].lost_plane, false);

      if (verbose_on) {
        sprintf(buffer, "turn=%05i element=%05i  %20s", grid[idx].lost_turn, grid[idx].lost_element, string_error_messages[lstatus].c_str());
        std::cout << buffer << std::endl;
      }

      idx++;
    }
  }

  if (verbose_on) std::cout << get_timestamp() << " end of dynap_ex loop" << std::endl;

  return Status::success;

}

Status::type dynap_ma(
    const Accelerator& accelerator,
    std::vector<Pos<double> >& cod,
    unsigned int nr_turns,
    const Pos<double>& p0,
    const double& e0,
    const double& e_tol,
    const double& s_min, const double& s_max,
    const std::vector<std::string>& fam_names,
    bool calculate_closed_orbit,
    std::vector<DynApGridPoint>& grid
  ) {

  Status::type status = Status::success;

  const std::vector<Element>& the_ring = accelerator.lattice;

  // finds 6D closed-orbit
  if (calculate_closed_orbit) {
    status = calc_closed_orbit(accelerator, cod, __FUNCTION__);
    if (status != Status::success) {
      cod.clear();
      for(unsigned int i=0; i<1+accelerator.lattice.size(); ++i) cod.push_back(Pos<double>(nan("")));
    }
  }

  // finds out which in which elements tracking is to be performed
  std::vector<unsigned int> elements;
  double s = 0.0;
  for(unsigned int i=0; i<the_ring.size(); ++i) {
    if ((s >= s_min) and (s <= s_max)) {
      if (std::find(fam_names.begin(), fam_names.end(), the_ring[i].fam_name) != fam_names.end()) {
        elements.push_back(i);                      // calcs at start of element
      }
    }
    s += the_ring[i].length;
  }
  if (verbose_on) std::cout << get_timestamp() << " number of elements within range is " << elements.size() << std::endl;

  // main loop
  if (verbose_on) std::cout << get_timestamp() << " looping over dynap_ma points" << std::endl;
  grid.resize(2*elements.size());
  int idx = 0;
  for(unsigned int i=0; i<elements.size(); ++i) {
    if (verbose_on) {
      std::string timestamp = get_timestamp();
      fprintf(stdout, "%s %05i %-20s: ", timestamp.c_str(), elements[i], the_ring[elements[i]].fam_name.c_str());
      fflush(stdout);
    }

    DynApGridPoint p;

    // finds negative energy acceptance
    p = find_momentum_acceptance(accelerator, cod, nr_turns, p0, -1.0 * e0, 5 * e_tol, elements[i]);
    grid[idx++] = find_fine_momentum_acceptance(accelerator, cod, nr_turns, p0, p.p.de, e_tol, elements[i]);
    if (status != Status::success) grid[idx-1].p.de = 0;
    if (verbose_on) {
      fprintf(stdout, "%+11.4E ", grid[idx-1].p.de); fflush(stdout);
    }
    // finds positive energy acceptance
    p = find_momentum_acceptance(accelerator, cod, nr_turns, p0, +1.0 * e0, 5 * e_tol, elements[i]);
    grid[idx++] = find_fine_momentum_acceptance(accelerator, cod, nr_turns, p0, p.p.de, e_tol, elements[i]);
    if (status != Status::success) grid[idx-1].p.de = 0;
    if (verbose_on) {
      fprintf(stdout, "%+11.4E ", grid[idx-1].p.de); fflush(stdout);
    }
    std::cout << std::endl;

  }

  return Status::success;

}


void dynap_fmap_run(ThreadSharedData* thread_data, int thread_id, long task_id) {

  //std::vector<DynApGridPoint>& grid = *((std::vector<DynApGridPoint>*) thread_data->input);
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
                            thread_nr_turns/2,
                            grid[task_id].lost_turn,
                            grid[task_id].lost_element,
                            grid[task_id].lost_plane,
                            true);
  if (lstatus == Status::success) naff_run(new_pos, grid[task_id].nux1, grid[task_id].nuy1);
  if (lstatus == Status::success) {
    p = new_pos.back();
    lstatus = track_ringpass (*thread_accelerator,
                              p,
                              new_pos,
                              thread_nr_turns/2,
                              grid[task_id].lost_turn,
                              grid[task_id].lost_element,
                              grid[task_id].lost_plane,
                              true);
    if (lstatus == Status::success) naff_run(new_pos, grid[task_id].nux2, grid[task_id].nuy2);
  }

  pthread_mutex_lock(thread_data->mutex);
  printf("thread:%02i|task:%06lu/%06lu  rx:%+.4e|ry:%+.4e  nu1:%.4e|%.4e  nu2:%.4e|%.4e  dnu:%.4e|%.4e\n", thread_id, (1+task_id), thread_data->nr_tasks, grid[task_id].p.rx, grid[task_id].p.ry, grid[task_id].nux1, grid[task_id].nuy1, grid[task_id].nux2, grid[task_id].nuy2, fabs(grid[task_id].nux2-grid[task_id].nux1), fabs(grid[task_id].nuy2-grid[task_id].nuy1));
  pthread_mutex_unlock(thread_data->mutex);

}

Status::type dynap_fmap(
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
      grid[idx].nux1 = grid[idx].nuy1 = 0.0;
      grid[idx].nux2 = grid[idx].nuy2 = 0.0;
      idx++;
    }
  }

  if (status == Status::success) {
    std::vector<double> output;
    ThreadSharedData thread_data;
    thread_data.nr_tasks = grid.size();
    thread_data.func =  dynap_fmap_run;
    thread_nr_turns = nr_turns;
    thread_accelerator = &accelerator;
    thread_cod = &cod;
    thread_grid = &grid;
    start_all_threads(thread_data, nr_threads);
  }

  return Status::success;

}

static DynApGridPoint find_momentum_acceptance(
    const Accelerator& accelerator,
    std::vector<Pos<double> >& cod,
    unsigned int nr_turns,
    const Pos<double>& p0,
    const double& e0,
    const double& e_tol,
    unsigned int element_idx
  ) {

  DynApGridPoint point;

  double e_stable   = 0;
  // search initial unstable energy offset
  double e_unstable = e0;
  while (true) {
    point.p = p0;     // offset
    point.p.de += e_unstable;
    point.start_element = element_idx; point.lost_turn = 0; point.lost_element = element_idx; point.lost_plane = Plane::no_plane;
    std::vector<Pos<double> > new_pos;
    Pos<double> p = point.p + cod[element_idx];  // p initial condition for tracking
    if (fabs(p.ry) < tiny_y_amp) p.ry = sgn(p.ry) * tiny_y_amp;
    Status::type status = track_ringpass (accelerator, p, new_pos, nr_turns, point.lost_turn, point.lost_element, point.lost_plane, false);
    if (status == Status::success) {
      e_stable    = e_unstable;
      e_unstable *= 2.0;
    } else break;
  }

  unsigned int lost_element = point.lost_element;
  unsigned int lost_turn    = point.lost_turn;
  Plane::type  lost_plane   = point.lost_plane;
  while (fabs(e_unstable - e_stable) > e_tol) {
    double e = 0.5 * (e_unstable + e_stable);
    point.p = p0;     // offset
    point.p.de += e;
    point.start_element = element_idx; point.lost_turn = 0; point.lost_element = element_idx; point.lost_plane = Plane::no_plane;
    std::vector<Pos<double> > new_pos;
    Pos<double> p = point.p + cod[element_idx];  // p initial condition for tracking
    if (fabs(p.ry) < tiny_y_amp) p.ry = sgn(p.ry) * tiny_y_amp;
    Status::type status = track_ringpass (accelerator, p, new_pos, nr_turns, point.lost_turn, point.lost_element, point.lost_plane, false);
    if (status == Status::success) {
      e_stable = e;
    } else {
      e_unstable = e;
      lost_element = point.lost_element;
      lost_turn    = point.lost_turn;
      lost_plane   = point.lost_plane;
    }
  }

  // records solution
  point.start_element = element_idx;
  point.p = p0;           // offset
  point.p.de += e_stable; // conservative estimate within [e_stable, e_unstable] interval
  point.lost_element = lost_element;
  point.lost_plane   = lost_plane;
  point.lost_turn    = lost_turn;

  return point;

}

static DynApGridPoint find_fine_momentum_acceptance(
    const Accelerator& accelerator,
    std::vector<Pos<double> >& cod,
    unsigned int nr_turns,
    const Pos<double>& p0,
    const double& e_init,
    const double& e_tol,
    unsigned int element_idx
    ) {

  double e_step = e_init > 0 ? -e_tol : e_tol;
  DynApGridPoint point; point.p = p0;
  double last_unstable = e_init-e_step;
  for(unsigned int i=1; i<5; ++i) {
    point.p.de = e_init + e_step * i;
    point.start_element = element_idx; point.lost_turn = 0; point.lost_element = element_idx; point.lost_plane = Plane::no_plane;
    std::vector<Pos<double> > new_pos;
    Pos<double> p = point.p + cod[element_idx];  // p initial condition for tracking
    if (fabs(p.ry) < tiny_y_amp) p.ry = sgn(p.ry) * tiny_y_amp;
    Status::type status = track_ringpass (accelerator, p, new_pos, nr_turns, point.lost_turn, point.lost_element, point.lost_plane, false);
    if (status != Status::success) {
      //e_stable    = e_unstable;
      last_unstable = point.p.de;
    };
  }
  point.p.de = last_unstable + e_step;
  return point;
}
