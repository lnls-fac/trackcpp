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
#include <trackcpp/output.h>
#include <trackcpp/dynap.h>
#include <trackcpp/tracking.h>
#include <trackcpp/lattice.h>
#include <trackcpp/pos.h>
#include <trackcpp/auxiliary.h>
#include <algorithm>
#include <vector>
#include <cfloat>

extern void naff_run(const std::vector<Pos<double>>& data, double& tunex, double& tuney);
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

// Status::type dynap_ma(
//     const Accelerator& accelerator,
//     std::vector<Pos<double> >& cod,
//     unsigned int nr_turns,
//     const Pos<double>& p0,
//     const double& e0,
//     const double& e_tol,
//     const double& s_min, const double& s_max,
//     const std::vector<std::string>& fam_names,
//     bool calculate_closed_orbit,
//     std::vector<DynApGridPoint>& grid,
//     unsigned int nr_threads
//   ) {
//
//   Status::type status = Status::success;
//
//   // finds 6D closed-orbit
//   if (calculate_closed_orbit) {
//     status = calc_closed_orbit(accelerator, cod, __FUNCTION__);
//     if (status != Status::success) {
//       cod.clear();
//       for(unsigned int i=0; i<1+accelerator.lattice.size(); ++i) cod.push_back(Pos<double>(nan("")));
//     }
//   }
//
//
//   // finds out which in which elements tracking is to be performed
//   grid.clear();
//   std::vector<unsigned int> elements;
//   double s = 0.0;
//   for(unsigned int i=0; i<accelerator.lattice.size(); ++i) {
//     if ((s >= s_min) and (s <= s_max)) {
//       if (std::find(fam_names.begin(), fam_names.end(), accelerator.lattice[i].fam_name) != fam_names.end()) {
//         elements.push_back(i);            // calcs at start of element
//         DynApGridPoint p;
//         p.start_element = 0; p.lost_turn = 0; p.lost_element = 0; p.lost_plane = Plane::no_plane;
//         grid.push_back(p); // for negative energy acceptance
//         grid.push_back(p); // for positive energy acceptance
//       }
//     }
//     s += accelerator.lattice[i].length;
//   }
//   if (verbose_on) std::cout << get_timestamp() << " number of elements within range is " << elements.size() << std::endl;
//
//   if (status == Status::success) {
//     //std::vector<double> output;
//     ThreadSharedData thread_data;
//     thread_type = "ma";
//     thread_data.nr_tasks = grid.size();
//     thread_data.func =  thread_dynap_ma;
//     thread_nr_turns = nr_turns;
//     thread_accelerator = &accelerator;
//     thread_cod = &cod;
//     thread_grid = &grid;
//     thread_ma_e0 = &e0;
//     thread_ma_e_tol = &e_tol;
//     thread_ma_p0 = &p0;
//     thread_elements = &elements;
//     start_all_threads(thread_data, nr_threads);
//   }
//
//   return Status::success;
//
// }

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
    if ((s >= s_min) and (s <= s_max)) {
      if (std::find(fam_names.begin(), fam_names.end(), accelerator.lattice[i].fam_name) != fam_names.end()) {
        elements.push_back(i);            // calcs at start of element
        DynApGridPoint p;
        p.start_element = 0; p.lost_turn = 0; p.lost_element = 0; p.lost_plane = Plane::no_plane;
        grid.push_back(p);        // for positive acceptance
        if (calc_type == "dynap_ma") {
          //std::cout << e_init << std::endl;
          grid.push_back(p); // for negative acceptance
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


// Status::type dynap_ma(
//     const Accelerator& accelerator,
//     std::vector<Pos<double> >& cod,
//     unsigned int nr_turns,
//     const Pos<double>& p0,
//     const double& e_init,
//     const double& e_delta,
//     unsigned int nr_steps_back,
//     double rescale,
//     unsigned int nr_iterations,
//     const double& s_min, const double& s_max,
//     const std::vector<std::string>& fam_names,
//     bool calculate_closed_orbit,
//     std::vector<DynApGridPoint>& grid,
//     unsigned int nr_threads
//   ) {
//
//   Status::type status = Status::success;
//
//   // finds 6D closed-orbit
//   if (calculate_closed_orbit) {
//     status = calc_closed_orbit(accelerator, cod, __FUNCTION__);
//     if (status != Status::success) {
//       cod.clear();
//       for(unsigned int i=0; i<1+accelerator.lattice.size(); ++i) cod.push_back(Pos<double>(nan("")));
//     }
//   }
//
//
//   // finds out which in which elements tracking is to be performed
//   grid.clear();
//   std::vector<unsigned int> elements;
//   double s = 0.0;
//   for(unsigned int i=0; i<accelerator.lattice.size(); ++i) {
//     if ((s >= s_min) and (s <= s_max)) {
//       if (std::find(fam_names.begin(), fam_names.end(), accelerator.lattice[i].fam_name) != fam_names.end()) {
//         elements.push_back(i);            // calcs at start of element
//         DynApGridPoint p;
//         p.start_element = 0; p.lost_turn = 0; p.lost_element = 0; p.lost_plane = Plane::no_plane;
//         grid.push_back(p); // for negative energy acceptance
//         grid.push_back(p); // for positive energy acceptance
//       }
//     }
//     s += accelerator.lattice[i].length;
//   }
//   if (verbose_on) std::cout << get_timestamp() << " number of elements within range is " << elements.size() << std::endl;
//
//   if (status == Status::success) {
//     //std::vector<double> output;
//     ThreadSharedData thread_data;
//     thread_type = "ma";
//     thread_data.nr_tasks = grid.size();
//     thread_data.func =  thread_dynap_ma;
//     thread_nr_turns = nr_turns;
//     thread_accelerator = &accelerator;
//     thread_cod = &cod;
//     thread_grid = &grid;
//     thread_ma_e_init = &e_init;
//     thread_ma_e_delta = &e_delta;
//     thread_ma_nr_steps_back = &nr_steps_back;
//     thread_ma_rescale = &rescale;
//     thread_ma_nr_iterations = &nr_iterations;
//     thread_ma_p0 = &p0;
//     thread_elements = &elements;
//     start_all_threads(thread_data, nr_threads);
//   }
//
//   return Status::success;
//
// }

// Status::type dynap_pxa_old(
//     const Accelerator& accelerator,
//     std::vector<Pos<double> >& cod,
//     unsigned int nr_turns,
//     const Pos<double>& p0,
//     const double& px0,
//     const double& px_tol,
//     const double& s_min, const double& s_max,
//     const std::vector<std::string>& fam_names,
//     bool calculate_closed_orbit,
//     std::vector<DynApGridPoint>& grid,
//     unsigned int nr_threads
//   ) {
//
//   Status::type status = Status::success;
//
//   // finds 6D closed-orbit
//   if (calculate_closed_orbit) {
//     status = calc_closed_orbit(accelerator, cod, __FUNCTION__);
//     if (status != Status::success) {
//       cod.clear();
//       for(unsigned int i=0; i<1+accelerator.lattice.size(); ++i) cod.push_back(Pos<double>(nan("")));
//     }
//   }
//
//
//   // finds out which in which elements tracking is to be performed
//   grid.clear();
//   std::vector<unsigned int> elements;
//   double s = 0.0;
//   for(unsigned int i=0; i<accelerator.lattice.size(); ++i) {
//     if ((s >= s_min) and (s <= s_max)) {
//       if (std::find(fam_names.begin(), fam_names.end(), accelerator.lattice[i].fam_name) != fam_names.end()) {
//         elements.push_back(i);            // calcs at start of element
//         DynApGridPoint p;
//         p.start_element = 0; p.lost_turn = 0; p.lost_element = 0; p.lost_plane = Plane::no_plane;
//         grid.push_back(p);
//       }
//     }
//     s += accelerator.lattice[i].length;
//   }
//   if (verbose_on) std::cout << get_timestamp() << " number of elements within range is " << elements.size() << std::endl;
//
//   if (status == Status::success) {
//     //std::vector<double> output;
//     ThreadSharedData thread_data;
//     thread_type = "pxa";
//     thread_data.nr_tasks = grid.size();
//     thread_data.func =  thread_dynap_pxa;
//     thread_nr_turns = nr_turns;
//     thread_accelerator = &accelerator;
//     thread_cod = &cod;
//     thread_grid = &grid;
//     thread_ma_e0 = &px0;
//     thread_ma_e_tol = &px_tol;
//     thread_ma_p0 = &p0;
//     thread_elements = &elements;
//     start_all_threads(thread_data, nr_threads);
//   }
//
//   return Status::success;
//
// }

// Status::type dynap_pxa(
//     const Accelerator& accelerator,
//     std::vector<Pos<double> >& cod,
//     unsigned int nr_turns,
//     const Pos<double>& p0,
//     const double& px_init,
//     const double& px_delta,
//     unsigned int nr_steps_back,
//     double rescale,
//     unsigned int nr_iterations,
//     const double& s_min, const double& s_max,
//     const std::vector<std::string>& fam_names,
//     bool calculate_closed_orbit,
//     std::vector<DynApGridPoint>& grid,
//     unsigned int nr_threads
//   ) {
//
//     Status::type status = Status::success;
//
//     // finds 6D closed-orbit
//     if (calculate_closed_orbit) {
//       status = calc_closed_orbit(accelerator, cod, __FUNCTION__);
//       if (status != Status::success) {
//         cod.clear();
//         for(unsigned int i=0; i<1+accelerator.lattice.size(); ++i) cod.push_back(Pos<double>(nan("")));
//       }
//     }
//
//
//     // finds out which in which elements tracking is to be performed
//     grid.clear();
//     std::vector<unsigned int> elements;
//     double s = 0.0;
//     for(unsigned int i=0; i<accelerator.lattice.size(); ++i) {
//       if ((s >= s_min) and (s <= s_max)) {
//         if (std::find(fam_names.begin(), fam_names.end(), accelerator.lattice[i].fam_name) != fam_names.end()) {
//           elements.push_back(i);            // calcs at start of element
//           DynApGridPoint p;
//           p.start_element = 0; p.lost_turn = 0; p.lost_element = 0; p.lost_plane = Plane::no_plane;
//           grid.push_back(p);
//         }
//       }
//       s += accelerator.lattice[i].length;
//     }
//     if (verbose_on) std::cout << get_timestamp() << " number of elements within range is " << elements.size() << std::endl;
//
//     if (status == Status::success) {
//       //std::vector<double> output;
//       ThreadSharedData thread_data;
//       thread_type = "pxa";
//       thread_data.nr_tasks = grid.size();
//       thread_data.func =  thread_dynap_pxa;
//       thread_nr_turns = nr_turns;
//       thread_accelerator = &accelerator;
//       thread_cod = &cod;
//       thread_grid = &grid;
//       thread_ma_e0 = &px_init;
//       thread_ma_e_tol = &px_delta;
//       thread_ma_nr_steps_back = &nr_steps_back;
//       thread_ma_rescale = &rescale;
//       thread_ma_nr_iterations = &nr_iterations;
//       thread_ma_p0 = &p0;
//       thread_elements = &elements;
//       start_all_threads(thread_data, nr_threads);
//     }
//
//     return Status::success;
//
//   }

// Status::type dynap_pya(
//     const Accelerator& accelerator,
//     std::vector<Pos<double> >& cod,
//     unsigned int nr_turns,
//     const Pos<double>& p0,
//     const double& py0,
//     const double& py_tol,
//     const double& s_min, const double& s_max,
//     const std::vector<std::string>& fam_names,
//     bool calculate_closed_orbit,
//     std::vector<DynApGridPoint>& grid,
//     unsigned int nr_threads
//   ) {
//
//   Status::type status = Status::success;
//
//   // finds 6D closed-orbit
//   if (calculate_closed_orbit) {
//     status = calc_closed_orbit(accelerator, cod, __FUNCTION__);
//     if (status != Status::success) {
//       cod.clear();
//       for(unsigned int i=0; i<1+accelerator.lattice.size(); ++i) cod.push_back(Pos<double>(nan("")));
//     }
//   }
//
//
//   // finds out which in which elements tracking is to be performed
//   grid.clear();
//   std::vector<unsigned int> elements;
//   double s = 0.0;
//   for(unsigned int i=0; i<accelerator.lattice.size(); ++i) {
//     if ((s >= s_min) and (s <= s_max)) {
//       if (std::find(fam_names.begin(), fam_names.end(), accelerator.lattice[i].fam_name) != fam_names.end()) {
//         elements.push_back(i);            // calcs at start of element
//         DynApGridPoint p;
//         p.start_element = 0; p.lost_turn = 0; p.lost_element = 0; p.lost_plane = Plane::no_plane;
//         grid.push_back(p);
//       }
//     }
//     s += accelerator.lattice[i].length;
//   }
//   if (verbose_on) std::cout << get_timestamp() << " number of elements within range is " << elements.size() << std::endl;
//
//   if (status == Status::success) {
//     //std::vector<double> output;
//     ThreadSharedData thread_data;
//     thread_type = "pya";
//     thread_data.nr_tasks = grid.size();
//     thread_data.func =  thread_dynap_pya;
//     thread_nr_turns = nr_turns;
//     thread_accelerator = &accelerator;
//     thread_cod = &cod;
//     thread_grid = &grid;
//     thread_ma_e0 = &py0;
//     thread_ma_e_tol = &py_tol;
//     thread_ma_p0 = &p0;
//     thread_elements = &elements;
//     start_all_threads(thread_data, nr_threads);
//   }
//
//   return Status::success;
//
// }

// Status::type dynap_pya(
//     const Accelerator& accelerator,
//     std::vector<Pos<double> >& cod,
//     unsigned int nr_turns,
//     const Pos<double>& p0,
//     const double& py_init,
//     const double& py_delta,
//     unsigned int nr_steps_back,
//     double rescale,
//     unsigned int nr_iterations,
//     const double& s_min, const double& s_max,
//     const std::vector<std::string>& fam_names,
//     bool calculate_closed_orbit,
//     std::vector<DynApGridPoint>& grid,
//     unsigned int nr_threads) {
//
//   Status::type status = Status::success;
//
//   // finds 6D closed-orbit
//   if (calculate_closed_orbit) {
//     status = calc_closed_orbit(accelerator, cod, __FUNCTION__);
//     if (status != Status::success) {
//       cod.clear();
//       for(unsigned int i=0; i<1+accelerator.lattice.size(); ++i) cod.push_back(Pos<double>(nan("")));
//     }
//   }
//
//
//   // finds out which in which elements tracking is to be performed
//   grid.clear();
//   std::vector<unsigned int> elements;
//   double s = 0.0;
//   for(unsigned int i=0; i<accelerator.lattice.size(); ++i) {
//     if ((s >= s_min) and (s <= s_max)) {
//       if (std::find(fam_names.begin(), fam_names.end(), accelerator.lattice[i].fam_name) != fam_names.end()) {
//         elements.push_back(i);            // calcs at start of element
//         DynApGridPoint p;
//         p.start_element = 0; p.lost_turn = 0; p.lost_element = 0; p.lost_plane = Plane::no_plane;
//         grid.push_back(p);
//       }
//     }
//     s += accelerator.lattice[i].length;
//   }
//   if (verbose_on) std::cout << get_timestamp() << " number of elements within range is " << elements.size() << std::endl;
//
//   if (status == Status::success) {
//     //std::vector<double> output;
//     ThreadSharedData thread_data;
//     thread_type = "pya";
//     thread_data.nr_tasks = grid.size();
//     thread_data.func =  thread_dynap_pya;
//     thread_nr_turns = nr_turns;
//     thread_accelerator = &accelerator;
//     thread_cod = &cod;
//     thread_grid = &grid;
//     thread_ma_e0 = &py_init;
//     thread_ma_e_tol = &py_delta;
//     thread_ma_nr_steps_back = &nr_steps_back;
//     thread_ma_rescale = &rescale;
//     thread_ma_nr_iterations = &nr_iterations;
//     thread_ma_p0 = &p0;
//     thread_elements = &elements;
//     start_all_threads(thread_data, nr_threads);
//   }
//
//   return Status::success;
//
// }


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
    // turns vchamber to original state and checks if found closed_orbit survives
    the_ring.vchamber_on = vchamber_state;
    status = track_findorbit6(the_ring, cod, cod[0]);
  }
  if (verbose_on) std::cout << string_error_messages[status] <<  std::endl;
  return status;
}

// static DynApGridPoint find_momentum_acceptance(
//   const Accelerator& accelerator,
//   const std::vector<Pos<double> >& cod,
//   unsigned int nr_turns,
//   const Pos<double>& p0,
//   double e0,
//   double e_tol,
//   unsigned int element_idx) {
//
//   DynApGridPoint point;
//
//   double e_stable   = 0;
//   // search initial unstable energy offset
//   double e_unstable = e0;
//   while (true) {
//     point.p = p0;     // offset
//     point.p.de += e_unstable;
//     point.start_element = element_idx; point.lost_turn = 0; point.lost_element = element_idx; point.lost_plane = Plane::no_plane;
//     std::vector<Pos<double> > new_pos;
//     Pos<double> p = point.p + cod[element_idx];  // p initial condition for tracking
//     if (fabs(p.ry) < tiny_y_amp) p.ry = sgn(p.ry) * tiny_y_amp;
//     Status::type status = track_ringpass (accelerator, p, new_pos, nr_turns, point.lost_turn, point.lost_element, point.lost_plane, false);
//     if (status == Status::success) {
//       e_stable    = e_unstable;
//       e_unstable *= 2.0;
//     } else break;
//   }
//
//   unsigned int lost_element = point.lost_element;
//   unsigned int lost_turn    = point.lost_turn;
//   Plane::type  lost_plane   = point.lost_plane;
//   while (fabs(e_unstable - e_stable) > e_tol) {
//     double e = 0.5 * (e_unstable + e_stable);
//     point.p = p0;     // offset
//     point.p.de += e;
//     point.start_element = element_idx; point.lost_turn = 0; point.lost_element = element_idx; point.lost_plane = Plane::no_plane;
//     std::vector<Pos<double> > new_pos;
//     Pos<double> p = point.p + cod[element_idx];  // p initial condition for tracking
//     if (fabs(p.ry) < tiny_y_amp) p.ry = sgn(p.ry) * tiny_y_amp;
//     Status::type status = track_ringpass (accelerator, p, new_pos, nr_turns, point.lost_turn, point.lost_element, point.lost_plane, false);
//     if (status == Status::success) {
//       e_stable = e;
//     } else {
//       e_unstable = e;
//       lost_element = point.lost_element;
//       lost_turn    = point.lost_turn;
//       lost_plane   = point.lost_plane;
//     }
//   }
//
//   // records solution
//   point.start_element = element_idx;
//   point.p = p0;           // offset
//   point.p.de += e_stable; // conservative estimate within [e_stable, e_unstable] interval
//   point.lost_element = lost_element;
//   point.lost_plane   = lost_plane;
//   point.lost_turn    = lost_turn;
//
//   return point;
//
// }

// static DynApGridPoint find_px_acceptance(
//   const Accelerator& accelerator,
//   const std::vector<Pos<double> >& cod,
//   unsigned int nr_turns,
//   const Pos<double>& p0,
//   double px0,
//   double px_tol,
//   unsigned int element_idx) {
//
//   DynApGridPoint point;
//
//   double px_stable   = 0;
//   // search initial unstable energy offset
//   double px_unstable = px0;
//   while (true) {
//     point.p = p0;     // offset
//     point.p.px += px_unstable;
//     point.start_element = element_idx; point.lost_turn = 0; point.lost_element = element_idx; point.lost_plane = Plane::no_plane;
//     std::vector<Pos<double> > new_pos;
//     Pos<double> p = point.p + cod[element_idx];  // p initial condition for tracking
//     if (fabs(p.ry) < tiny_y_amp) p.ry = sgn(p.ry) * tiny_y_amp;
//     Status::type status = track_ringpass (accelerator, p, new_pos, nr_turns, point.lost_turn, point.lost_element, point.lost_plane, false);
//     if (status == Status::success) {
//       px_stable    = px_unstable;
//       px_unstable *= 2.0;
//     } else break;
//   }
//
//   unsigned int lost_element = point.lost_element;
//   unsigned int lost_turn    = point.lost_turn;
//   Plane::type  lost_plane   = point.lost_plane;
//   while (fabs(px_unstable - px_stable) > px_tol) {
//     double px = 0.5 * (px_unstable + px_stable);
//     point.p = p0;     // offset
//     point.p.px += px;
//     point.start_element = element_idx; point.lost_turn = 0; point.lost_element = element_idx; point.lost_plane = Plane::no_plane;
//     std::vector<Pos<double> > new_pos;
//     Pos<double> p = point.p + cod[element_idx];  // p initial condition for tracking
//     if (fabs(p.ry) < tiny_y_amp) p.ry = sgn(p.ry) * tiny_y_amp;
//     Status::type status = track_ringpass (accelerator, p, new_pos, nr_turns, point.lost_turn, point.lost_element, point.lost_plane, false);
//     if (status == Status::success) {
//       px_stable = px;
//     } else {
//       px_unstable = px;
//       lost_element = point.lost_element;
//       lost_turn    = point.lost_turn;
//       lost_plane   = point.lost_plane;
//     }
//   }
//
//   // records solution
//   point.start_element = element_idx;
//   point.p = p0;           // offset
//   point.p.px += px_stable; // conservative estimate within [e_stable, e_unstable] interval
//   point.lost_element = lost_element;
//   point.lost_plane   = lost_plane;
//   point.lost_turn    = lost_turn;
//
//   return point;
//
// }

// static DynApGridPoint find_py_acceptance(
//   const Accelerator& accelerator,
//   const std::vector<Pos<double> >& cod,
//   unsigned int nr_turns,
//   const Pos<double>& p0,
//   double py0,
//   double py_tol,
//   unsigned int element_idx) {
//
//   DynApGridPoint point;
//
//   double py_stable   = 0;
//   // search initial unstable energy offset
//   double py_unstable = py0;
//   while (true) {
//     point.p = p0;     // offset
//     point.p.py += py_unstable;
//     point.start_element = element_idx; point.lost_turn = 0; point.lost_element = element_idx; point.lost_plane = Plane::no_plane;
//     std::vector<Pos<double> > new_pos;
//     Pos<double> p = point.p + cod[element_idx];  // p initial condition for tracking
//     if (fabs(p.ry) < tiny_y_amp) p.ry = sgn(p.ry) * tiny_y_amp;
//     Status::type status = track_ringpass (accelerator, p, new_pos, nr_turns, point.lost_turn, point.lost_element, point.lost_plane, false);
//     if (status == Status::success) {
//       py_stable    = py_unstable;
//       py_unstable *= 2.0;
//     } else break;
//   }
//
//   unsigned int lost_element = point.lost_element;
//   unsigned int lost_turn    = point.lost_turn;
//   Plane::type  lost_plane   = point.lost_plane;
//   while (fabs(py_unstable - py_stable) > py_tol) {
//     double py = 0.5 * (py_unstable + py_stable);
//     point.p = p0;     // offset
//     point.p.py += py;
//     point.start_element = element_idx; point.lost_turn = 0; point.lost_element = element_idx; point.lost_plane = Plane::no_plane;
//     std::vector<Pos<double> > new_pos;
//     Pos<double> p = point.p + cod[element_idx];  // p initial condition for tracking
//     if (fabs(p.ry) < tiny_y_amp) p.ry = sgn(p.ry) * tiny_y_amp;
//     Status::type status = track_ringpass (accelerator, p, new_pos, nr_turns, point.lost_turn, point.lost_element, point.lost_plane, false);
//     if (status == Status::success) {
//       py_stable = py;
//     } else {
//       py_unstable = py;
//       lost_element = point.lost_element;
//       lost_turn    = point.lost_turn;
//       lost_plane   = point.lost_plane;
//     }
//   }
//
//   // records solution
//   point.start_element = element_idx;
//   point.p = p0;           // offset
//   point.p.py += py_stable; // conservative estimate within [e_stable, e_unstable] interval
//   point.lost_element = lost_element;
//   point.lost_plane   = lost_plane;
//   point.lost_turn    = lost_turn;
//
//   return point;
//
// }

// static DynApGridPoint find_fine_momentum_acceptance(
//   const Accelerator& accelerator,
//   const std::vector<Pos<double> >& cod,
//   unsigned int nr_turns,
//   const Pos<double>& p0,
//   double e_init,
//   double e_tol,
//   unsigned int element_idx) {
//
//   double e_step = e_init > 0 ? -e_tol : e_tol;
//   DynApGridPoint point; point.p = p0;
//   double last_unstable = e_init-e_step;
//   for(unsigned int i=1; i<5; ++i) {
//     point.p.de = e_init + e_step * i;
//     point.start_element = element_idx; point.lost_turn = 0; point.lost_element = element_idx; point.lost_plane = Plane::no_plane;
//     std::vector<Pos<double> > new_pos;
//     Pos<double> p = point.p + cod[element_idx];  // p initial condition for tracking
//     if (fabs(p.ry) < tiny_y_amp) p.ry = sgn(p.ry) * tiny_y_amp;
//     Status::type status = track_ringpass (accelerator, p, new_pos, nr_turns, point.lost_turn, point.lost_element, point.lost_plane, false);
//     if (status != Status::success) {
//       //e_stable    = e_unstable;
//       last_unstable = point.p.de;
//     };
//   }
//   point.p.de = last_unstable + e_step;
//   return point;
// }

// static DynApGridPoint find_fine_px_acceptance(
//   const Accelerator& accelerator,
//   const std::vector<Pos<double> >& cod,
//   unsigned int nr_turns,
//   const Pos<double>& p0,
//   double px_init,
//   double px_tol,
//   unsigned int element_idx) {
//
//   double px_step = px_init > 0 ? -px_tol : px_tol;
//   DynApGridPoint point; point.p = p0;
//   double last_unstable = px_init-px_step;
//   for(unsigned int i=1; i<5; ++i) {
//     point.p.px = px_init + px_step * i;
//     point.start_element = element_idx; point.lost_turn = 0; point.lost_element = element_idx; point.lost_plane = Plane::no_plane;
//     std::vector<Pos<double> > new_pos;
//     Pos<double> p = point.p + cod[element_idx];  // p initial condition for tracking
//     if (fabs(p.ry) < tiny_y_amp) p.ry = sgn(p.ry) * tiny_y_amp;
//     Status::type status = track_ringpass (accelerator, p, new_pos, nr_turns, point.lost_turn, point.lost_element, point.lost_plane, false);
//     if (status != Status::success) {
//       //e_stable    = e_unstable;
//       last_unstable = point.p.px;
//     };
//   }
//   point.p.px = last_unstable + px_step;
//   return point;
// }

// static DynApGridPoint find_fine_py_acceptance(
//   const Accelerator& accelerator,
//   const std::vector<Pos<double> >& cod,
//   unsigned int nr_turns,
//   const Pos<double>& p0,
//   double py_init,
//   double py_tol,
//   unsigned int element_idx) {
//
//   double py_step = py_init > 0 ? -py_tol : py_tol;
//   DynApGridPoint point; point.p = p0;
//   double last_unstable = py_init-py_step;
//   for(unsigned int i=1; i<5; ++i) {
//     point.p.py = py_init + py_step * i;
//     point.start_element = element_idx; point.lost_turn = 0; point.lost_element = element_idx; point.lost_plane = Plane::no_plane;
//     std::vector<Pos<double> > new_pos;
//     Pos<double> p = point.p + cod[element_idx];  // p initial condition for tracking
//     if (fabs(p.ry) < tiny_y_amp) p.ry = sgn(p.ry) * tiny_y_amp;
//     Status::type status = track_ringpass (accelerator, p, new_pos, nr_turns, point.lost_turn, point.lost_element, point.lost_plane, false);
//     if (status != Status::success) {
//       //e_stable    = e_unstable;
//       last_unstable = point.p.py;
//     };
//   }
//   point.p.py = last_unstable + py_step;
//   return point;
// }

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
  if (lstatus == Status::success) naff_run(new_pos, grid[task_id].nux1, grid[task_id].nuy1);
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
    if (lstatus == Status::success) naff_run(new_pos, grid[task_id].nux2, grid[task_id].nuy2);
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

// static void thread_dynap_ma(ThreadSharedData* thread_data, int thread_id, long task_id) {
//
//   std::vector<DynApGridPoint>& grid = *thread_grid;
//   const std::vector<unsigned int>& elements = *thread_elements;
//
//   DynApGridPoint p = (*thread_grid)[task_id];
//
//   unsigned int element_nr = task_id / 2;
//   double e0 = (*thread_ma_e0) * ((task_id % 2) ? 1.0 : -1.0);
//
//   p = find_momentum_acceptance(*thread_accelerator,
//                                *thread_cod,
//                                thread_nr_turns,
//                                *thread_ma_p0,
//                                e0,
//                                5 * (*thread_ma_e_tol),
//                                (*thread_elements)[element_nr]);
//   p = find_fine_momentum_acceptance(*thread_accelerator,
//                                     *thread_cod,
//                                     thread_nr_turns,
//                                     *thread_ma_p0,
//                                     p.p.de,
//                                     *thread_ma_e_tol,
//                                     (*thread_elements)[element_nr]);
//   grid[task_id] = p;
//
//   pthread_mutex_lock(thread_data->mutex);
//   //printf("thread:%02i|task:%06lu/%06lu  %+.4e\n", thread_id, (1+task_id), thread_data->nr_tasks, p.p.de);
//   printf("thread:%02i|task:%06lu/%06lu  element:%04i|de:%+.4e  %s\n", thread_id, (1+task_id), thread_data->nr_tasks, element_nr, grid[task_id].p.de, thread_accelerator->lattice[(*thread_elements)[element_nr]].fam_name.c_str());
//   pthread_mutex_unlock(thread_data->mutex);
//
// }

// static void thread_dynap_ma(ThreadSharedData* thread_data, int thread_id, long task_id) {
//
//   std::vector<DynApGridPoint>& grid = *thread_grid;
//   const std::vector<unsigned int>& elements = *thread_elements;
//
//   DynApGridPoint p = (*thread_grid)[task_id];
//   DynApGridPoint point;
//   unsigned int element_nr = task_id / 2;
//   unsigned int start_element = elements[element_nr];
//
//   double e_init = (*thread_ma_e_init) * ((task_id % 2) ? 1.0 : -1.0);
//   double e_delta = (*thread_ma_e_delta) * ((task_id % 2) ? 1.0 : -1.0);
//   double nr_iterations = (*thread_ma_nr_iterations);
//   double nr_steps_back = (*thread_ma_nr_steps_back);
//   double rescale = (*thread_ma_rescale);
//
//   double e = e_init;
//   while (true) {
//     while (true) {
//       //if (task_id == 0) std::cout << e << std::endl;
//       point.p = *thread_ma_p0;     // offset
//       point.p.de += e;             // sets trial energy deviation
//       point.start_element = start_element; point.lost_turn = 0; point.lost_element = start_element; point.lost_plane = Plane::no_plane;
//       std::vector<Pos<double> > new_pos;
//       Pos<double> p = point.p + (*thread_cod)[start_element];  // p initial condition for tracking
//       if (fabs(p.ry) < tiny_y_amp) p.ry = sgn(p.ry) * tiny_y_amp;
//       Status::type status = track_ringpass (*thread_accelerator, p, new_pos, thread_nr_turns, point.lost_turn, point.lost_element, point.lost_plane, false);
//       if (status != Status::success) { e -= e_delta; point.p.de = e; break; }
//       e += e_delta;
//     }
//     if (nr_iterations < 1) break; else nr_iterations--;
//     e -= e_delta * nr_steps_back;
//     e_delta *= rescale;
//     e += e_delta;
//   }
//
//   grid[task_id] = point;
//
//   pthread_mutex_lock(thread_data->mutex);
//   printf("thread:%02i|task:%06lu/%06lu  element:%04i|de:%+.4e  %s\n", thread_id, (1+task_id), thread_data->nr_tasks, element_nr, grid[task_id].p.de, thread_accelerator->lattice[(*thread_elements)[element_nr]].fam_name.c_str());
//   pthread_mutex_unlock(thread_data->mutex);
//
// }

static void thread_dynap_acceptance(ThreadSharedData* thread_data, int thread_id, long task_id) {

  std::vector<DynApGridPoint>& grid = *thread_grid;
  const std::vector<unsigned int>& elements = *thread_elements;

  DynApGridPoint p = (*thread_grid)[task_id];
  DynApGridPoint point;
  unsigned int element_nr = task_id / 2;
  unsigned int start_element = elements[element_nr];
  double p_init, p_delta;



  // checks the rype of calculation
  const int ma = 0; const int pxa = 1; const int pya = 2;
  int calc_type;
  if (thread_type == "dynap_ma") {
    calc_type = ma;
    p_init = (*thread_ma_e_init) * ((task_id % 2) ? -1.0 : 1.0);
    p_delta = (*thread_ma_e_delta) * ((task_id % 2) ? -1.0 : 1.0);
  } else if (thread_type == "dynap_pxa") {
    calc_type = pxa;
    p_init = (*thread_ma_e_init);
    p_delta = (*thread_ma_e_delta);
  } else if (thread_type == "dynap_pya") {
    calc_type = pya;
    p_init = (*thread_ma_e_init);
    p_delta = (*thread_ma_e_delta);
  }

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


// static void thread_dynap_pxa(ThreadSharedData* thread_data, int thread_id, long task_id) {
//
//   std::vector<DynApGridPoint>& grid = *thread_grid;
//   const std::vector<unsigned int>& elements = *thread_elements;
//
//   DynApGridPoint p = (*thread_grid)[task_id];
//   DynApGridPoint point;
//   unsigned int element_nr = task_id;
//   unsigned int start_element = elements[element_nr];
//
//
//   double px_init = (*thread_ma_e_init);
//   double px_delta = (*thread_ma_e_delta);
//   double nr_iterations = (*thread_ma_nr_iterations);
//   double nr_steps_back = (*thread_ma_nr_steps_back);
//   double rescale = (*thread_ma_rescale);
//
//   double px = px_init;
//   while (true) {
//     while (true) {
//       point.p = *thread_ma_p0;     // offset
//       point.p.px += px;            // sets trial energy deviation
//       point.start_element = start_element; point.lost_turn = 0; point.lost_element = start_element; point.lost_plane = Plane::no_plane;
//       std::vector<Pos<double> > new_pos;
//       Pos<double> p = point.p + (*thread_cod)[start_element];  // p initial condition for tracking
//       if (fabs(p.ry) < tiny_y_amp) p.ry = sgn(p.ry) * tiny_y_amp;
//       Status::type status = track_ringpass (*thread_accelerator, p, new_pos, thread_nr_turns, point.lost_turn, point.lost_element, point.lost_plane, false);
//       if (status != Status::success) { px -= px_delta; point.p.px = px; break; }
//       px += px_delta;
//     }
//     if (nr_iterations < 1) break; else nr_iterations--;
//     px -= px_delta * nr_steps_back;
//     px_delta *= rescale;
//     px += px_delta;
//   }
//
//   grid[task_id] = point;
//
//   pthread_mutex_lock(thread_data->mutex);
//   printf("thread:%02i|task:%06lu/%06lu  element:%04i|px:%+.4e  %s\n", thread_id, (1+task_id), thread_data->nr_tasks, element_nr, grid[task_id].p.px, thread_accelerator->lattice[(*thread_elements)[element_nr]].fam_name.c_str());
//   pthread_mutex_unlock(thread_data->mutex);
//
//
// }

// static void thread_dynap_pxa_old(ThreadSharedData* thread_data, int thread_id, long task_id) {
//
//   std::vector<DynApGridPoint>& grid = *thread_grid;
//   const std::vector<unsigned int>& elements = *thread_elements;
//
//   DynApGridPoint p = (*thread_grid)[task_id];
//
//   unsigned int element_nr = task_id;
//   double px0 = (*thread_ma_e0);
//
//   p = find_px_acceptance(*thread_accelerator,
//                          *thread_cod,
//                          thread_nr_turns,
//                          *thread_ma_p0,
//                          px0,
//                          5 * (*thread_ma_e_tol),
//                          (*thread_elements)[element_nr]);
//   p = find_fine_px_acceptance(*thread_accelerator,
//                               *thread_cod,
//                               thread_nr_turns,
//                               *thread_ma_p0,
//                               p.p.px,
//                               *thread_ma_e_tol,
//                               (*thread_elements)[element_nr]);
//   grid[task_id] = p;
//
//   pthread_mutex_lock(thread_data->mutex);
//   printf("thread:%02i|task:%06lu/%06lu  element:%04i|px:%+.4e  %s\n", thread_id, (1+task_id), thread_data->nr_tasks, element_nr, grid[task_id].p.px, thread_accelerator->lattice[(*thread_elements)[element_nr]].fam_name.c_str());
//   pthread_mutex_unlock(thread_data->mutex);
//
// }

// static void thread_dynap_pya(ThreadSharedData* thread_data, int thread_id, long task_id) {
//
//   std::vector<DynApGridPoint>& grid = *thread_grid;
//   const std::vector<unsigned int>& elements = *thread_elements;
//
//   DynApGridPoint p = (*thread_grid)[task_id];
//
//   unsigned int element_nr = task_id;
//   double py0 = (*thread_ma_e0);
//
//   p = find_py_acceptance(*thread_accelerator,
//                          *thread_cod,
//                          thread_nr_turns,
//                          *thread_ma_p0,
//                          py0,
//                          5 * (*thread_ma_e_tol),
//                          (*thread_elements)[element_nr]);
//   p = find_fine_py_acceptance(*thread_accelerator,
//                               *thread_cod,
//                               thread_nr_turns,
//                               *thread_ma_p0,
//                               p.p.py,
//                               *thread_ma_e_tol,
//                               (*thread_elements)[element_nr]);
//   grid[task_id] = p;
//
//   pthread_mutex_lock(thread_data->mutex);
//   printf("thread:%02i|task:%06lu/%06lu  element:%04i|py:%+.4e  %s\n", thread_id, (1+task_id), thread_data->nr_tasks, element_nr, grid[task_id].p.py, thread_accelerator->lattice[(*thread_elements)[element_nr]].fam_name.c_str());
//   pthread_mutex_unlock(thread_data->mutex);
//
// }


// static void thread_dynap_pya(ThreadSharedData* thread_data, int thread_id, long task_id) {
//
//   std::vector<DynApGridPoint>& grid = *thread_grid;
//   const std::vector<unsigned int>& elements = *thread_elements;
//
//   DynApGridPoint p = (*thread_grid)[task_id];
//   DynApGridPoint point;
//   unsigned int element_nr = task_id;
//   unsigned int start_element = elements[element_nr];
//
//
//   double py_init = (*thread_ma_e_init);
//   double py_delta = (*thread_ma_e_delta);
//   double nr_iterations = (*thread_ma_nr_iterations);
//   double nr_steps_back = (*thread_ma_nr_steps_back);
//   double rescale = (*thread_ma_rescale);
//
//   double py = py_init;
//   while (true) {
//     while (true) {
//       point.p = *thread_ma_p0;     // offset
//       point.p.py += py;            // sets trial energy deviation
//       point.start_element = start_element; point.lost_turn = 0; point.lost_element = start_element; point.lost_plane = Plane::no_plane;
//       std::vector<Pos<double> > new_pos;
//       Pos<double> p = point.p + (*thread_cod)[start_element];  // p initial condition for tracking
//       if (fabs(p.ry) < tiny_y_amp) p.ry = sgn(p.ry) * tiny_y_amp;
//       Status::type status = track_ringpass (*thread_accelerator, p, new_pos, thread_nr_turns, point.lost_turn, point.lost_element, point.lost_plane, false);
//       if (status != Status::success) { py -= py_delta; point.p.py = py; break; }
//       py += py_delta;
//     }
//     if (nr_iterations < 1) break; else nr_iterations--;
//     py -= py_delta * nr_steps_back;
//     py_delta *= rescale;
//     py += py_delta;
//   }
//
//   grid[task_id] = point;
//
//   pthread_mutex_lock(thread_data->mutex);
//   printf("thread:%02i|task:%06lu/%06lu  element:%04i|py:%+.4e  %s\n", thread_id, (1+task_id), thread_data->nr_tasks, element_nr, grid[task_id].p.py, thread_accelerator->lattice[(*thread_elements)[element_nr]].fam_name.c_str());
//   pthread_mutex_unlock(thread_data->mutex);
//
//
// }
