#include "trackc++.h"
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

static Status::type calc_closed_orbit(const Accelerator& accelerator, std::vector<Pos<double> >& cod, const char* function_name) {
	if (verbose_on) {
		std::cout << get_timestamp() << " <" << function_name << ">: calculating closed-orbit...";
		std::cout.flush();
	}
	Status::type status = track_findorbit6(accelerator, cod);
	if (verbose_on) std::cout << string_error_messages[status] <<  std::endl;
	return Status::success;
}

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


	// finds 6D closed-orbit
	if (calculate_closed_orbit) {
		Status::type status = calc_closed_orbit(accelerator, cod, __FUNCTION__);
		if (status != Status::success) return status;
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
			Status::type status = track_ringpass (accelerator, p, new_pos, nr_turns, grid[idx].lost_turn, grid[idx].lost_element, grid[idx].lost_plane, false);

			if (verbose_on) {
				sprintf(buffer, "turn=%05i element=%05i  %20s", grid[idx].lost_turn, grid[idx].lost_element, string_error_messages[status].c_str());
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

	// finds 6D closed-orbit
	if (calculate_closed_orbit) {
		Status::type status = calc_closed_orbit(accelerator, cod, __FUNCTION__);
		if (status != Status::success) return status;
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
			Status::type status = track_ringpass (accelerator, p, new_pos, nr_turns, grid[idx].lost_turn, grid[idx].lost_element, grid[idx].lost_plane, false);

			if (verbose_on) {
				sprintf(buffer, "turn=%05i element=%05i  %20s", grid[idx].lost_turn, grid[idx].lost_element, string_error_messages[status].c_str());
				std::cout << buffer << std::endl;
			}

			idx++;
		}
	}

	if (verbose_on) std::cout << get_timestamp() << " end of dynap_ex loop" << std::endl;

	return Status::success;

}


DynApGridPoint find_fine_momentum_acceptance(
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

DynApGridPoint find_momentum_acceptance(
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

Status::type dynap_ma_track(
		const Accelerator& accelerator,
		std::vector<Pos<double> >& cod,
		unsigned int nr_turns,
		const Pos<double>& p0,
		const double& energy_offset,
		unsigned int element_idx
	) {
	Pos<double> p = p0;     // offset
	p.de += energy_offset;
	//unsigned int start_element = element_idx;
	unsigned int lost_turn = 0;
	unsigned int lost_element = element_idx;
	Plane::type lost_plane = Plane::no_plane;
	std::vector<Pos<double> > new_pos;
	Pos<double> pos = p + cod[element_idx];  // p initial condition for tracking
	if (fabs(pos.ry) < tiny_y_amp) pos.ry = sgn(pos.ry) * tiny_y_amp;
	return track_ringpass (accelerator, pos, new_pos, nr_turns, lost_turn, lost_element, lost_plane, false);

}

Status::type dynap_ma2(
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

	const std::vector<Element>& the_ring = accelerator.lattice;

	// finds 6D closed-orbit
	if (calculate_closed_orbit) {
		Status::type status = calc_closed_orbit(accelerator, cod, __FUNCTION__);
		if (status != Status::success) return status;
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



	for(unsigned int i = 0; i<elements.size(); ++i) {

		auto element_idx = i;

		unsigned int this_nr_turns = (nr_turns > 16) ? 16 : nr_turns;

		// goes up with small number of turns
		DynApGridPoint point;
		double energy_offset = 0;
		while (true) {
			point.p = p0;     // offset
			point.p.de += energy_offset;
			point.start_element = element_idx; point.lost_turn = 0; point.lost_element = element_idx; point.lost_plane = Plane::no_plane;
			std::vector<Pos<double> > new_pos;
			Pos<double> p = point.p + cod[element_idx];  // p initial condition for tracking
			if (fabs(p.ry) < tiny_y_amp) p.ry = sgn(p.ry) * tiny_y_amp;
			Status::type status = track_ringpass (accelerator, p, new_pos, this_nr_turns, point.lost_turn, point.lost_element, point.lost_plane, false);
			if (status == Status::success) {
				energy_offset += e_tol;
			} else {
				energy_offset = energy_offset - e_tol;
				break;
			}
		}

		// goes down increasing number of turns
		this_nr_turns = (2 * this_nr_turns > nr_turns) ? nr_turns : 2 * this_nr_turns;


		while (true) {

			this_nr_turns = (2 * this_nr_turns > nr_turns) ? nr_turns : 2 * this_nr_turns;

		}

	}

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

	const std::vector<Element>& the_ring = accelerator.lattice;

	// finds 6D closed-orbit
	if (calculate_closed_orbit) {
		Status::type status = calc_closed_orbit(accelerator, cod, __FUNCTION__);
		if (status != Status::success) return status;
	}

	// finds out which in which elements tracking is to be performed
	std::vector<unsigned int> elements;
	double s = 0.0;
	for(unsigned int i=0; i<the_ring.size(); ++i) {
		if ((s >= s_min) and (s <= s_max)) {
			if (std::find(fam_names.begin(), fam_names.end(), the_ring[i].fam_name) != fam_names.end()) {
				elements.push_back(i);                      // calcs at start of element
//#pragma message("temporary: elements.push_back(i+1) !!!!")
//				elements.push_back((i+1) % the_ring.size()); // calcs at end of element
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
		//grid[idx++] = find_momentum_acceptance(accelerator, cod, nr_turns, p0, -1.0 * e0, e_tol, elements[i]);
		if (verbose_on) {
			fprintf(stdout, "%+11.4E ", grid[idx-1].p.de); fflush(stdout);
		}
		// finds positive energy acceptance
		p = find_momentum_acceptance(accelerator, cod, nr_turns, p0, +1.0 * e0, 5 * e_tol, elements[i]);
		grid[idx++] = find_fine_momentum_acceptance(accelerator, cod, nr_turns, p0, p.p.de, e_tol, elements[i]);
		//grid[idx++] = find_momentum_acceptance(accelerator, cod, nr_turns, p0, +1.0 * e0, e_tol, elements[i]);
		if (verbose_on) {
			fprintf(stdout, "%+11.4E ", grid[idx-1].p.de); fflush(stdout);
		}
		std::cout << std::endl;

	}

	return Status::success;

}
