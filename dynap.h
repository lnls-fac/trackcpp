#ifndef _DYNAP_H
#define _DYNAP_H

#include "accelerator.h"
#include "elements.h"
#include "pos.h"
#include "auxiliary.h"
#include <vector>


struct DynApGridPoint {
	Pos<double>  p;              // Dislocation around closed-orbit
	unsigned int start_element;
	unsigned int lost_turn;
	unsigned int lost_element;
	Plane::type  lost_plane;     // Plane::no_plane,Plane::x,Plane::y,Plane::z
	double       nux1, nuy1;     // tunes at first half number of turns
	double       nux2, nuy2;     // tunes at second half number of turns
};

Status::type dynap_xy(
		const Accelerator& accelerator,
		std::vector<Pos<double> >& cod,
		unsigned int nr_turns,
		const Pos<double>& p0,
		unsigned int nrpts_x, double x_min, double x_max,
		unsigned int nrpts_y, double y_min, double y_max,
		bool calculate_closed_orbit,
		std::vector<DynApGridPoint>& grid
	);

Status::type dynap_ex(
		const Accelerator& accelerator,
		std::vector<Pos<double> >& cod,
		unsigned int nr_turns,
		const Pos<double>& p0,
		unsigned int nrpts_e, double e_min, double e_max,
		unsigned int nrpts_x, double x_min, double x_max,
		bool calculate_closed_orbit,
		std::vector<DynApGridPoint>& grid
	);

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
	);

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
	);

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
	);


#endif
