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

#ifndef _TRACKING_H
#define _TRACKING_H

#include "passmethods.h"
#include "accelerator.h"
#include "elements.h"
#include "auxiliary.h"
#include "linalg.h"
#include <vector>
#include <limits>
#include <utility>
#include <cmath>

Status::type track_findm66(
	const Accelerator& accelerator, const Pos<double>& fixed_point,
	std::vector<Matrix>& tm, Matrix& m66, Pos<double>& v0);

Status::type track_findm66(
	const Accelerator& accelerator, const Pos<double>& fixed_point,
	std::vector<Matrix>& tm, Matrix& m66, Pos<double>& v0,
	std::vector<unsigned int >& indices);

Status::type track_findorbit4(
	const Accelerator& accelerator, std::vector<Pos<double> >& closed_orbit,
	const Pos<double>& fixed_point_guess = Pos<double>(0));

Status::type track_findorbit6(
	const Accelerator& accelerator, std::vector<Pos<double> >& closed_orbit,
	const Pos<double>& fixed_point_guess = Pos<double>(0));

Pos<double>  linalg_solve4_posvec(
	const std::vector<Pos<double> >& M, const Pos<double>& b);

Pos<double>  linalg_solve6_posvec(
	const std::vector<Pos<double> >& M, const Pos<double>& b);


template <typename T>
Status::type track_elementpass (
		     const Element& el,                 // element through which to track particle
		     Pos<T> &orig_pos,                  // initial electron coordinates
		     const Accelerator& accelerator) {

	Status::type status = Status::success;

	switch (el.pass_method) {
	case PassMethod::pm_identity_pass:
		if ((status = pm_identity_pass<T>(orig_pos, el, accelerator)) != Status::success) return status;
		break;
	case PassMethod::pm_drift_pass:
		if ((status = pm_drift_pass<T>(orig_pos, el, accelerator)) != Status::success) return status;
		break;
	case PassMethod::pm_str_mpole_symplectic4_pass:
		if ((status = pm_str_mpole_symplectic4_pass<T>(orig_pos, el, accelerator)) != Status::success) return status;
		break;
	case PassMethod::pm_bnd_mpole_symplectic4_pass:
		if ((status = pm_bnd_mpole_symplectic4_pass<T>(orig_pos, el, accelerator)) != Status::success) return status;
		break;
	case PassMethod::pm_corrector_pass:
		if ((status = pm_corrector_pass<T>(orig_pos, el, accelerator)) != Status::success) return status;
		break;
	case PassMethod::pm_cavity_pass:
		if ((status = pm_cavity_pass<T>(orig_pos, el, accelerator)) != Status::success) return status;
		break;
	case PassMethod::pm_thinquad_pass:
		if ((status = pm_thinquad_pass<T>(orig_pos, el, accelerator)) != Status::success) return status;
		break;
	case PassMethod::pm_thinsext_pass:
		if ((status = pm_thinsext_pass<T>(orig_pos, el, accelerator)) != Status::success) return status;
		break;
	case PassMethod::pm_kickmap_pass:
		if ((status = pm_kickmap_pass<T>(orig_pos, el, accelerator)) != Status::success) return status;
		break;
	case PassMethod::pm_matrix_pass:
		if ((status = pm_matrix_pass<T>(orig_pos, el, accelerator)) != Status::success) return status;
		break;
	default:
		return Status::passmethod_not_defined;
	}

	return status;

}

template <typename T>
Status::type track_elementpass (
		     const Element& el,  // element through which to track particle
		     std::vector<Pos<T> >& orig_pos,  // initial electron coordinates
		     const Accelerator& accelerator) {

	Status::type status  = Status::success;

	for(auto&& pos: orig_pos) {
		Status::type status2 = track_elementpass(el, pos, accelerator);
		if (status2 != Status::success) status = status2;
	}
	return status;
}



// linepass
// --------
// tracks particles along a beam transport line
//
// inputs:
//		line:	 		Element vector representing the beam transport line
//		orig_pos:		Pos vector representing initial positions of particles
//		element_offset:	equivalent to shifting the lattice so that '*element_offset' is the index for the first element
//		trajectory:		flag indicating that trajectory is to be recorded at entrance of all elements
//						(otherwise only the coordinates at the exit of last element is recorded)
// outputs:
//		pos:			Pos vector of particles' final coordinates (or trajetory)
//		element_offset:	in case of problems with passmethods, '*element_offset' is the index of the corresponding element
//		RETURN:			status do tracking (see 'auxiliary.h')

template <typename T>
Status::type track_linepass (
		const Accelerator& accelerator,
		Pos<T>& orig_pos,              // initial electron coordinates
		std::vector<Pos<T> >& pos,     // vector with tracked electron coordinates at start of every element and at the end of last one.
		unsigned int& element_offset,  // index of starting element for tracking
		Plane::type& lost_plane,       // return plane in which particle was lost, if the case.
		std::vector<unsigned int >& indices) {// indices to return;

	Status::type status = Status::success;
	lost_plane = Plane::no_plane;

	const std::vector<Element>& line = accelerator.lattice;
	int nr_elements  = line.size();

	//pos.clear(); other functions assume pos is not clearedin linepass!
	pos.reserve(pos.size() + indices.size());

	// create vector of booleans help determine when to store position
	std::vector<bool> indcs;
	indcs.reserve(nr_elements+1);
	for (unsigned int i=0; i<=nr_elements; ++i) indcs[i] = false;
	for (auto&& i: indices) if (i<=nr_elements) indcs[i] = true;

	for(int i=0; i<nr_elements; ++i) {

		const Element& element = line[element_offset];  // syntactic-sugar for read-only access to element object parameters

		// stores trajectory at entrance of each element
		if (indcs[i]) pos.push_back(orig_pos);

		status = track_elementpass (element, orig_pos, accelerator);

		const T& rx = orig_pos.rx;
		const T& ry = orig_pos.ry;

		// checks if particle is lost

		if (not isfinite(rx)) {
			lost_plane = Plane::x;
			status = Status::particle_lost;
		}
		if (not isfinite(ry)) {
			if (status != Status::particle_lost) {
				lost_plane = Plane::y;
				status = Status::particle_lost;
			} else {
				lost_plane = Plane::xy;
			}
		}
		if ((status != Status::particle_lost) and accelerator.vchamber_on) {
			if (element.vchamber == VChamberShape::rectangle) {
				// rectangular vaccum chamber
				if (((rx < element.hmin) or (rx > element.hmax))) {
					lost_plane = Plane::x;
					status = Status::particle_lost;
				}
				if (((ry < element.vmin) or (ry > element.vmax))) {
					if (status != Status::particle_lost) {
						lost_plane = Plane::y;
						status = Status::particle_lost;
					} else {
						lost_plane = Plane::xy;
					}
				}
			} else if (element.vchamber == VChamberShape::kite) {
				// kite-shaped vaccum chamber
				if (((rx < element.hmin) or (rx > element.hmax))) {
					lost_plane = Plane::xy;
					status = Status::particle_lost;
				}
				if ((ry > get_kite_ry(element, 0, rx)) or  // upper right
				    (ry > get_kite_ry(element, 1, rx)) or  // upper left
					(ry < get_kite_ry(element, 2, rx)) or  // lower left
					(ry < get_kite_ry(element, 3, rx))) {  // lower right
					lost_plane = Plane::xy;
					status = Status::particle_lost;
				}
			} else if (element.vchamber == VChamberShape::ellipse) {
				// elliptical vaccum chamber
				double lx = (element.hmax - element.hmin) / 2;
				double ly = (element.vmax - element.vmin) / 2;
				double xc = (element.hmax + element.hmin) / 2;
				double yc = (element.vmax + element.vmin) / 2;
				T xn = (rx - xc)/lx;
				T yn = (ry - yc)/ly;
				T amplitude = xn*xn + yn*yn;
				if (amplitude > 1) {
					lost_plane = Plane::xy;
					status = Status::particle_lost;
				}
			} else {
				// any other shape not implemented safely signals lost particle
				lost_plane = Plane::xy;
				status = Status::particle_lost;
			}
		}

		if (status != Status::success) {
			// fill rest of vector with nans
			for(int ii=i+1; ii<=nr_elements; ++ii) {
				if (indcs[ii]) pos.emplace_back(
					nan(""),nan(""),nan(""),nan(""),nan(""),nan(""));
			}
			return status;
		}
		// moves to next element index
		element_offset = (element_offset + 1) % nr_elements;
	}

	// stores final particle position at the end of the line
	if (indcs[nr_elements]) pos.push_back(orig_pos);

	return (status == Status::success) ? status: Status::particle_lost;
}


template <typename T>
T get_kite_ry(const Element& elem, int side, const T& rx) {
	double x1, y1, x2, y2;
	if (side == 0) {
		// upper-right
		x1 = 0; y1 = elem.vmax;
		x2 = elem.hmax; y2 = 0;
	} else if (side == 1) {
		// upper-left
		x1 = elem.hmin; y1 = 0;
		x2 = 0; y1 = elem.vmax;
	} else if (side == 2) {
		// lower-left
		x1 = elem.hmin; y1 = 0;
		x2 = 0; y1 = elem.vmin;
	} else {
		// lower-right
		x1 = 0; y1 = elem.vmin;
		x2 = elem.hmax; y1 = 0;
	}
	return y1 + (rx - x1) * (y2 - y1) / (x2 - x1);
}


template <typename T>
Status::type track_linepass (
		const Accelerator& accelerator,
		Pos<T>& orig_pos,              // initial electron coordinates
		std::vector<Pos<T> >& pos,     // vector with tracked electron coordinates at start of every element and at the end of last one.
		unsigned int& element_offset,  // index of starting element for tracking
		Plane::type& lost_plane,       // return plane in which particle was lost, if the case.
		bool trajectory) {             // whether function should return coordinates at all elements

	std::vector<unsigned int> indices;
	unsigned int nr_elements = accelerator.lattice.size();

	if (trajectory){
		indices.reserve(nr_elements + 1);
		for (unsigned int i=0; i<=nr_elements; ++i) indices.push_back(i);
	}else{
		indices.push_back(nr_elements);
	}

	return track_linepass (
		accelerator, orig_pos, pos, element_offset, lost_plane, indices);
}


template <typename T>
Status::type track_linepass (
		const Accelerator& accelerator,
		std::vector<Pos<T> > &orig_pos,
		std::vector<Pos<T> > &pos,
		unsigned int element_offset,
		std::vector<unsigned int >& lost_plane,
		std::vector<unsigned int >& lost_element,
		std::vector<unsigned int >& indices) {

	int nr_elements = accelerator.lattice.size();
	Status::type status  = Status::success;
	Status::type status2  = Status::success;
	std::vector<Pos<T> > final_pos;
	Plane::type lp;

	pos.reserve(indices.size() * orig_pos.size());
    lost_plane.reserve(orig_pos.size());
	lost_element.reserve(orig_pos.size());

	for(unsigned int i=0; i<orig_pos.size(); ++i) {
		unsigned int le = element_offset;

		status2 = track_linepass (
			accelerator, orig_pos[i], final_pos, le, lp, indices);

		if (status2 != Status::success) status = status2;

		lost_plane.push_back(lp);
		lost_element.push_back(le);
		for (auto&& p: final_pos) pos.push_back(p);
		final_pos.clear();
	}
	return status;
}


// ringpass
// --------
// tracks particles around a ring
//
// inputs:
//		ring: 			Element vector representing the ring
//		orig_pos:		Pos vector representing initial positions of particles
//		element_offset:	equivalent to shifting the lattice so that '*element_offset' is the index for the first element
//		nr_turns:		number of turns for tracking
//
// outputs:
//		pos:			Pos vector of particles' coordinates of the turn_by_turn data (at end of the ring at each turn)
//		turn_idx:		in case of problems with passmethods, '*turn_idx' is the index of the corresponding turn
//		element_offset:	in case of problems with passmethods, '*element_offset' is the index of the corresponding element
//		RETURN:			status do tracking (see 'auxiliary.h')


template <typename T>
Status::type track_ringpass (
		const Accelerator& accelerator,
		Pos<T> &orig_pos,
		std::vector<Pos<T> >& pos,
		const unsigned int nr_turns,
		unsigned int& lost_turn,
		unsigned int& element_offset,
		Plane::type& lost_plane,
		bool trajectory) {

	Status::type status  = Status::success;
	std::vector<Pos<T> > final_pos;

	if (trajectory) pos.reserve(nr_turns+1);

	for(lost_turn=0; lost_turn<nr_turns; ++lost_turn) {

		// stores trajectory at beggining of each turn
		if (trajectory) pos.push_back(orig_pos);

		if ((status = track_linepass (accelerator, orig_pos, final_pos, element_offset, lost_plane, false)) != Status::success) {

			// fill last of vector with nans
			pos.emplace_back(
				nan(""),nan(""),nan(""),nan(""),nan(""),nan(""));
			if (trajectory) for(int i=lost_turn+1; i<nr_turns; ++i) {
					pos.emplace_back(
						nan(""),nan(""),nan(""),nan(""),nan(""),nan(""));
				}
			return status;
		}
		final_pos.clear();
	}
	pos.push_back(orig_pos);

	return status;
}


template <typename T>
Status::type track_ringpass (
		const Accelerator& accelerator,
		std::vector<Pos<T> > &orig_pos,
		std::vector<Pos<T> > &pos,
		const unsigned int nr_turns,
		std::vector<unsigned int > &lost_turn,
		unsigned int element_offset,
		std::vector<unsigned int >& lost_plane,
		std::vector<unsigned int >& lost_element,
		bool trajectory) {

	Status::type status  = Status::success;
	Status::type status2  = Status::success;
	std::vector<Pos<T> > final_pos;
	unsigned int lt;
	Plane::type lp;

	if (trajectory) pos.reserve((nr_turns + 1) * orig_pos.size());
    else pos.reserve(orig_pos.size());
	lost_turn.reserve(orig_pos.size());
	lost_plane.reserve(orig_pos.size());
	lost_element.reserve(orig_pos.size());

	for(unsigned int i=0; i<orig_pos.size(); ++i) {
		unsigned int le = element_offset;

		status2 = track_ringpass(
			accelerator, orig_pos[i], final_pos, nr_turns,
			lt, le, lp, trajectory);

		if (status2 != Status::success) status = status2;

		lost_turn.push_back(lt);
		lost_plane.push_back(lp);
		lost_element.push_back(le);
		for (auto&& p: final_pos) pos.push_back(p);
		final_pos.clear();
	}
	return status;
}

#endif
