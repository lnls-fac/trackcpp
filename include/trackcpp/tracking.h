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
    Accelerator& accelerator, const Pos<double>& fixed_point,
    std::vector<Matrix>& tm, Matrix& m66, Pos<double>& v0);

Status::type track_findm66(
    Accelerator& accelerator, const Pos<double>& fixed_point,
    std::vector<Matrix>& tm, Matrix& m66, Pos<double>& v0,
    std::vector<unsigned int >& indices);

Status::type track_findorbit4(
    Accelerator& accelerator, std::vector<Pos<double> >& closed_orbit,
    const Pos<double>& fixed_point_guess = Pos<double>(0));

Status::type track_findorbit6(
    Accelerator& accelerator, std::vector<Pos<double> >& closed_orbit,
    const Pos<double>& fixed_point_guess = Pos<double>(0));

Pos<double>  linalg_solve4_posvec(
    const std::vector<Pos<double> >& M, const Pos<double>& b);

Pos<double>  linalg_solve6_posvec(
    const std::vector<Pos<double> >& M, const Pos<double>& b);


template <typename T>
Status::type track_elementpass (
    const Accelerator& accelerator,
    const Element& el, // element through which to track particle
    Pos<T> &orig_pos // initial electron coordinates
) {

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
    case PassMethod::pm_drift_g2l_pass:
        if ((status = pm_drift_g2l_pass<T>(orig_pos, el, accelerator)) != Status::success) return status;
        break;
    case PassMethod::pm_kickpoly_pass:
        if ((status = pm_kickpoly_pass<T>(orig_pos, el, accelerator)) != Status::success) return status;
        break;
    default:
        return Status::passmethod_not_defined;
    }
    return status;
}

template <typename T>
Status::type track_elementpass (
    const Accelerator& accelerator,
    const Element& el, // element through which to track particle
    std::vector<Pos<T> >& orig_pos // initial electron coordinates
){
    Status::type status = Status::success;
    for(auto&& pos: orig_pos) {
        Status::type status2 = track_elementpass(accelerator, el, pos);
        if (status2 != Status::success) status = status2;
    }
    return status;
}

// linepass
// --------
// tracks particles along a beam transport line
//
// inputs:
//     accelerator: Element vector representing the beam transport line.
//     orig_pos: Pos vector representing initial positions of particles.
//     indices: vector of integers defining at which elements to return.
//     element_offset: equivalent to shifting the lattice so that
//         '*element_offset' is the index for the first element.
//
// outputs:
//     element_offset: in case of problems with passmethods, such as particle
//         loss '*element_offset' is the index of the corresponding element.
//     pos: Pos vector of particles coordinates at the start of every element
//         selected by indices and at the end of the last element.
//     lost_plane: which plane the particle was lost.
//     lost_pos: return coordinates of lost particle.
//
// RETURN:
//      status do tracking (see 'auxiliary.h')
//
template <typename T>
Status::type track_linepass (
    const Accelerator& accelerator,
    Pos<T>& orig_pos,
    const std::vector<unsigned int>& indices,
    unsigned int& element_offset,
    std::vector<Pos<T> >& pos,
    Plane::type& lost_plane
) {

    Status::type status = Status::success;
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

        // syntactic-sugar for read-only access to element object parameters
        const Element& element = line[element_offset];

        // stores trajectory at entrance of each element
        if (indcs[i]) pos.push_back(orig_pos);

        status = track_elementpass(accelerator, element, orig_pos);
        lost_plane = check_particle_loss(accelerator, element, orig_pos);
        if (lost_plane != Plane::no_plane) status = Status::particle_lost;

        if (status != Status::success) {
            // fill rest of vector with nans
            for(int ii=i+1; ii<=nr_elements; ++ii) {
                if (indcs[ii]) pos.emplace_back(nan(""));
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

// linepass
// --------
// tracks particles along a beam transport line
//
// inputs:
//     accelerator: Element vector representing the beam transport line.
//     orig_pos: Pos vector representing initial positions of particles.
//     trajectory: whether to return coordinates at all elements.
//     element_offset: equivalent to shifting the lattice so that
//         '*element_offset' is the index for the first element.
//
// outputs:
//     element_offset: in case of problems with passmethods, such as particle
//         loss '*element_offset' is the index of the corresponding element.
//     pos: Pos vector of particles coordinates at the start of every element,
//         if trajectory == true, and at the end of the last element.
//     lost_plane: which plane the particle was lost.
//     lost_pos: return coordinates of lost particle.
//
// RETURN:
//      status do tracking (see 'auxiliary.h')
//
template <typename T>
Status::type track_linepass (
    const Accelerator& accelerator,
    Pos<T>& orig_pos,
    const bool trajectory,
    unsigned int& element_offset,
    std::vector<Pos<T> >& pos,
    Plane::type& lost_plane
) {
    std::vector<unsigned int> indices;
    unsigned int nr_elements = accelerator.lattice.size();
    if (trajectory){
        indices.reserve(nr_elements + 1);
        for (unsigned int i=0; i<=nr_elements; ++i) indices.push_back(i);
    }else{
        indices.push_back(nr_elements);
    }

    return track_linepass (
        accelerator,
        orig_pos,
        indices,
        element_offset,
        pos,
        lost_plane
    );
}


// linepass
// --------
// tracks particles along a beam transport line
//
// inputs:
//     accelerator: Element vector representing the beam transport line.
//     orig_pos: Pos vector representing initial positions of particles.
//     indices: vector of integers defining at which elements to return.
//     element_offset: equivalent to shifting the lattice so that
//         '*element_offset' is the index for the first element.
//
// outputs:
//     element_offset: in case of problems with passmethods, such as particle
//         loss '*element_offset' is the index of the corresponding element.
//     pos: Pos vector of particles coordinates at the start of every element
//         selected by indices and at the end of the last element.
//     lost_plane: which plane each particle was lost.
//     lost_flag: whether or not each particle was lost.
//     lost_element: which element each particle was lost.
//
// RETURN:
//      status do tracking (see 'auxiliary.h')
//
template <typename T>
Status::type track_linepass (
    const Accelerator& accelerator,
    std::vector<Pos<T>> &orig_pos,
    const std::vector<unsigned int >& indices,
    const unsigned int element_offset,
    std::vector<Pos<T>> &pos,
    std::vector<unsigned int >& lost_plane,
    std::vector<bool>& lost_flag,
    std::vector<int>& lost_element
) {

    int nr_elements = accelerator.lattice.size();
    Status::type status  = Status::success;
    Status::type status2  = Status::success;
    std::vector<Pos<T> > final_pos;
    Plane::type lp;

    pos.reserve(indices.size() * orig_pos.size());
    lost_flag.reserve(orig_pos.size());
    lost_plane.reserve(orig_pos.size());
    lost_element.reserve(orig_pos.size());

    for(unsigned int i=0; i<orig_pos.size(); ++i) {
        unsigned int le = element_offset;

        status2 = track_linepass(
            accelerator, orig_pos[i], indices, le, final_pos, lp
        );

        if (status2 != Status::success){
            status = status2;
            lost_element.push_back(le);
            lost_flag.push_back(true);
        } else {
            lost_element.push_back(-1);
            lost_flag.push_back(false);
        }
        lost_plane.push_back(lp);
        for (auto&& p: final_pos) pos.push_back(p);
        final_pos.clear();
    }
    return status;
}


// linepass
// --------
// tracks particles along a beam transport line
//
// inputs:
//     accelerator: Element vector representing the beam transport line.
//     orig_pos: Pos vector representing initial positions of particles.
//     indices: vector of integers defining at which elements to return.
//     element_offset: equivalent to shifting the lattice so that
//         '*element_offset' is the index for the first element.
//
// outputs:
//     element_offset: in case of problems with passmethods, such as particle
//         loss '*element_offset' is the index of the corresponding element.
//     pos: Pos vector of particles coordinates at the start of every element,
//         if trajectory == true, and at the end of the last element.
//     lost_plane: which plane each particle was lost.
//     lost_flag: whether or not each particle was lost.
//     lost_element: which element each particle was lost.
//
// RETURN:
//      status do tracking (see 'auxiliary.h')
//

// ringpass
// --------
// tracks particles around a ring
//
// inputs:
//     accelerator: Element vector representing the beam transport line.
//     orig_pos: Pos vector representing initial positions of particles.
//     nr_turns: number of turns for tracking
//     turn_by_turn: whether to return coordinates at all elements.
//     element_offset: equivalent to shifting the lattice so that
//         '*element_offset' is the index for the first element
//
// outputs:
//     element_offset: in case of problems with passmethods, such as particle
//         loss '*element_offset' is the index of the corresponding element.
//     pos: Pos vector of particles coordinates of the turn_by_turn data
//         (at end of the ring at each turn)
//     lost_plane: which plane each particle was lost.
//     lost_turn: which turn each particle was lost.
//
// RETURN:
//     status do tracking (see 'auxiliary.h')
template <typename T>
Status::type track_ringpass (
    const Accelerator& accelerator,
    Pos<T>& orig_pos,
    const unsigned int nr_turns,
    const bool turn_by_turn,
    unsigned int& element_offset,
    std::vector<Pos<T> >& pos,
    Plane::type& lost_plane,
    unsigned int& lost_turn
) {

    Status::type status  = Status::success;
    std::vector<Pos<T> > final_pos;

    if (turn_by_turn) pos.reserve(nr_turns+1);

    for(lost_turn=0; lost_turn<nr_turns; ++lost_turn) {

        // stores turn_by_turn at beggining of each turn
        if (turn_by_turn) pos.push_back(orig_pos);

        if ((status = track_linepass(
            accelerator,
            orig_pos,
            false,
            element_offset,
            final_pos,
            lost_plane
        )) != Status::success) {

            // fill last of vector with nans
            pos.emplace_back(nan(""));
            if (turn_by_turn) for(int i=lost_turn+1; i<nr_turns; ++i) {
                pos.emplace_back(nan(""));
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
    const unsigned int nr_turns,
    bool turn_by_turn,
    unsigned int element_offset,
    std::vector<Pos<T> > &pos,
    std::vector<unsigned int>& lost_plane,
    std::vector<int>& lost_turn,
    std::vector<bool>& lost_flag,
    std::vector<int>& lost_element
) {

    Status::type status  = Status::success;
    Status::type status2  = Status::success;
    std::vector<Pos<T> > final_pos;
    unsigned int lt;
    Plane::type lp;

    if (turn_by_turn) pos.reserve((nr_turns + 1) * orig_pos.size());
    else pos.reserve(orig_pos.size());
    lost_turn.reserve(orig_pos.size());
    lost_plane.reserve(orig_pos.size());
    lost_element.reserve(orig_pos.size());
    lost_flag.reserve(orig_pos.size());

    for(unsigned int i=0; i<orig_pos.size(); ++i) {
        unsigned int le = element_offset;

        status2 = track_ringpass(
            accelerator,
            orig_pos[i],
            nr_turns,
            turn_by_turn,
            le,
            final_pos,
            lp,
            lt
        );

        if (status2 != Status::success) status = status2;

        if (status2 != Status::success){
            status = status2;
            lost_turn.push_back(lt);
            lost_element.push_back(le);
            lost_flag.push_back(true);
        } else {
            lost_turn.push_back(-1);
            lost_element.push_back(-1);
            lost_flag.push_back(false);
        }
        lost_plane.push_back(lp);
        for (auto&& p: final_pos) pos.push_back(p);
        final_pos.clear();
    }
    return status;
}


// ------------------- auxiliary methods ----------------
template <typename T>
Plane::type check_particle_loss(
    const Accelerator& accelerator,
    const Element& ele,
    const Pos<T>& orig_pos
){

    Plane::type lost_plane = Plane::no_plane;
    const T& rx = orig_pos.rx;
    const T& ry = orig_pos.ry;

    if (not isfinite(rx)) lost_plane = Plane::x;
    if (not isfinite(ry)) {
        if (lost_plane == Plane::no_plane) return Plane::y;
        else return Plane::xy;
    }
    if (lost_plane != Plane::no_plane) return lost_plane;

    if (not accelerator.vchamber_on) return Plane::no_plane;

    // invalid p-norm shape (negative p). Safely signals lost particle
    if (ele.vchamber < 0) return Plane::xy;
    // rectangular vacuum chamber
    if (ele.vchamber == VChamberShape::rectangle) {
        if (((rx <= ele.hmin) or (rx >= ele.hmax))) lost_plane = Plane::x;
        if (((ry <= ele.vmin) or (ry >= ele.vmax))) {
            if (lost_plane == Plane::no_plane) return Plane::y;
            else return Plane::xy;
        }
        if (lost_plane != Plane::no_plane) return lost_plane;
    }
    // lost in rhombus, elliptic and all finite p-norm shapes
    else if (get_norm_amp_in_vchamber(ele, rx, ry) > 1) return Plane::xy;
    return lost_plane;
}


template <typename T>
T get_norm_amp_in_vchamber(const Element& elem, const T& rx, const T& ry) {
    double lx = (elem.hmax - elem.hmin) / 2;
    double ly = (elem.vmax - elem.vmin) / 2;
    double xc = (elem.hmax + elem.hmin) / 2;
    double yc = (elem.vmax + elem.vmin) / 2;
    T xn = abs((rx - xc)/lx);
    T yn = abs((ry - yc)/ly);
    T amplitude;
    if (elem.vchamber == VChamberShape::rhombus) {
        amplitude = xn + yn;
    } else if (elem.vchamber == VChamberShape::ellipse) {
        amplitude = xn*xn + yn*yn;
    } else {
        amplitude = pow(xn, elem.vchamber) + pow(yn, elem.vchamber);
    }
    return amplitude;
}

#endif
