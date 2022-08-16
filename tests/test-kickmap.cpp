#include <trackcpp/lattice.h>
#include <trackcpp/optics.h>
#include <trackcpp/elements.h>
#include <trackcpp/tracking.h>
#include <algorithm>
#include <cmath>


void sirius_model(std::vector<Element>& the_ring) {

    //""" --- drift spaces --- """

    Element d10      = Element::drift("d10", 1.00000);

    //""" --- quadrupoles --- """
    Element qf = Element::quadrupole("qf",  0.500, +1.0);
	Element qd = Element::quadrupole("qd",  0.500, -1.0);

    Element epu = Element::kickmap("epu", "test_kicktable.txt", 1);

    the_ring.clear();
    the_ring.push_back(d10);
    the_ring.push_back(qf);
    the_ring.push_back(epu);
    the_ring.push_back(d10);
    the_ring.push_back(qd);

}


int main() {

    std::vector<Element> the_ring;

    Accelerator accelerator;
    Pos<double> fixed_point;
    Matrix m66;
    std::vector<Twiss> twiss;
    Twiss twiss0;
    Status::type status;

    accelerator.energy = 2.99792462e9;
    accelerator.harmonic_number = 864;
    accelerator.cavity_on = true;
    accelerator.radiation_on = RadiationState::damping;
    accelerator.vchamber_on = true;


    sirius_model(accelerator.lattice);

    // return 0;

    fixed_point.rx = 0.001;
    std::cout << fixed_point << std::endl;
    track_elementpass(accelerator.lattice[2], fixed_point, accelerator);
    std::cout << fixed_point << std::endl;
    return 0;

    status = calc_twiss(accelerator, fixed_point, m66, twiss, twiss0);
    std::cout << status << std::endl;


    Pos<double> fp = fixed_point;
    std::vector<Pos<double>> closed_orbit;
    Plane::type lost_plane;
    unsigned int element_offset = 0;

    status = track_linepass(accelerator, fp, closed_orbit, element_offset, lost_plane, true);
    if (status != Status::success) return status;

    std::vector<Matrix> atm;
    Pos<double> v0;
    status = track_findm66(accelerator, closed_orbit[0], atm, m66, v0);
    if (status != Status::success) return status;
    std::cout << m66[0][0] + m66[1][1] << std::endl;
    std::cout << m66[2][2] + m66[3][3] << std::endl;

    const double dpp = 1e-8;
    Pos<double> fpp = fixed_point;
    fpp.de += dpp;
    std::vector<Pos<double>> codp;

    twiss.clear(); twiss.reserve(atm.size());
    if (twiss0.isundef()) { // case of periodic solution
        twiss0.spos = 0.0;
        // --- beta functions
        double sin_mux = sgn(m66[0][1]) * std::sqrt(-m66[0][1]*m66[1][0]-pow(m66[0][0]-m66[1][1],2)/4.0);
        double sin_muy = sgn(m66[2][3]) * std::sqrt(-m66[2][3]*m66[3][2]-pow(m66[2][2]-m66[3][3],2)/4.0);
        twiss0.alphax = (m66[0][0]-m66[1][1])/2.0/sin_mux;
        twiss0.alphay = (m66[2][2]-m66[3][3])/2.0/sin_muy;
        twiss0.betax  =  m66[0][1]/sin_mux;
        twiss0.betay  =  m66[2][3]/sin_muy;
        // --- closed orbit
        twiss0.co = closed_orbit[0];

        // // --- dispersion function based on eta = (1 - M)^(-1) D
        // Vector Dx({m66[0][4], m66[1][4]});
        // Vector Dy({m66[2][4], m66[3][4]});
        // Matrix eye2({{1,0},{0,1}});
        // Matrix mx; m66.getM(mx, 2, 2, 0, 0); mx.linear_combination(1.0,eye2,-1.0,mx); mx.inverse();
        // Matrix my; m66.getM(my, 2, 2, 2, 2); my.linear_combination(1.0,eye2,-1.0,my); my.inverse();
        // twiss0.etax.multiplication(mx, Dx);
        // twiss0.etay.multiplication(my, Dy);

        // Dispersion Function based on tracking:
        Status::type status = track_findorbit4(accelerator, codp, fpp);
        std::cout << "h1" << std::endl;
        if (status != Status::success) return status;

        twiss0.etax[0] = (codp[0].rx - closed_orbit[0].rx) / dpp;
        twiss0.etax[1] = (codp[0].px - closed_orbit[0].px) / dpp;
        twiss0.etay[0] = (codp[0].ry - closed_orbit[0].ry) / dpp;
        twiss0.etay[1] = (codp[0].py - closed_orbit[0].py) / dpp;
    } else {
        fpp.rx += twiss0.etax[0] * dpp;
        fpp.px += twiss0.etax[1] * dpp;
        fpp.ry += twiss0.etay[0] * dpp;
        fpp.py += twiss0.etay[1] * dpp;
        Status::type status = track_linepass(accelerator, fpp, codp, element_offset, lost_plane, true);
        std::cout << "h2" << std::endl;
        if (status != Status::success) return status;
    }

    std::cout << "OK" << std::endl;

    // std::vector<Matrix> atm;
    // Pos<double> v0;
    // status = track_findm66(accelerator, closed_orbit[0], atm, m66, v0);
    // if (status != Status::success) return status;

    // Pos<double> fp;
    // std::vector<Pos<double>> closed_orbit, traj;
    // Plane::type lost_plane;
    // unsigned int element_offset = 0;
    // Status::type status1 = track_linepass(sirius, fp, traj, element_offset, lost_plane, true);
    // for(auto i=0; i <= closed_orbit)
    // std::cout << status1 << std::endl;

    // Status::type status1 = track_findm66(sirius, fp, tm, m66, v0); //, indices);
    // std::cout << status1 << std::endl;

    // Status::type status2 = calc_twiss(sirius, fp, m66, twiss, twiss0);
    // std::cout << status2 << std::endl;

    // std::cout << twiss.size() << std::endl;
    // std::cout << twiss[0].betax << std::endl;

    return 0;

}
