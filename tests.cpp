// TRACKC++
// ========
// Author: 		Ximenes R. Resende
// email:  		xresende@gmail.com, ximenes.resende@lnls.br
// affiliation:	LNLS - Laboratorio Nacional de Luz Sincrotron
// Date: 		Tue Dec 10 17:57:20 BRST 2013

#include "trackc++.h"
#include <ctime>

int test_printlattice(const Accelerator& accelerator) {
	latt_print(accelerator.lattice);
	return 0;
}

int test_linepass(const Accelerator& accelerator) {

	const std::vector<Element>& the_ring = accelerator.lattice;

	Pos<> pos;
	pos.rx = 1*0.00100; pos.px = 0*0.00001;
	pos.ry = 1*0.00010; pos.py = 0*0.00001;

	std::vector<Pos<> > new_pos;
	unsigned int element_offset = 0;
	Plane::type lost_plane;
	Status::type status = track_linepass(accelerator, pos, new_pos, element_offset, lost_plane, true);
	std::cout << "status: " << string_error_messages[status] << std::endl;


	FILE *fp;
	fp = fopen("orbit_trackc++.txt", "w");
	for(unsigned int i=1-1; i<new_pos.size(); ++i) {
	//for(unsigned int i=2000; i<292; ++i) {
		const Pos<>& c = new_pos[i];
		fprintf(stdout, "%03i: %15s  %+23.16E %+23.16E %+23.16E %+23.16E %+23.16E %+23.16E \n", i+1, the_ring[i % the_ring.size()].fam_name.c_str(), c.rx, c.px, c.ry, c.py, c.de, c.dl);
		fprintf(fp, "%+23.16E %+23.16E %+23.16E %+23.16E %+23.16E %+23.16E \n", c.rx, c.px, c.ry, c.py, c.de, c.dl);
	}
	fclose(fp);

	return 0;

}


int test_linepass_tpsa(const Accelerator& accelerator, const std::vector<Element>& the_ring) {

	const int order = 3;
	Pos<Tpsa<6,order> > tpsa;

	tpsa.rx = Tpsa<6,order>(0, 0); tpsa.px = Tpsa<6,order>(0, 1);
	tpsa.ry = Tpsa<6,order>(0, 2); tpsa.py = Tpsa<6,order>(0, 3);
	tpsa.de = Tpsa<6,order>(0, 4); tpsa.dl = Tpsa<6,order>(0, 5);
	std::vector<Pos<Tpsa<6,order> > > new_tpsa;
	unsigned int element_offset = 0;
	Plane::type lost_plane;
	track_linepass(accelerator, tpsa, new_tpsa, element_offset, lost_plane, false);
	for(unsigned int i=0; i<new_tpsa.size(); ++i) {
		//const Pos<Tpsa<6,1> >& c = new_particles[i];
		//std::cout << c.rx << std::endl;
		//fprintf(stdout, "%03i: %15s  %+23.16E %+23.16E %+23.16E %+23.16E %+23.16E\n", i+1, the_ring[i].fam_name.c_str(), c.rx, c.px, c.ry, c.py, c.de);
	}

	return 0;

}


#include <cstdlib>
int test_ringpass(const Accelerator& accelerator) {


	Pos<> pos;
	pos.rx = 0.00100;
	pos.ry = 0.00010;

	std::vector<Pos<double> > new_pos;
	unsigned int element_offset = 0, lost_turn = 0;
	Plane::type lost_plane;

	clock_t begin, end;
	double time_spent;
	begin = clock();
	Status::type status = track_ringpass(accelerator, pos, new_pos, 5000, lost_turn, element_offset, lost_plane, true);
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

	if (status != Status::success) {
		std::cerr << "problem" << std::endl;
	}

	std::cout << "tracking_time: " << time_spent << std::endl;
	for(unsigned int i=0; i<new_pos.size(); ++i) {
		const Pos<>& c = new_pos[i];
		fprintf(stdout, "%03i: %+23.16E %+23.16E %+23.16E %+23.16E %+23.16E %+23.16E\n", i+1, c.rx, c.px, c.ry, c.py, c.de, c.dl);
	}


	return 0;

}

int test_findorbit6(const Accelerator& accelerator) {

	std::vector<Pos<double> > orbit;
	Status::type status = track_findorbit6(accelerator, orbit);
	if (status == Status::findorbit_not_converged) {
		std::cerr << "findorbit not converged!" << std::endl;
	} else {
		const Pos<>& c = orbit[0];
		fprintf(stdout, "closed_orbit: %+23.16E %+23.16E %+23.16E %+23.16E %+23.16E\n", c.rx, c.px, c.ry, c.py, c.de);
	}
	return 0;
}

// #include <cstdio>
// int test_findm66(const Accelerator& accelerator) {
//
// 	const std::vector<Element>& the_ring = accelerator.lattice;
//
// 	std::vector<Matrix> m66;
//
// 	track_findm66 (accelerator, cod, m66);
//
// 	for(unsigned int i=0; i<the_ring.size(); ++i) {
// 		std::cout << "element#     : " << i+1 << std::endl;
// 		std::cout << the_ring[i];
// 		for(unsigned int r=0; r<6; ++r) {
// 			for(unsigned int c=0; c<6; ++c) {
// 				printf("%+10.4E ", m66[i][r][c]);
// 			}
// 			std::cout << std::endl;
// 		}
// 		std::cout << std::endl;
// 	}
// 	return 0;
// }

int test_dynap_xy(const Accelerator& accelerator) {

	std::vector<Pos<double> > closed_orbit;
	unsigned int nr_turns = 5000;
	Pos<double> p0(0,0,0,0,0,0);
	unsigned int nrpts_x = 10;
	double x_min = -0.015, x_max = +0.015;
	unsigned int nrpts_y = 10;
	double y_min = 0, y_max = +0.0035;
	std::vector<DynApGridPoint> points;
	dynap_xy(accelerator, closed_orbit, nr_turns, p0, nrpts_x, x_min, x_max, nrpts_y, y_min, y_max, true, points);

	return 0;
}

int test_cmd_track_linepass() {

	std::vector<std::string> args = {
			"trackc++",
			"track_linepass",
			"/home/ximenes/pytrack/sirius_v500_ac10_5_bare_in_TRACY.txt",
			"3e9",
			"864",
			"off",
			"off",
			"off",
			"0",
			"0.001",
			"0.0",
			"0.0",
			"0.0",
			"0.0",
			"0.0"
	};

	return cmd_track_linepass(args);

}

int test_cmd_dynap_xy() {

	std::vector<std::string> args = {
			"trackc++",
			"dynap_xy",
			"/home/ximenes/pytrack/sirius_v500_ac10_5_bare_with_ids_in.txt",
			"3e9",
			"864",
			"on",
			"on",
			"on",
			"0.0",
			"5000",
			"4",
			"-0.015",
			"+0.015",
			"4",
			"0.0",
			"+0.0035"
	};
	return cmd_dynap_xy(args);

}


int test_cmd_dynap_ex() {

	std::vector<std::string> args = {
			"trackc++",
			"dynap_ex",
			"/home/ximenes/pytrack/sirius_v500_ac10_5_bare_in.txt",
			"3e9",
			"864",
			"on",
			"on",
			"on",
			"1e-6",
			"5000",
			"2",
			"-0.05",
			"+0.05",
			"2",
			"-0.015",
			"+0.015"
	};
	return cmd_dynap_xy(args);

}


int test_cmd_dynap_ma() {

	std::vector<std::string> args = {
			"trackc++",
			"dynap_ma",
			"/home/ximenes/pytrack/sirius_v500_ac10_5_bare_in.txt",
			"3e9",
			"864",
			"on",
			"on",
			"off",
			"5000",
			"30e-6",   // y0
			"0.01",    // e0
			"0.0005",  // tol_e
			"0",       // s_min [m]
			"30",      // s_max [m]
			"qaf",
			"qad"};

	return cmd_dynap_ma(args);

}

int test_kicktable(Accelerator& accelerator) {

	Kicktable t;
	const Kicktable *ptrKicktable = nullptr;
	add_kicktable("/home/fac_files/code/python/trackc++/pytrack/id_kicktable.txt", accelerator.kicktables, ptrKicktable);
	add_kicktable("/home/fac_files/code/python/trackc++/pytrack/id_kicktable2.txt", accelerator.kicktables, ptrKicktable);
	return 0;

}

int test_read_flat_file(Accelerator& accelerator) {

	read_flat_file("/home/ximenes/flatfile.txt", accelerator);
	std::cout << "nr_elements: " << accelerator.lattice.size() << std::endl;
	return 0;

}

int test_simple_drift() {

	Accelerator accelerator;

	accelerator.energy = 3e9; // [ev]
	accelerator.harmonic_number = 864;
	accelerator.radiation_on = false;
	accelerator.cavity_on = false;
	accelerator.vchamber_on = false;

	Element ds = Element::drift("ds", 1.0);

	accelerator.lattice.push_back(ds);

	Pos<double> pos(0.001,0.002,0.003,0.004,0.005,0.006);
	track_elementpass (ds, pos, accelerator);

	fprintf(stdout, "test_simple_drift\n");
	fprintf(stdout, "rx: %+.16f\n", pos.rx);
	fprintf(stdout, "px: %+.16f\n", pos.px);
	fprintf(stdout, "ry: %+.16f\n", pos.ry);
	fprintf(stdout, "py: %+.16f\n", pos.py);
	fprintf(stdout, "de: %+.16f\n", pos.de);
	fprintf(stdout, "dl: %+.16f\n", pos.dl);

}

int test_simple_quadrupole() {

	Accelerator accelerator;

	accelerator.energy = 3e9; // [ev]
	accelerator.harmonic_number = 864;
	accelerator.radiation_on = false;
	accelerator.cavity_on = false;
	accelerator.vchamber_on = false;

	Element ds = Element::quadrupole("qs", 1.0, 2.0);

	accelerator.lattice.push_back(ds);

	Pos<double> pos(0.001,0.002,0.003,0.004,0.005,0.006);
	track_elementpass (ds, pos, accelerator);

	fprintf(stdout, "test_simple_quadrupole\n");
	fprintf(stdout, "rx: %+.16f\n", pos.rx);
	fprintf(stdout, "px: %+.16f\n", pos.px);
	fprintf(stdout, "ry: %+.16f\n", pos.ry);
	fprintf(stdout, "py: %+.16f\n", pos.py);
	fprintf(stdout, "de: %+.16f\n", pos.de);
	fprintf(stdout, "dl: %+.16f\n", pos.dl);

}

int cmd_tests(const std::vector<std::string>& args) {


	//Accelerator accelerator;
	//sirius_v500(accelerator.lattice);
	//latt_read_flat_file("/home/ximenes/pytrack/sirius_v500_ac10_5_bare_in.txt", accelerator);
	//Status::type status = latt_read_flat_file("/home/fac_files/code/python/trackc++/pytrack/flat_file_ff.txt", accelerator);
	//if (status != Status::success) {
	//	return EXIT_FAILURE;
	//}
	//accelerator.lattice[15].nr_steps = 1;
	// accelerator.energy = 3e9; // [ev]
	// accelerator.harmonic_number = 864;
	// accelerator.radiation_on = true;
	// accelerator.cavity_on = true;
	// accelerator.vchamber_on = false;
	//accelerator.radiation_on = false;
	//accelerator.cavity_on = false;

	//latt_setcavity(the_ring, "on");
	//latt_setradiation(the_ring, "on", 3e9);
	//the_ring[13].hkick = 1e-4;

	//std::vector<Element> the_ring2(the_ring);
	//the_ring2.insert(the_ring2.begin(), the_ring.begin(), the_ring.end());
	//latt_print(the_ring);
	//std::cout << the_ring.size() << std::endl;

	//test_printlattice(accelerator);
	//test_findm66(accelerator);
	//test_linepass(accelerator);
	//test_ringpass(accelerator);
	//test_linepass_tpsa(the_ring);
	//test_findorbit6(accelerator);
	//test_dynap_xy(the_ring);
	//test_read_flat_file(accelerator);

	//test_cmd_dynap_xy();
	//test_cmd_dynap_ex();
	//test_cmd_dynap_ma();
	//test_cmd_track_linepass();
	//test_kicktable(accelerator);
	test_simple_drift();
	test_simple_quadrupole();

	return 0;

}
