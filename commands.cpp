#include "output.h"
#include "commands.h"
#include "dynap.h"
#include "flat_file.h"
#include "tracking.h"
#include "lattice.h"
#include "accelerator.h"
#include "elements.h"
#include "auxiliary.h"
#include <string>
#include <cstdlib>
#include <iostream>

int cmd_dynap_xy(const std::vector<std::string>& args) {

	if (args.size() != 16) {
		std::cerr << "dynap_xy: invalid number of arguments!" << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << std::endl;
	std::cout << "[cmd_dynap_xy]" << std::endl << std::endl;

	std::string  flat_filename(args[2]);
	double       ring_energy     = std::atof(args[3].c_str());
	unsigned int harmonic_number = std::atoi(args[4].c_str());
	std::string  cavity_state(args[5]);
	std::string  radiation_state(args[6]);
	std::string  vchamber_state(args[7]);
	double       de = std::atof(args[8].c_str());
	unsigned int nr_turns = std::atoi(args[9].c_str());
	unsigned int x_nrpts = std::atoi(args[10].c_str());
	double       x_min = std::atof(args[11].c_str());
	double       x_max = std::atof(args[12].c_str());
	unsigned int y_nrpts = std::atoi(args[13].c_str());
	double       y_min = std::atof(args[14].c_str());
	double       y_max = std::atof(args[15].c_str());

	print_header(stdout);
	std::cout << std::endl;
	std::cout << "flat_filename   : " << flat_filename << std::endl;
	std::cout << "energy[eV]      : " << ring_energy << std::endl;
	std::cout << "harmonic_number : " << harmonic_number << std::endl;
	std::cout << "cavity_state    : " << cavity_state << std::endl;
	std::cout << "radiation_state : " << radiation_state << std::endl;
	std::cout << "vchamber_state  : " << vchamber_state << std::endl;
	std::cout << "de              : " << de << std::endl;
	std::cout << "nr_turns        : " << nr_turns << std::endl;
	std::cout << "x_nrpts         : " << x_nrpts << std::endl;
	std::cout << "x_min[m]        : " << x_min << std::endl;
	std::cout << "x_max[m]        : " << x_max << std::endl;
	std::cout << "y_nrpts         : " << y_nrpts << std::endl;
	std::cout << "y_min[m]        : " << y_min << std::endl;
	std::cout << "y_max[m]        : " << y_max << std::endl;

	std::cout << std::endl;
	std::cout << get_timestamp() << " begin timestamp" << std::endl;

	Accelerator accelerator;

	// reads flat file
	Status::type status = read_flat_file(flat_filename, accelerator);
	if (status == Status::file_not_found) {
		std::cerr << "dynap_xy: flat file not found!" << std::endl;
		return EXIT_FAILURE;
	}
	std::cout << get_timestamp() << " input file with flat lattice read." << std::endl;

	// builds accelerator
	accelerator.energy = ring_energy;
	accelerator.harmonic_number = harmonic_number;
	accelerator.cavity_on = (cavity_state == "on");
	accelerator.radiation_on = (radiation_state == "on");
	accelerator.vchamber_on = (vchamber_state == "on");

	// calcs dynamical aperture
	std::vector<Pos<double> > cod;
	Pos<double> p0(0,0,0,0,de,0);
	std::vector<DynApGridPoint> grid;
	dynap_xy(accelerator, cod, nr_turns, p0, x_nrpts, x_min, x_max, y_nrpts, y_min, y_max, true, grid);

	// generates output files
	std::cout << get_timestamp() << " saving closed-orbit to file" << std::endl;
	status = print_closed_orbit(accelerator, cod);
	if (status == Status::file_not_opened) return status;
	std::cout << get_timestamp() << " saving dynap_xy grid to file" << std::endl;
	status = print_dynapgrid (accelerator, grid, "[dynap_xy]", "dynap_xy_out.txt");
	if (status == Status::file_not_opened) return status;

	std::cout << get_timestamp() << " end timestamp" << std::endl;
	return EXIT_SUCCESS;

}

int cmd_dynap_ex(const std::vector<std::string>& args) {

	if (args.size() != 16) {
		std::cerr << "dynap_ex: invalid number of arguments!" << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << std::endl;
	std::cout << "[cmd_dynap_ex]" << std::endl << std::endl;

	std::string  flat_filename(args[2]);
	double       ring_energy     = std::atof(args[3].c_str());
	unsigned int harmonic_number = std::atoi(args[4].c_str());
	std::string  cavity_state(args[5]);
	std::string  radiation_state(args[6]);
	std::string  vchamber_state(args[7]);
	double       y = std::atof(args[8].c_str());
	unsigned int nr_turns = std::atoi(args[9].c_str());
	unsigned int e_nrpts = std::atoi(args[10].c_str());
	double       e_min = std::atof(args[11].c_str());
	double       e_max = std::atof(args[12].c_str());
	unsigned int x_nrpts = std::atoi(args[13].c_str());
	double       x_min = std::atof(args[14].c_str());
	double       x_max = std::atof(args[15].c_str());

	print_header(stdout);
	std::cout << std::endl;
	std::cout << "flat_filename   : " << flat_filename << std::endl;
	std::cout << "energy[eV]      : " << ring_energy << std::endl;
	std::cout << "harmonic_number : " << harmonic_number << std::endl;
	std::cout << "cavity_state    : " << cavity_state << std::endl;
	std::cout << "radiation_state : " << radiation_state << std::endl;
	std::cout << "vchamber_state  : " << vchamber_state << std::endl;
	std::cout << "y[m]            : " << y << std::endl;
	std::cout << "nr_turns        : " << nr_turns << std::endl;
	std::cout << "e_nrpts         : " << e_nrpts << std::endl;
	std::cout << "e_min           : " << e_min << std::endl;
	std::cout << "e_max           : " << e_max << std::endl;
	std::cout << "x_nrpts         : " << x_nrpts << std::endl;
	std::cout << "x_min[m]        : " << x_min << std::endl;
	std::cout << "x_max[m]        : " << x_max << std::endl;

	std::cout << std::endl;
	std::cout << get_timestamp() << " begin timestamp" << std::endl;

	Accelerator accelerator;

	// reads flat file
	Status::type status = read_flat_file(flat_filename, accelerator);
	if (status == Status::file_not_found) {
		std::cerr << "dynap_ex: flat file not found!" << std::endl;
		return EXIT_FAILURE;
	}

	// builds accelerator
	accelerator.energy = ring_energy;
	accelerator.harmonic_number = harmonic_number;
	accelerator.cavity_on = (cavity_state == "on");
	accelerator.radiation_on = (radiation_state == "on");
	accelerator.vchamber_on = (vchamber_state == "on");

	// calcs dynamical aperture
	std::vector<Pos<double> > cod;
	Pos<double> p0(0,0,y,0,0,0);
	std::vector<DynApGridPoint> grid;
	dynap_ex(accelerator, cod, nr_turns, p0, e_nrpts, e_min, e_max, x_nrpts, x_min, x_max, true, grid);

	// generates output files
	std::cout << get_timestamp() << " saving closed-orbit to file" << std::endl;
	status = print_closed_orbit(accelerator, cod);
	if (status == Status::file_not_opened) return status;
	std::cout << get_timestamp() << " saving dynap_ex grid to file" << std::endl;
	status = print_dynapgrid (accelerator, grid, "[dynap_ex]", "dynap_ex_out.txt");
	if (status == Status::file_not_opened) return status;

	std::cout << get_timestamp() << " end timestamp" << std::endl;
	return EXIT_SUCCESS;

}

int cmd_dynap_ma(const std::vector<std::string>& args) {

	if (args.size() < 15) {
		std::cerr << "dynap_ex: invalid number of arguments!" << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << std::endl;
	std::cout << "[cmd_dynap_ma]" << std::endl << std::endl;

	std::string  flat_filename(args[2]);
	double       ring_energy      = std::atof(args[3].c_str());
	unsigned int harmonic_number = std::atoi(args[4].c_str());
	std::string  cavity_state(args[5].c_str());
	std::string  radiation_state(args[6].c_str());
	std::string  vchamber_state(args[7].c_str());
	unsigned int nr_turns = std::atoi(args[8].c_str());
	double       y0    = std::atof(args[9].c_str());
	double       e0    = std::atof(args[10].c_str());
	double       e_tol = std::atof(args[11].c_str());
	double       s_min = std::atof(args[12].c_str());
	double       s_max = std::atof(args[13].c_str());
	std::vector<std::string> fam_names;
	for(unsigned int i=14; i<args.size(); ++i) fam_names.push_back(args[i]);


	print_header(stdout);
	std::cout << std::endl;
	std::cout << "flat_filename   : " << flat_filename << std::endl;
	std::cout << "energy[eV]      : " << ring_energy << std::endl;
	std::cout << "harmonic_number : " << harmonic_number << std::endl;
	std::cout << "cavity_state    : " << cavity_state << std::endl;
	std::cout << "radiation_state : " << radiation_state << std::endl;
	std::cout << "vchamber_state  : " << vchamber_state << std::endl;
	std::cout << "nr_turns        : " << nr_turns << std::endl;
	std::cout << "y0[m]           : " << y0 << std::endl;
	std::cout << "e0              : " << e0 << std::endl;
	std::cout << "e_tol           : " << e_tol << std::endl;
	std::cout << "s_min[m]        : " << s_min << std::endl;
	std::cout << "s_max[m]        : " << s_max << std::endl;
	std::cout << "fam_names       : ";
	for(unsigned int i=0; i<fam_names.size(); ++i) std::cout << fam_names[i] << " "; std::cout << std::endl;

	std::cout << std::endl;
	std::cout << get_timestamp() << " begin timestamp" << std::endl;

	Accelerator accelerator;

	// reads flat file
	Status::type status = read_flat_file(flat_filename, accelerator);
	if (status == Status::file_not_found) {
		std::cerr << "dynap_ma: flat file not found!" << std::endl;
		return EXIT_FAILURE;
	}

	// builds accelerator
	accelerator.energy = ring_energy;
	accelerator.harmonic_number = harmonic_number;
	accelerator.cavity_on = (cavity_state == "on");
	accelerator.radiation_on = (radiation_state == "on");
	accelerator.vchamber_on = (vchamber_state == "on");

	// calcs dynamical aperture
	std::vector<Pos<double> > cod;
	Pos<double> p0(0,0,y0,0,0,0);
	std::vector<DynApGridPoint> grid;
	dynap_ma(accelerator, cod, nr_turns, p0, e0, e_tol, s_min, s_max, fam_names, true, grid);

	// generates output files
	std::cout << get_timestamp() << " saving closed-orbit to file" << std::endl;
	status = print_closed_orbit(accelerator, cod);
	if (status == Status::file_not_opened) return status;
	std::cout << get_timestamp() << " saving dynap_ex grid to file" << std::endl;
	status = print_dynapgrid (accelerator, grid, "[dynap_ma]", "dynap_ma_out.txt");
	if (status == Status::file_not_opened) return status;

	std::cout << get_timestamp() << " end timestamp" << std::endl;
	return EXIT_SUCCESS;

}

int cmd_track_linepass(const std::vector<std::string>& args) {

	if (args.size() != 15) {
		std::cerr << "cmd_track_linepass: invalid number of arguments!" << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << std::endl;
	std::cout << "[cmd_track_linepass]" << std::endl << std::endl;

	std::string  flat_filename(args[2]);
	double       ring_energy     = std::atof(args[3].c_str());
	unsigned int harmonic_number = std::atoi(args[4].c_str());
	std::string  cavity_state(args[5]);
	std::string  radiation_state(args[6]);
	std::string  vchamber_state(args[7]);
	unsigned int start_element = std::atoi(args[8].c_str());
	double       rx0 = std::atof(args[9].c_str());
	double       px0 = std::atof(args[10].c_str());
	double       ry0 = std::atof(args[11].c_str());
	double       py0 = std::atof(args[12].c_str());
	double       de0 = std::atof(args[13].c_str());
	double       dl0 = std::atof(args[14].c_str());

	print_header(stdout);
	std::cout << std::endl;
	std::cout << "flat_filename   : " << flat_filename << std::endl;
	std::cout << "energy[eV]      : " << ring_energy << std::endl;
	std::cout << "harmonic_number : " << harmonic_number << std::endl;
	std::cout << "cavity_state    : " << cavity_state << std::endl;
	std::cout << "radiation_state : " << radiation_state << std::endl;
	std::cout << "vchamber_state  : " << vchamber_state << std::endl;
	std::cout << "start_element   : " << start_element << std::endl;
	std::cout << "rx0[m]          : " << rx0 << std::endl;
	std::cout << "px0[rad]        : " << px0 << std::endl;
	std::cout << "ry0[m]          : " << ry0 << std::endl;
	std::cout << "py0[rad]        : " << py0 << std::endl;
	std::cout << "de0             : " << de0 << std::endl;
	std::cout << "dl0[m]          : " << dl0 << std::endl;

	std::cout << std::endl;
	std::cout << get_timestamp() << " begin timestamp" << std::endl;

	Accelerator accelerator;

	// reads flat file
	Status::type status = read_flat_file(flat_filename, accelerator);
	if (status == Status::file_not_found) {
		std::cerr << "track_linepass: flat file not found!" << std::endl;
		return EXIT_FAILURE;
	}
	if (status == Status::flat_file_error) {
		std::cerr << "track_linepass: flat file with incorrect syntax!" << std::endl;
		return EXIT_FAILURE;
	}
	std::cout << get_timestamp() << " input file with flat lattice read." << std::endl;

	// builds accelerator
	accelerator.energy = ring_energy;
	accelerator.harmonic_number = harmonic_number;
	accelerator.cavity_on = (cavity_state == "on");
	accelerator.radiation_on = (radiation_state == "on");
	accelerator.vchamber_on = (vchamber_state == "on");

	// does tracking
	Pos<double> pos(rx0,px0,ry0,py0,de0,dl0);
	std::vector<Pos<double>> pos_list;
	Plane::type lost_plane;
	unsigned int offset_element = start_element;
	track_linepass(accelerator, pos, pos_list, offset_element, lost_plane, true);

	std::cout << get_timestamp() << " saving track_linepass data to file" << std::endl;
	status = print_tracking_linepass(accelerator, pos_list, start_element, "track_linepass_out.txt");
	if (status == Status::file_not_opened) return status;

	std::cout << get_timestamp() << " end timestamp" << std::endl;
	return EXIT_SUCCESS;

}
