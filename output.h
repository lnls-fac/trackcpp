#ifndef _OUTPUT_H
#define _OUTPUT_H

#include "dynap.h"
#include "accelerator.h"
#include "elements.h"
#include <vector>
#include <string>
#include <cstdlib>

Status::type print_closed_orbit      (const Accelerator& accelerator, const std::vector<Pos<double>>&    cod,  const std::string& filename = "cod_out.txt");
Status::type print_dynapgrid         (const Accelerator& accelerator, const std::vector<DynApGridPoint>& grid, const std::string& label, const std::string& filename = "dynap_out.txt");
Status::type print_tracking_ringpass (const Accelerator& accelerator, const std::vector<Pos<double>>& points, const std::string& filename = "track_linepass_out.txt");
Status::type print_tracking_linepass (const Accelerator& accelerator, const std::vector<Pos<double>>& points, const unsigned int start_element, const std::string& filename);
void         print_header            (FILE* fp);

#endif
