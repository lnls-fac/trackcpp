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

#include "commands.h"
#include <trackcpp/trackcpp.h>

int main(int argc, char *argv[]) {
    if (argc == 1) {
        print_header (stdout);
        return EXIT_SUCCESS;
    };

    std::vector<std::string> args;
    for(int i=0; i<argc; ++i) args.push_back(std::string(argv[i]));

    std::string cmd(args[1]);
    if (cmd == "help") return cmd_help(args);
    if (cmd == "tests") return cmd_tests(args);
    if (cmd == "dynap_xy") return cmd_dynap_xy(args);
    if (cmd == "dynap_ex") return cmd_dynap_ex(args);
    // if (cmd == "dynap_ma") return cmd_dynap_ma(args);
    // if (cmd == "dynap_pxa") return cmd_dynap_pxa(args);
    // if (cmd == "dynap_pya") return cmd_dynap_pya(args);
    if (cmd == "dynap_ma") return cmd_dynap_acceptance(args);
    if (cmd == "dynap_pxa") return cmd_dynap_acceptance(args);
    if (cmd == "dynap_pya") return cmd_dynap_acceptance(args);
    if (cmd == "dynap_xyfmap") return cmd_dynap_xyfmap(args);
    if (cmd == "dynap_exfmap") return cmd_dynap_exfmap(args);
    if (cmd == "track_linepass") return cmd_track_linepass(args);
    std::cerr << "trackcpp: invalid command!" << std::endl;
    return EXIT_FAILURE;
}
