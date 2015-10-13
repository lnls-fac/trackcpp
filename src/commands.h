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

#ifndef _COMMANDS_H
#define _COMMANDS_H

#include <vector>
#include <string>

int cmd_tests          (const std::vector<std::string>& args);
int cmd_dynap_xy       (const std::vector<std::string>& args);
int cmd_dynap_ex       (const std::vector<std::string>& args);
int cmd_dynap_acceptance (const std::vector<std::string>& args);
//int cmd_dynap_ma       (const std::vector<std::string>& args);
//int cmd_dynap_pxa      (const std::vector<std::string>& args);
//int cmd_dynap_pya      (const std::vector<std::string>& args);
int cmd_dynap_xyfmap   (const std::vector<std::string>& args);
int cmd_dynap_exfmap   (const std::vector<std::string>& args);
int cmd_track_linepass (const std::vector<std::string>& args);
int cmd_help           (const std::vector<std::string>& args);

#endif
