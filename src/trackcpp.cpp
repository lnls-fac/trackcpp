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

#include <trackcpp/trackcpp.h>
#include <ctime>
#include <cstdlib>
#include <string>


bool verbose_on = true;

bool isfinite(const double& v) {
	return std::isfinite(v);
}

std::string get_timestamp() {

	char buffer[30];
	time_t t = time(NULL);   // get time now
	struct tm* now = localtime(&t);
	sprintf(buffer, "[%04i-%02i-%02i %02i:%02i:%02i]", 1900+now->tm_year, 1+now->tm_mon, now->tm_mday, now->tm_hour, now->tm_min, now->tm_sec);
	return std::string(buffer);
}
