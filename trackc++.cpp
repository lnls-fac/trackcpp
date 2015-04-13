#include "trackc++.h"
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
