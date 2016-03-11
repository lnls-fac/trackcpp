#ifndef _NAFF_H
#define _NAFF_H

#include <trackcpp/pos.h>
#include <vector>

/* Frequency Map Analysis */
void naff_run(const std::vector<Pos<double>>& data, double& tunex, double& tuney);

#endif
