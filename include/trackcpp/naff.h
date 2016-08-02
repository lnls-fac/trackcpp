#ifndef _NAFF_H
#define _NAFF_H

#include <trackcpp/pos.h>
#include <vector>

/* Frequency Map Analysis */
void naff_traj(const std::vector<Pos<double>>& data, double& tunex, double& tuney);
void naff_general(const std::vector<double>& re,
                  const std::vector<double>& im,
                  bool is_real,
                  int nr_ff,
                  int win,
                  std::vector<double>& ff_out,
                  std::vector<double>& re_out,
                  std::vector<double>& im_out);

#endif
