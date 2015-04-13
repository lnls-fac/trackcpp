#ifndef _FLAT_FILE_H
#define _FLAT_FILE_H

#include "accelerator.h"
#include "elements.h"
#include "auxiliary.h"
#include <string>
#include <vector>

struct FlatFileType {
		enum type_ {
			marker    = -1,
			drift     =  0,
			mpole     =  1,
			cavity    =  2,
			corrector =  3,
			thin_kick =  3,
			kicktable =  6
		};
	};


Status::type read_flat_file_tracy(const std::string& filename, Accelerator& accelerator);
Status::type read_flat_file_trackcpp(const std::string& filename, Accelerator& accelerator);
Status::type read_flat_file(const std::string& filename, Accelerator& accelerator);

#endif
