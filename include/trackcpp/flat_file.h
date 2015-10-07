// The MIT License (MIT)
//
// Copyright (c) 2015 LNLS Accelerator Division
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#ifndef _FLAT_FILE_H
#define _FLAT_FILE_H

#include "accelerator.h"
#include "elements.h"
#include "auxiliary.h"
#include <string>
#include <vector>
#include <iomanip>

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

//Status::type read_flat_file(const std::string& filename, Accelerator& accelerator);

Status::type read_flat_file(const std::string& filename, Accelerator& accelerator, bool file_flag = true);
Status::type read_flat_file(std::string& filename, Accelerator& accelerator, bool file_flag = true);

Status::type write_flat_file(const std::string& filename, const Accelerator& accelerator, bool file_flag = true);
Status::type write_flat_file(std::string& filename, const Accelerator& accelerator, bool file_flag = true);

#endif
