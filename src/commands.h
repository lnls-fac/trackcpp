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

#ifndef _COMMANDS_H
#define _COMMANDS_H

#include <vector>
#include <string>

int cmd_tests          (const std::vector<std::string>& args);
int cmd_dynap_xy       (const std::vector<std::string>& args);
int cmd_dynap_ex       (const std::vector<std::string>& args);
int cmd_dynap_ma       (const std::vector<std::string>& args);
int cmd_dynap_pxa      (const std::vector<std::string>& args);
int cmd_dynap_pya      (const std::vector<std::string>& args);
int cmd_dynap_xyfmap   (const std::vector<std::string>& args);
int cmd_dynap_exfmap   (const std::vector<std::string>& args);
int cmd_track_linepass (const std::vector<std::string>& args);

#endif
