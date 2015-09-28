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
    if (cmd == "tests")    return cmd_tests(args);
    if (cmd == "dynap_xy") return cmd_dynap_xy(args);
    if (cmd == "dynap_ex") return cmd_dynap_ex(args);
    if (cmd == "dynap_ma") return cmd_dynap_ma(args);
    if (cmd == "dynap_ma2") return cmd_dynap_ma2(args);
    if (cmd == "dynap_pxa") return cmd_dynap_pxa(args);
    if (cmd == "dynap_pya") return cmd_dynap_pya(args);
    if (cmd == "dynap_xyfmap") return cmd_dynap_xyfmap(args);
    if (cmd == "dynap_exfmap") return cmd_dynap_exfmap(args);
    if (cmd == "track_linepass") return cmd_track_linepass(args);
    std::cerr << "trackcpp: invalid command!" << std::endl;
    return EXIT_FAILURE;
}
