#include "trackcpp.h"


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
    //if (cmd == "dynap_ma_threads") return cmd_dynap_ma_threads(args);
    //if (cmd == "dynap_xy") return cmd_dynap_xy(args);
    //if (cmd == "dynap_ex") return cmd_dynap_ex(args);
    if (cmd == "dynap_ma") return cmd_dynap_ma(args);
    if (cmd == "dynap_xyfmap") return cmd_dynap_xyfmap(args);
    if (cmd == "dynap_exfmap") return cmd_dynap_exfmap(args);
    if (cmd == "track_linepass") return cmd_track_linepass(args);
    std::cerr << "trackcpp: invalid command!" << std::endl;
    return EXIT_FAILURE;
}
