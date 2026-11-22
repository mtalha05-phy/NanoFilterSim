// Usage: nfs_run case.txt out.csv
#include <iostream>
#include "nanofilter/case_io.hpp"
int main(int argc, char** argv) {
    if (argc != 3) { std::cerr << "usage: nfs_run case.txt out.csv\n"; return 2; }
    try {
        nfs::writeCsv(argv[2], nfs::runCase(nfs::readCase(argv[1])));
    } catch (const std::exception& e) { std::cerr << "error: " << e.what() << "\n"; return 1; }
    return 0;
}
