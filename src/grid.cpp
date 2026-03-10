#include "grid.hpp"
#include <fstream>

namespace grid {
    double xmin, ymin, xmax, ymax, dx, dy, an0;

    void read(std::string const& fname) {
        std::ifstream ifs(fname);
        // parse file and assign to xmin, … etc.
    }
}
