#pragma once
#include <string>

namespace grid {
    extern double xmin, ymin, xmax, ymax, dx, dy, an0;
    void read(std::string const& fname);
    // Other grid funcs
}
