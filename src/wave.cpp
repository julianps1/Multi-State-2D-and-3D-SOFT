#include "wave.hpp"
#include <fstream>

namespace wave {
    double alfa,beta,q0,y0,px0,py0,dt;

    void read(std::string const& fname) {
        std::ifstream ifs(fname);
        // parse file and assign to alfa, … etc.
    }
}
