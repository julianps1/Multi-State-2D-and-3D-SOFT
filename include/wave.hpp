#pragma once
#include <string>

namespace wave {
    extern double alfa, beta, q0, y0, px0, py0, dt;
    void read(std::string const& fname);
    // … other wf‑related functions …
}
