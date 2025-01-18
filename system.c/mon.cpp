#include "mon.h"
#include <iostream>

void Mon::monitor() {
    std::cout << "At time " << sc_time_stamp() << ": Multiplying " << a.read()
              << " * " << b.read() << " results in " << f.read() << std::endl;
    assert((a.read() * b.read() == f.read()) && "Error: Incorrect multiplication result.");
}

