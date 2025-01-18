#ifndef MON_H
#define MON_H

#include "systemc.h"

SC_MODULE(Mon) {
    sc_in<int> a;   // Input port for the first operand
    sc_in<int> b;   // Input port for the second operand
    sc_in<int> f;   // Input port for the result from the multiplier
    sc_in<bool> Clk;  // Clock signal

    void monitor();

    SC_CTOR(Mon) {
        SC_METHOD(monitor);
        sensitive << Clk.pos();  // Monitor on the negative clock edge
    }
};

#endif // MON_H

