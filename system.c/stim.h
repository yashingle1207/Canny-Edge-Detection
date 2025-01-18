#ifndef STIM_H
#define STIM_H

#include "systemc.h"

SC_MODULE(Stim) {
    sc_out<int> a;  // Output port for the first integer
    sc_out<int> b;  // Output port for the second integer
    sc_in<bool> Clk;  // Clock input

    void stimulus();

    SC_CTOR(Stim) {
        SC_THREAD(stimulus);
        sensitive << Clk.pos();  // Trigger on the positive edge of the clock
    }
};

#endif // STIM_H

