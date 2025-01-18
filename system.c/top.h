#ifndef TOP_H
#define TOP_H

#include "systemc.h"
#include "stim.h"
#include "mult.h"
#include "mon.h"

SC_MODULE(Top) {
    // Signal declarations
    sc_signal<int> asig, bsig, fsig;
    sc_clock testclk;

    // Module instantiations
    Stim stim1;
    Mult uut;
    Mon mon1;

    // Constructor
    SC_CTOR(Top) : stim1("stim1"), uut("uut"), mon1("mon1"), testclk("testclk", 10, SC_NS) {
        // Connect the modules
        stim1.a(asig);
        stim1.b(bsig);
        stim1.Clk(testclk);

        uut.a(asig);
        uut.b(bsig);
        uut.f(fsig);

        mon1.a(asig);
        mon1.b(bsig);
        mon1.f(fsig);
        mon1.Clk(testclk);
    }
};

#endif // TOP_H

