#ifndef MULT_H
#define MULT_H

#include "systemc.h"

SC_MODULE(Mult) {
    sc_in<int> a;
    sc_in<int> b;
    sc_out<int> f;

    void action() {
        f.write(a.read() * b.read());
    }

    SC_CTOR(Mult) {
        SC_METHOD(action);
        sensitive << a << b;
    }
};

#endif // MULT_H

