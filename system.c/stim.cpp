#include "stim.h"

void Stim::stimulus() {
    wait();  // Wait for a stable state at the beginning of simulation

    // Series of test vectors
    a.write(1); b.write(42); wait(1, SC_NS);
    a.write(2); b.write(21); wait(1, SC_NS);
    a.write(3); b.write(14); wait(1, SC_NS);
    a.write(6); b.write(7);  wait(1, SC_NS);
    a.write(7); b.write(6);  wait(1, SC_NS);
    a.write(14); b.write(3); wait(1, SC_NS);
    a.write(21); b.write(2); wait(1, SC_NS);
    a.write(42); b.write(1); wait(1, SC_NS);

    sc_stop();  // Stop the simulation after sending all vectors
}

