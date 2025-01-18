#include <systemc.h>
#include "top.h"

int sc_main(int argc, char* argv[]) {
    Top top("Top");  // Instantiate your top module

    // Setup trace file
    sc_trace_file *tf = sc_create_vcd_trace_file("simulation_trace");  // Create a VCD trace file
    sc_trace(tf, top.asig, "asig");
    sc_trace(tf, top.bsig, "bsig");
    sc_trace(tf, top.fsig, "fsig");
    // Add more signals as needed

    // Start the simulation
    sc_start(100, SC_NS);  // Run for a certain time or until simulation stops

    // Close the trace file
    sc_close_vcd_trace_file(tf);

    return 0;
}

