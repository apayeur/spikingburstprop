/* sim_direct_transmit.cpp
*
* Figure S2: Comparison of direct event or burst transmission and
* transmission via STP.
*
* To plot with gnuplot, do: p 'output0.0.mem' u 1:($2*1000) w l, 'transbr.0.ras_input' u 1:($2-70) w dots
*
* Created by Alexandre Payeur.
*/

#include "auryn.h"
#include "TransmitBurstConnection.h"
#include "TransmitEventConnection.h"
#include "STPeTMConnection.h"

#define N 1

using namespace auryn;

namespace po = boost::program_options;

int main(int ac, char* av[]) 
{

	int errcode = 0;
	string simname = "transbr";
	string logfile = simname;
	string tmpstr;

	// BEGIN Global definitions
	auryn_init( ac, av );
	sys->set_simulation_name(simname);
	// END Global definitions

    // define bursting input group
    NaudGroup * input_neurons = new NaudGroup(10);
    input_neurons->random_mem(-70e-3, 5.e-3);
    
    // define Poisson noise groups
    PoissonGroup * poisson_soma = new PoissonGroup(4000, 1);
    PoissonGroup * poisson_dend = new PoissonGroup(100, 1.);
    
	// define receiving group for FileCurrent
	IFGroup * output_neuron = new IFGroup(1);

    // define burst-mediated connection
    // Note: comment or uncomments the different blocks below according to what should be tested. Lazy way to do that...
    
    //TransmitBurstConnection * conn = new TransmitBurstConnection(input_neurons, output_neuron, 0.1, 0.1, GLUT);
    //TransmitEventConnection * conn = new TransmitEventConnection(input_neurons, output_neuron, 0.1, 0.1, GLUT);
    
    STPeTMConnection * conn = new STPeTMConnection(input_neurons, output_neuron, 0.1, 0.1, GLUT);
    //STD
    //conn->set_target("g_ampa");
    //conn->set_tau_d(20.e-3);
    //conn->set_tau_f(1.);
    //conn->set_ujump(0.9);
    //conn->set_f(0.1);
    //STF
    conn->set_target("g_ampa");
    conn->set_tau_d(100e-3);
    conn->set_tau_f(100e-3);
    conn->set_ujump(0.02);
    conn->set_f(0.1);
    
    // define Poisson-to-NaudGroup connections
    SparseConnection * con_ext_dend = new SparseConnection(poisson_dend, input_neurons, 0.1, 0.1);
    con_ext_dend->set_target("g_ampa_dend");
    
    SparseConnection * con_ext_soma = new SparseConnection(poisson_soma, input_neurons, 0.15, 0.1);
    con_ext_soma->set_target("g_ampa");
    
	// define monitors
	VoltageMonitor * vmon0 = new VoltageMonitor( output_neuron, 0, sys->fn("output",0,"memSTF"), 1e-3 );
    VoltageMonitor * vmon_input = new VoltageMonitor( input_neurons, 0, sys->fn("input",0,"mem"), 1e-3 );
    SpikeMonitor * spkmon = new SpikeMonitor ( input_neurons, sys->fn("ras_input") );
    
    // current injectors
    CurrentInjector * curr_soma = new CurrentInjector(input_neurons, "mem");
    CurrentInjector * curr_dend = new CurrentInjector(input_neurons, "Vd");
    
    curr_soma->set_all_currents(250e-12/input_neurons[0].get_Cs());
    curr_dend->set_all_currents(50e-12/input_neurons[0].get_Cd());

    // get presynaptic partner of postsynaptic neuron
    NeuronID pre = 0;
    std::vector< neuron_pair > presyn_partners = conn->get_pre_partners(pre);
    for (int i=0;i<presyn_partners.size();i++){
        std::cout << presyn_partners[i].i << std::endl;
    }
    
	// run simulation
	logger->msg("Running ...",PROGRESS);
	sys->run(5.0);

	if (errcode)
		auryn_abort(errcode);

	logger->msg("Freeing ...",PROGRESS,true);
	auryn_free();
	return errcode;
}
