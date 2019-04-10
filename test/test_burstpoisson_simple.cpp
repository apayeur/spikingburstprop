//
//  test_burstpoisson_simple.cpp
//  
//
//  Created by Alexandre Payeur on 4/9/19.
//

/*
 * Simpler test of the simplified two-compartment bursting neuron
 * defined
 */

#include "auryn.h"
#include "BurstPoissonGroup.h"

using namespace auryn;

namespace po = boost::program_options;

int main(int ac, char* av[])
{
    int errcode = 0;
    char strbuf [255];
    unsigned int seed = 1;
    string dir = "./test-output/";
    string simname = "simpleburstpoisson";
    
    //***          Options            ***//
    try {
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help", "produce help message")
        ("dir", po::value<string>(), "output directory")
        ("seed", po::value<int>(), "random seed")
        ;
        
        po::variables_map vm;
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);
        
        if (vm.count("help")) {
            std::cout << desc << "\n";
            return 1;
        }
        
        if (vm.count("dir")) {
            std::cout << "dir set to "
            << vm["dir"].as<string>() << ".\n";
            dir = vm["dir"].as<string>();
        }
        
        if (vm.count("seed")) {
            std::cout << "seed set to "
            << vm["seed"].as<int>() << ".\n";
            seed = vm["seed"].as<int>();
        }
    }
    catch(std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        std::cerr << "Exception of unknown type!\n";
    }
    auryn_init( ac, av, dir, simname );
    
    // define 2-comp neuron group
    NeuronID number_of_neurons = 1;
    BurstPoissonGroup* pyr = new BurstPoissonGroup(number_of_neurons);
 
    // define current input
    CurrentInjector * curr_inject_soma = new CurrentInjector(pyr, "mem");
    CurrentInjector * curr_inject_dend = new CurrentInjector(pyr, "Vd");
    
 
    
    // define monitors
    VoltageMonitor * vmon   = new VoltageMonitor( pyr, 0,   sys->fn("mem") );
    StateMonitor * smon_vd  = new StateMonitor( pyr, 0, "Vd", sys->fn("Vd") );
    
    // run simulation
    logger->msg("Running ...",PROGRESS);
    
    // The alternating currents switched between a maximum of 175 pA and a minimum of -210 pA in the dendrites and between 230 pA and 130 pA in the soma. Dendritic input during 75e-3 s and somatic input during 50.e-3 seconds.
    const double max_dendritic_current = 1000e-12; //90e-12;
    const double min_dendritic_current = 800e-12;
    const double max_somatic_current = 400.e-12; //650e-12;
    const double min_somatic_current = 300e-12;
    
    const double period = 200e-3;
    const double simtime = 1000e-3;
    const double segtime_maxsoma = period/2.;
    const double segtime_minsoma = period/2;
    const double segtime_maxdend = 0.7*period;
    const double segtime_mindend = 0.3*period;
    const double small_overlap = 0.1*period;
    
    // simulate
    // pre stim
    
    curr_inject_soma->set_all_currents(min_somatic_current/pyr[0].get_Cs());
    curr_inject_dend->set_all_currents(min_dendritic_current/pyr[0].get_Cd());
    sys->run(simtime);
    
    
    // simulate current
    
    for (int i=0;i<10;i++)
    {
        curr_inject_soma->set_all_currents(max_somatic_current/pyr[0].get_Cs());
        curr_inject_dend->set_all_currents(min_dendritic_current/pyr[0].get_Cd());
        sys->run(segtime_mindend - small_overlap);
        
        curr_inject_soma->set_all_currents(max_somatic_current/pyr[0].get_Cs());
        curr_inject_dend->set_all_currents(max_dendritic_current/pyr[0].get_Cd());
        sys->run(segtime_maxsoma  - (segtime_mindend - small_overlap) );
        
        curr_inject_soma->set_all_currents(min_somatic_current/pyr[0].get_Cs());
        curr_inject_dend->set_all_currents(max_dendritic_current/pyr[0].get_Cd());
        sys->run(segtime_maxdend + segtime_mindend - small_overlap - segtime_maxsoma);
        
        curr_inject_soma->set_all_currents(min_somatic_current/pyr[0].get_Cs());
        curr_inject_dend->set_all_currents(min_dendritic_current/pyr[0].get_Cd());
        sys->run(small_overlap);
    }
    
    
    
    if (errcode)
        auryn_abort(errcode);
    
    logger->msg("Freeing ...",PROGRESS,true);
    auryn_free();
    return errcode;
}
