/*
 * Test of the simplified two-compartment bursting neuron
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
    string simname = "burstpoisson";
    
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
    NeuronID number_of_neurons = 4000;
    BurstPoissonGroup* pyr = new BurstPoissonGroup(number_of_neurons);
    pyr->syn_exc_soma->set_nmda_ampa_current_ampl_ratio(0.);
    pyr->syn_exc_dend->set_nmda_ampa_current_ampl_ratio(0.);
    pyr->random_mem(-70e-3, 5.e-3);
    
    // define external population of Poisson neurons
    NeuronID nb_exc_Poisson_neurons = 25000;
    NeuronID nb_inh_Poisson_neurons = nb_exc_Poisson_neurons;
    float poisson_rate_exc = 5.;
    float poisson_rate_inh = 5.;
    PoissonGroup* exc_Poisson = new PoissonGroup(nb_exc_Poisson_neurons, poisson_rate_exc);
    PoissonGroup* inh_Poisson = new PoissonGroup(nb_inh_Poisson_neurons, poisson_rate_inh);
    exc_Poisson->seed(seed);
    inh_Poisson->seed(seed);
    
    // define current input
    CurrentInjector * curr_inject_soma = new CurrentInjector(pyr, "mem");
    CurrentInjector * curr_inject_dend = new CurrentInjector(pyr, "Vd");
    
    // specify connectivity
    // These connections yield an avg ER = 4-5 Hz without injected current, with subthreshold
    // voltage fluctuations around 4 mV (Destexhe et al.). The BP at zero injected dendritic and somatic
    // currents is around 12%, and Vd fluctuations are about 2.5 mV when disregarding bAP and Ca spikes (thresholding the distribution std(Vd[Vd<-60mV])).
    int number_of_ext_conn[2] = {100, 30}; //exc, inh
    
    float p_exc = float(number_of_ext_conn[0])/float(nb_exc_Poisson_neurons);
    float p_inh = float(number_of_ext_conn[1])/float(nb_inh_Poisson_neurons);
    
    const float w_exc = 0.17;        // conductance amplitude in units of leak conductance
    const float w_inh = 0.26;

    SparseConnection * con_ext_exc_soma = new SparseConnection(exc_Poisson, pyr, w_exc, p_exc, GLUT);
    con_ext_exc_soma->set_target("g_ampa");
    SparseConnection * con_ext_inh_soma = new SparseConnection(inh_Poisson, pyr, w_inh, p_inh, GABA);
    con_ext_inh_soma->set_target("g_gaba");
    
    const float w_exc_dend = 0.0425;        // conductance amplitude in units of leak conductance
    const float w_inh_dend = 0.065;
    
    SparseConnection * con_ext_exc_dend = new SparseConnection(exc_Poisson, pyr, w_exc_dend, p_exc, GLUT);
    con_ext_exc_dend->set_target("g_ampa_dend");
    SparseConnection * con_ext_inh_dend = new SparseConnection(inh_Poisson, pyr, w_inh_dend, p_inh, GABA);
    con_ext_inh_dend->set_target("g_gaba_dend");

    
    // define monitors
    double binSize_rate = 20.e-3; // ms
    auto seed_str = std::to_string(seed);
    SpikeMonitor * smon = new SpikeMonitor( pyr, sys->fn("ras_seed"+seed_str) );
    VoltageMonitor * vmon   = new VoltageMonitor( pyr, 0,   sys->fn("mem") );
    VoltageMonitor * vmon1   = new VoltageMonitor( pyr, 1,   sys->fn("mem1") );
    StateMonitor * smon_vd  = new StateMonitor( pyr, 0, "Vd", sys->fn("Vd") );
    StateMonitor * smon_vd1  = new StateMonitor( pyr, 1, "Vd", sys->fn("Vd1") );
    StateMonitor * smon_m   = new StateMonitor( pyr, 0, "thr", sys->fn("thr") );
    PopulationRateMonitor * pmon = new PopulationRateMonitor( pyr, sys->fn("prate"), binSize_rate );
    BurstRateMonitor * brmon = new BurstRateMonitor( pyr, sys->fn("brate_seed"+seed_str), binSize_rate);
    //StateMonitor * smon_ws  = new StateMonitor( pyr, 0, "wsoma", sys->fn("wsoma") );
    
    // run simulation
    logger->msg("Running ...",PROGRESS);
    
    // The alternating currents switched between a maximum of 175 pA and a minimum of -210 pA in the dendrites and between 230 pA and 130 pA in the soma. Dendritic input during 75e-3 s and somatic input during 50.e-3 seconds.
    const double max_dendritic_current = 40e-12; //90e-12;
    const double min_dendritic_current = -10e-12;
    const double max_somatic_current = 100.e-12; //650e-12;
    const double min_somatic_current = 0e-12;

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

    // alternating current to somata only (constant dendritic current):
    //double dendritic_current[1] = {min_dendritic_current, (min_dendritic_current + max_dendritic_current)/2., max_dendritic_current};
    /*
    double dendritic_current[1] = {max_dendritic_current};
    
    for (int j=0;j<1;j++){
        for (int i=0;i<10;i++)
        {
            //double somatic_current = min_somatic_current + (i + 1)*(max_somatic_current - min_somatic_current)/11;
            curr_inject_soma->set_all_currents(max_somatic_current/pyr[0].get_Cs());
            curr_inject_dend->set_all_currents(dendritic_current[j]/pyr[0].get_Cd());
            
            sys->run(segtime_mindend - small_overlap); //20 ms
            
            curr_inject_soma->set_all_currents(max_somatic_current/pyr[0].get_Cs());
            curr_inject_dend->set_all_currents(dendritic_current[j]/pyr[0].get_Cd());
            sys->run(segtime_maxsoma  - (segtime_mindend - small_overlap) ); //30 ms
            
            curr_inject_soma->set_all_currents(min_somatic_current/pyr[0].get_Cs());
            curr_inject_dend->set_all_currents(dendritic_current[j]/pyr[0].get_Cd());
            sys->run(segtime_maxdend + segtime_mindend - small_overlap - segtime_maxsoma); //40 ms
            
            curr_inject_soma->set_all_currents(min_somatic_current/pyr[0].get_Cs());
            curr_inject_dend->set_all_currents(dendritic_current[j]/pyr[0].get_Cd());
            sys->run(small_overlap); // 10 ms
        }
    }*/
    
    
    if (errcode)
        auryn_abort(errcode);
    
    logger->msg("Freeing ...",PROGRESS,true);
    auryn_free();
    return errcode;
}
