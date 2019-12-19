/***********************************************************
 Description: This program tests the burst-and-
 event-dependent plasticity rule.
 
 - Two populations of 2-compartment pyramidal neurons are
 each subjected to Poisson spiking noise.
 - Weights of the dendritic Poisson inputs onto the postsynaptic
 population modulate the burst fraction of that population.
 - The plastic weights are those that connect the presynaptic
 pyr population to the soma of the postsynaptic population.
************************************************************/

#include "auryn.h"
#include "EBCPConnection.h"
#include "BCPConnection.h"
#include "AdaptiveEBCPConnection.h"
#include "BurstPoissonGroup.h"

using namespace auryn;

namespace po = boost::program_options;

int main(int ac, char* av[])
{
    int errcode = 0;
    int seed = 1;
    double w0 = 0.8;
    double d0 = 0.1;
    double kappa = 1.0;
    double tau_pre = 20.e-3;
    double simtime = 1800;
    double moving_average_time_constant = 4.;
    std::string connect_type = "BCP";
    string simname = "plasticity_rule";
    std::string dir = ".";
    
    
    /**************************************************/
    /********              OPTIONS            *********/
    /**************************************************/
    try {
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help", "produce help message")
        ("d0", po::value<double>(), "dendritic input initial strength")
        ("w0", po::value<double>(), "total input initial strength")
        ("simtime", po::value<double>(), "simulation time")
        ("kappa", po::value<double>(), "poisson group rate")
        ("dir", po::value<string>(), "output directory")
        ("seed", po::value<int>(), "random seed")
        ("connect_type", po::value<string>(), "connection type: BCP or EBCP")
        ("tau_pre", po::value<double>(), "presynaptic trace time constant")
        ("alpha", po::value<double>(), "moving average time constant")
        ("sim_name", po::value<string>(), "simulation name")
        ;
        
        po::variables_map vm;
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);
        
        if (vm.count("help")) {
            std::cout << desc << "\n";
            return 1;
        }
        
        if (vm.count("d0")) {
            std::cout << "d0 set to "
            << vm["d0"].as<double>() << ".\n";
            d0 = vm["d0"].as<double>();
        }
        
        if (vm.count("w0")) {
            std::cout << "w0 set to "
            << vm["w0"].as<double>() << ".\n";
            w0 = vm["w0"].as<double>();
        }
        
        if (vm.count("kappa")) {
            std::cout << "kappa set to "
            << vm["kappa"].as<double>() << ".\n";
            kappa = vm["kappa"].as<double>();
        }
        
        if (vm.count("dir")) {
            std::cout << "dir set to "
            << vm["dir"].as<string>() << ".\n";
            dir = vm["dir"].as<string>();
        }
        
        if (vm.count("simtime")) {
            std::cout << "simtime set to "
            << vm["simtime"].as<double>() << ".\n";
            simtime = vm["simtime"].as<double>();
        }
        
        if (vm.count("seed")) {
            std::cout << "seed set to "
            << vm["seed"].as<int>() << ".\n";
            seed = vm["seed"].as<int>();
        }
        
        if (vm.count("connect_type")) {
            std::cout << "connection type set to "
            << vm["connect_type"].as<string>() << ".\n";
            connect_type = vm["connect_type"].as<string>();
        }
        
        if (vm.count("tau_pre")) {
            std::cout << "presynaptic trace time constant set to "
            << vm["tau_pre"].as<double>() << ".\n";
            tau_pre = vm["tau_pre"].as<double>();
        }
        
        if (vm.count("alpha")) {
            std::cout << "moving average time constant set to "
            << vm["alpha"].as<double>() << ".\n";
            moving_average_time_constant = vm["alpha"].as<double>();
        }
        
        if (vm.count("sim_name")) {
            std::cout << "simulation name set to "
            << vm["sim_name"].as<string>() << ".\n";
            simname = vm["sim_name"].as<string>();
        }
    }
    catch(std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        std::cerr << "Exception of unknown type!\n";
    }
    
    
    // INITIALIZE AURYN AND SET THE SEED FOR REPRODUCIBILITY
    auryn_init( ac, av, dir, simname);
    sys->set_master_seed(seed);
    
    
    /**************************************************/
    /******          NEURON POPULATIONS       *********/
    /**************************************************/
    
    // Main neuron group
    NaudGroup * neurons_exc = new NaudGroup(1000);
    neurons_exc->syn_exc_soma->set_nmda_ampa_current_ampl_ratio(0.);
    neurons_exc->syn_exc_dend->set_nmda_ampa_current_ampl_ratio(0.);
    //neurons_exc->syn_exc_soma->set_ampa_nmda_ratio(1);
    //neurons_exc->syn_exc_dend->set_ampa_nmda_ratio(100);
    
    // Input group
    //PoissonGroup * poisson_soma = new PoissonGroup(4000, kappa);
    NaudGroup * input = new NaudGroup(1000);
    input->syn_exc_soma->set_nmda_ampa_current_ampl_ratio(0.);
    input->syn_exc_dend->set_nmda_ampa_current_ampl_ratio(0.);
    input->random_mem(-70e-3, 5.e-3);
    
    // Noise groups
        // dendritic noise for main neuron group
    PoissonGroup * poisson_dend = new PoissonGroup(100, kappa);
    
        // noise for input neuron group
    NeuronID nb_exc_Poisson_neurons = 25000;
    NeuronID nb_inh_Poisson_neurons = nb_exc_Poisson_neurons;
    float poisson_rate = 5.;
    PoissonGroup* exc_Poisson = new PoissonGroup(nb_exc_Poisson_neurons, poisson_rate);
    PoissonGroup* inh_Poisson = new PoissonGroup(nb_inh_Poisson_neurons, poisson_rate);
    
    
    /**************************************************/
    /******           CONNECTIVITY            *********/
    /**************************************************/
    //-- Connect input group to main group's somata with plastic synapses
    const double we_soma = w0;
    const double we_dend = d0;
    const double sparseness = 0.2;
    const double learning_rate = 5e-3;
    const double max_weight = 1.0;
    
    BCPConnection * con_ext_soma;
    
    if (connect_type == "BCP") {
        std::cout << "WARNING! CHANGING TAU_PRE NOT SUPPORTED!!! tau_pre = 20 ms"<<std::endl;
        con_ext_soma = new BCPConnection(input, neurons_exc, we_soma,
                                         sparseness, learning_rate, max_weight, GLUT);
    }
    else if (connect_type == "EBCP"){
        con_ext_soma = new EBCPConnection(input, neurons_exc, we_soma,
                                          sparseness, learning_rate, max_weight, tau_pre, GLUT);
    }
    else if (connect_type == "AdaptiveEBCP"){
        con_ext_soma = new AdaptiveEBCPConnection(input, neurons_exc, we_soma,
                                                  sparseness, learning_rate, max_weight, tau_pre, GLUT);
        con_ext_soma->set_post_trace_event_tau(moving_average_time_constant);
        con_ext_soma->set_post_trace_burst_tau(moving_average_time_constant);

    }
    
    con_ext_soma->set_target("g_ampa");
    
    
    //-- Connect first Poisson noise group to main group's dendrites
    SparseConnection * con_ext_dend = new SparseConnection(poisson_dend, neurons_exc, we_dend, sparseness);
    con_ext_dend->set_target("g_ampa_dend");
    
    //-- Connect noise to input population
    int number_of_ext_conn[2] = {100, 30}; //exc, inh
    
    float p_exc = float(number_of_ext_conn[0])/float(nb_exc_Poisson_neurons);
    float p_inh = float(number_of_ext_conn[1])/float(nb_inh_Poisson_neurons);
    
    const float w_exc = 0.17;   // conductance amplitude in units of leak conductance
    const float w_inh = 0.26;
    
    SparseConnection * con_ext_exc_soma1 = new SparseConnection(exc_Poisson, input, w_exc, p_exc, GLUT);
    con_ext_exc_soma1->set_target("g_ampa");
    SparseConnection * con_ext_inh_soma1 = new SparseConnection(inh_Poisson, input, w_inh, p_inh, GABA);
    con_ext_inh_soma1->set_target("g_gaba");
    
    const float w_exc_dend = 0.0425;        // conductance amplitude in units of leak conductance
    const float w_inh_dend = 0.065;
    
    SparseConnection * con_ext_exc_dend1 = new SparseConnection(exc_Poisson, input, w_exc_dend, p_exc, GLUT);
    con_ext_exc_dend1->set_target("g_ampa_dend");
    SparseConnection * con_ext_inh_dend1 = new SparseConnection(inh_Poisson, input, w_inh_dend, p_inh, GABA);
    con_ext_inh_dend1->set_target("g_gaba_dend");
    
    /**************************************************/
    /******              MONITORS             *********/
    /**************************************************/
    const double binsize = 4.; //std::max(moving_average_time_constant/5., 5.);
    sys->set_online_rate_monitor_target(neurons_exc);
    sys->set_online_rate_monitor_tau(binsize);
    SpikeMonitor * smon = new SpikeMonitor( neurons_exc, sys->fn("ras") );
    PopulationRateMonitor * pmon = new PopulationRateMonitor( neurons_exc, sys->fn("prate"), binsize );
    BurstRateMonitor * brmon = new BurstRateMonitor( neurons_exc, sys->fn("brate"), binsize );
    BurstRateMonitor * brmon_input = new BurstRateMonitor( input, sys->fn("brate_input"), binsize );

    VoltageMonitor * vmon   = new VoltageMonitor( neurons_exc, 0, sys->fn("mem"), 1e-3);
    StateMonitor * smon_vd  = new StateMonitor( neurons_exc, 0, "Vd", sys->fn("Vd") );
    StateMonitor * smon_m   = new StateMonitor( neurons_exc, 0, "thr", sys->fn("thr") );
    smon_m->record_for(10);
    StateMonitor * smon_ws  = new StateMonitor( neurons_exc, 0, "wsoma", sys->fn("wsoma") );
    smon_ws->record_for(10);
    StateMonitor * smon_wd  = new StateMonitor( neurons_exc, 0, "wdend", sys->fn("wdend") );
    smon_wd->record_for(10);
    
    WeightMonitor * wmon = new WeightMonitor( con_ext_soma, sys->fn("consyn"), binsize);
    wmon->add_equally_spaced(100);
    
    WeightSumMonitor * wsmon = new WeightSumMonitor( con_ext_soma, sys->fn("wsum"), binsize );
    
    // Monitors for  estimating burst probability
    std::vector< StateMonitor* > smon_tr_post_ev;
    std::vector< StateMonitor* > smon_tr_post_b;
    
    for (int i=0;i<50;i++){
        StateMonitor * ev = new StateMonitor( con_ext_soma->get_tr_event(), i, sys->fn(simname,i,"trevent"), binsize);
        smon_tr_post_ev.push_back(ev);
        StateMonitor * b = new StateMonitor( con_ext_soma->get_tr_burst(), i, sys->fn(simname,i,"trburst"), binsize);
        smon_tr_post_b.push_back(b);
        
    }
    /**************************************************/
    /******         CURRENT INJECTORS         *********/
    /**************************************************/
    CurrentInjector *curr_input_pop = new CurrentInjector(input, "Vd");
    CurrentInjector *curr_dend_main = new CurrentInjector(neurons_exc, "Vd");
    curr_dend_main->set_all_currents(50.e-12/neurons_exc[0].get_Cd());
    
    /**************************************************/
    /******             SIMULATION            *********/
    /**************************************************/
    logger->msg("Running ...",PROGRESS);
    
    // simulate
    /*

    double testtime = 10;
    con_ext_soma->stdp_active = false;
    sys->run(5*moving_average_time_constant);
    con_ext_soma->stdp_active = true;
    sys->run(4000.);
     */
    /*
    vmon->record_for(20);
    smon_vd->record_for(20);
    sys->run(100.);
    vmon->record_for(20);
    smon_vd->record_for(20);
    sys->run(100.);
    vmon->record_for(20);
    smon_vd->record_for(20);
    sys->run(100.);
    vmon->record_for(20);
    smon_vd->record_for(20);
    sys->run(100.);
    smon_vd->record_for(20);
    vmon->record_for(20);
    sys->run(100.);
*/
    double testtime = 10;
    
    con_ext_soma->stdp_active = false;
    con_ext_dend->set_all(1.0*we_dend);
    sys->run(moving_average_time_constant*4);
    
    con_ext_soma->stdp_active = true;
    sys->run(simtime);

    con_ext_dend->set_all(1.5*we_dend);
    sys->run(simtime);

    con_ext_dend->set_all(1.0*we_dend);
    sys->run(simtime);
    
    con_ext_dend->set_all(0.55*we_dend);
    sys->run(simtime);
    
    con_ext_dend->set_all(1.0*we_dend);
    sys->run(simtime);
    
    // test that plasticity is not affected by burst prob in input population
    curr_input_pop->set_all_currents(40.e-12/input[0].get_Cd());
    sys->run(simtime);
    
    if (errcode)
        auryn_abort(errcode);
    
    logger->msg("Freeing ...",PROGRESS,true);
    auryn_free();
    return errcode;
}
