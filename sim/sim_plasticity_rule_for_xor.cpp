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
    int seed = 123;
    double w0 = 0.8;
    double d0 = 0.;
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
    const NeuronID number_of_neurons = 500;
    BurstPoissonGroup * main_neurons = new BurstPoissonGroup(number_of_neurons);
    main_neurons->syn_exc_soma->set_nmda_ampa_current_ampl_ratio(0.);
    main_neurons->syn_exc_dend->set_nmda_ampa_current_ampl_ratio(0.);
    //main_neurons->syn_exc_soma->set_ampa_nmda_ratio(1);
    //main_neurons->syn_exc_dend->set_ampa_nmda_ratio(100);
    
    // Input group
    BurstPoissonGroup * input = new BurstPoissonGroup(number_of_neurons);
    input->syn_exc_soma->set_nmda_ampa_current_ampl_ratio(0.);
    input->syn_exc_dend->set_nmda_ampa_current_ampl_ratio(0.);
    input->random_mem(-70e-3, 5.e-3);
    
    // Noise groups
        // dendritic noise for main neuron group
    float poisson_rate = 5.;
    PoissonGroup* exc_Poisson      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* exc_Poisson_dend = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_dend = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    
        // noise for input neuron group
    PoissonGroup* exc_Poisson_input      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_input      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* exc_Poisson_dend_input = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_dend_input = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    
    
    /**************************************************/
    /******           CONNECTIVITY            *********/
    /**************************************************/
    //-- Connect input group to main group's somata with plastic synapses
    const double we_soma = w0;
    const double sparseness = 0.05;
    const double learning_rate = 1e-3;
    const double max_weight = 1.0;
    
    BCPConnection * con_ext_soma;
    
    if (connect_type == "BCP") {
        con_ext_soma = new BCPConnection(input, main_neurons, we_soma,
                                         sparseness, learning_rate, max_weight, tau_pre, GLUT);
    }
    else if (connect_type == "EBCP"){
        con_ext_soma = new EBCPConnection(input, main_neurons, we_soma,
                                          sparseness, learning_rate, max_weight, tau_pre, GLUT);
    }
    else if (connect_type == "AdaptiveEBCP"){
        con_ext_soma = new AdaptiveEBCPConnection(input, main_neurons, we_soma,
                                                  sparseness, learning_rate, max_weight, tau_pre, GLUT);
        con_ext_soma->set_post_trace_event_tau(moving_average_time_constant);
        con_ext_soma->set_post_trace_burst_tau(moving_average_time_constant);

    }
    con_ext_soma->set_target("g_ampa");
    
    // noise
    // External Poisson neurons -> main neuron group
    float w_epois_to_soma = 0.35;
    float ratio_ie_soma = 1.;
    IdentityConnection * epois_to_soma = new IdentityConnection(exc_Poisson, main_neurons, w_epois_to_soma, GLUT);
    epois_to_soma->set_target("g_ampa");
    IdentityConnection * ipois_to_soma = new IdentityConnection(inh_Poisson, main_neurons, ratio_ie_soma*w_epois_to_soma, GABA);
    ipois_to_soma->set_target("g_gaba");
    
    float w_epois_to_dend = 0.05;
    float ratio_ie_dend = 1.;
    IdentityConnection * epois_to_dend = new IdentityConnection(exc_Poisson_dend, main_neurons, w_epois_to_dend, GLUT);
    epois_to_dend->set_target("g_ampa_dend");
    IdentityConnection * ipois_to_dend = new IdentityConnection(inh_Poisson_dend, main_neurons, ratio_ie_dend*w_epois_to_dend, GABA);
    ipois_to_dend->set_target("g_gaba_dend");
    
    // External Poisson neurons -> input neurons
    IdentityConnection * epois_to_soma_input = new IdentityConnection(exc_Poisson_input, input, w_epois_to_soma, GLUT);
    epois_to_soma_input->set_target("g_ampa");
    IdentityConnection * ipois_to_soma_input = new IdentityConnection(inh_Poisson_input, input, ratio_ie_soma*w_epois_to_soma, GABA);
    ipois_to_soma_input->set_target("g_gaba");
    IdentityConnection * epois_to_dend_input = new IdentityConnection(exc_Poisson_dend_input, input, w_epois_to_dend, GLUT);
    epois_to_dend_input->set_target("g_ampa_dend");
    IdentityConnection * ipois_to_dend_input = new IdentityConnection(inh_Poisson_dend_input, input, ratio_ie_dend*w_epois_to_dend, GABA);
    ipois_to_dend_input->set_target("g_gaba_dend");

    
    /**************************************************/
    /******              MONITORS             *********/
    /**************************************************/
    const double binsize = 5.;//moving_average_time_constant/5.;
    sys->set_online_rate_monitor_target(main_neurons);
    sys->set_online_rate_monitor_tau(binsize);
    SpikeMonitor * smon = new SpikeMonitor( main_neurons, sys->fn("ras") );
    PopulationRateMonitor * pmon = new PopulationRateMonitor( main_neurons, sys->fn("prate"), binsize );
    BurstRateMonitor * brmon = new BurstRateMonitor( main_neurons, sys->fn("brate"), binsize );
    BurstRateMonitor * brmon_input = new BurstRateMonitor( input, sys->fn("brate_input"), binsize );

    VoltageMonitor * vmon   = new VoltageMonitor( main_neurons, 0, sys->fn("mem"), 1e-3);
    vmon->record_for(10);
    StateMonitor * smon_vd  = new StateMonitor( main_neurons, 0, "Vd", sys->fn("Vd") );
    smon_vd->record_for(10);
    StateMonitor * smon_m   = new StateMonitor( main_neurons, 0, "thr", sys->fn("thr") );
    smon_m->record_for(10);
    StateMonitor * smon_ws  = new StateMonitor( main_neurons, 0, "wsoma", sys->fn("wsoma") );
    smon_ws->record_for(10);
    StateMonitor * smon_wd  = new StateMonitor( main_neurons, 0, "wdend", sys->fn("wdend") );
    smon_wd->record_for(10);
    
    //WeightMonitor * wmon = new WeightMonitor( con_ext_soma, sys->fn("consyn"), binsize);
    //wmon->add_equally_spaced(100);
    
    WeightSumMonitor * wsmon = new WeightSumMonitor( con_ext_soma, sys->fn("wsum") );
    
    StateMonitor * smon_tr_post_ev  = new StateMonitor( con_ext_soma->get_tr_event(), 0, sys->fn("trevent"), binsize);
    StateMonitor * smon_tr_post_b   = new StateMonitor( con_ext_soma->get_tr_burst(), 0, sys->fn("trburst"), binsize);
    
    /**************************************************/
    /******         CURRENT INJECTORS         *********/
    /**************************************************/
    CurrentInjector *curr_input_pop = new CurrentInjector(input, "Vd");
    CurrentInjector *curr_dend_exc = new CurrentInjector(main_neurons, "Vd");
    const double bkg_dend_curr = 400e-12;
    curr_dend_exc->set_all_currents(bkg_dend_curr/main_neurons[0].get_Cd());
    
    /**************************************************/
    /******             SIMULATION            *********/
    /**************************************************/
    logger->msg("Running ...",PROGRESS);
    
    // simulate
    double testtime = 10;
    sys->run(5000.);
    //con_ext_soma->set_all(0.25);
    //curr_dend_exc->set_all_currents(1.5*bkg_dend_curr/main_neurons[0].get_Cd());
    //sys->run(10.);
    //curr_dend_exc->set_all_currents(bkg_dend_curr/main_neurons[0].get_Cd());
    //sys->run(200.);

    /*
    //con_ext_soma->stdp_active = false;
    sys->run(moving_average_time_constant*3);
    
    con_ext_soma->stdp_active = true;
    sys->run(simtime);
    
    curr_dend_exc->set_all_currents(1.75*bkg_dend_curr/main_neurons[0].get_Cd());
    sys->run(simtime);
    
    curr_dend_exc->set_all_currents(bkg_dend_curr/main_neurons[0].get_Cd());
    sys->run(simtime);
    
    curr_dend_exc->set_all_currents(0*bkg_dend_curr/main_neurons[0].get_Cd());
    sys->run(simtime);
    
    curr_dend_exc->set_all_currents(bkg_dend_curr/main_neurons[0].get_Cd());
    sys->run(simtime);
    
    // test that plasticity is not affected by burst prob in input population
    curr_input_pop->set_all_currents(30.e-12/input[0].get_Cd());
    sys->run(simtime);
    */
    
    if (errcode)
        auryn_abort(errcode);
    
    logger->msg("Freeing ...",PROGRESS,true);
    auryn_free();
    return errcode;
}
