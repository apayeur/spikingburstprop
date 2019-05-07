/***********************************************************
 Description: Simulation of feedforward propagation
 of events and backward propagation of burst probabilities.
************************************************************/

#include "auryn.h"
#include <cmath>
#include <fstream>
#include <string>
#include "STPeTMConnection.h"


using namespace auryn;

namespace po = boost::program_options;

void fix_parameters_pv_neurons(AdExGroup* pv);
void initialize_pyr_neurons(NaudGroup* pyr);
void fix_parameters_som_neurons(AdExGroup* som);
void set_Facilitating_connection(STPeTMConnection * pre_to_post);
void set_Depressing_connection(STPeTMConnection * pre_to_post);


int main(int ac, char* av[])
{
    int errcode = 0;
    char strbuf [255];
    unsigned int seed = 1;
    string dir = "./";
    string simname = "propagation";
    
    /**************************************************/
    /********              OPTIONS            *********/
    /**************************************************/
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
    
    // INITIALIZE AURYN
    auryn_init( ac, av, dir, simname );

    /**************************************************/
    /******          NEURON POPULATIONS       *********/
    /**************************************************/
    
    // First layer of 2-comp pyramidal neurons
    NeuronID number_of_neurons = 4000;
    NaudGroup* pyr1 = new NaudGroup(number_of_neurons);
    initialize_pyr_neurons(pyr1);
    
    // Second layer of 2-comp pyramidal neurons
    NaudGroup* pyr2 = new NaudGroup(number_of_neurons);
    initialize_pyr_neurons(pyr2);
    
    // PV neurons
    AdExGroup* pv = new AdExGroup(1000);
    fix_parameters_pv_neurons(pv);
    
    // SOM neurons
    AdExGroup* som = new AdExGroup(1000);
    fix_parameters_som_neurons(som);
    
    // External population of Poisson neurons
    NeuronID nb_exc_Poisson_neurons = 25000;
    NeuronID nb_inh_Poisson_neurons = nb_exc_Poisson_neurons;
    float poisson_rate = 5.;
    PoissonGroup* exc_Poisson = new PoissonGroup(nb_exc_Poisson_neurons, poisson_rate);
    PoissonGroup* inh_Poisson = new PoissonGroup(nb_inh_Poisson_neurons, poisson_rate);
    exc_Poisson->seed(seed);
    inh_Poisson->seed(seed);
    
    
    /**************************************************/
    /******           CONNECTIVITY            *********/
    /**************************************************/
    
    //-- CONNECT BACKGROUND POISSON
    
    // External Poisson neurons -> Pyr1
    int number_of_ext_conn[2] = {100, 30}; //exc, inh
    
    float p_exc = float(number_of_ext_conn[0])/float(nb_exc_Poisson_neurons);
    float p_inh = float(number_of_ext_conn[1])/float(nb_inh_Poisson_neurons);
    
    const float w_exc = 0.17;   // conductance amplitude in units of leak conductance
    const float w_inh = 0.26;
    
    SparseConnection * con_ext_exc_soma1 = new SparseConnection(exc_Poisson, pyr1, w_exc, p_exc, GLUT);
    con_ext_exc_soma1->set_target("g_ampa");
    SparseConnection * con_ext_inh_soma1 = new SparseConnection(inh_Poisson, pyr1, w_inh, p_inh, GABA);
    con_ext_inh_soma1->set_target("g_gaba");
    
    const float w_exc_dend = 0.0425;        // conductance amplitude in units of leak conductance
    const float w_inh_dend = 0.065;
    
    SparseConnection * con_ext_exc_dend1 = new SparseConnection(exc_Poisson, pyr1, w_exc_dend, p_exc, GLUT);
    con_ext_exc_dend1->set_target("g_ampa_dend");
    SparseConnection * con_ext_inh_dend1 = new SparseConnection(inh_Poisson, pyr1, w_inh_dend, p_inh, GABA);
    con_ext_inh_dend1->set_target("g_gaba_dend");
    
    // External Poisson neurons -> Pyr2
    SparseConnection * con_ext_exc_soma2 = new SparseConnection(exc_Poisson, pyr2, w_exc, p_exc, GLUT);
    con_ext_exc_soma2->set_target("g_ampa");
    SparseConnection * con_ext_inh_soma2 = new SparseConnection(inh_Poisson, pyr2, w_inh, p_inh, GABA);
    con_ext_inh_soma2->set_target("g_gaba");
    SparseConnection * con_ext_exc_dend2 = new SparseConnection(exc_Poisson, pyr2, w_exc_dend, p_exc, GLUT);
    con_ext_exc_dend2->set_target("g_ampa_dend");
    SparseConnection * con_ext_inh_dend2 = new SparseConnection(inh_Poisson, pyr2, w_inh_dend, p_inh, GABA);
    con_ext_inh_dend2->set_target("g_gaba_dend");
    
    // External Poisson neurons -> SOM
    float p_ext_to_som = 500./float(nb_exc_Poisson_neurons);
    SparseConnection * con_ext_exc_som = new SparseConnection(exc_Poisson, som, 0.08, p_ext_to_som, GLUT);
    SparseConnection * con_ext_inh_som = new SparseConnection(inh_Poisson, som, 0.1, p_ext_to_som, GABA);
    
    // External Poisson neurons -> PV
    float p_ext_to_pv = 100./float(nb_exc_Poisson_neurons);
    float w_extexc_to_pv = 0.12;
    float w_extinh_to_pv = 0.15;
    SparseConnection * con_ext_exc_pv = new SparseConnection(exc_Poisson, pv, w_extexc_to_pv, p_ext_to_pv, GLUT);
    SparseConnection * con_ext_inh_pv = new SparseConnection(inh_Poisson, pv, w_extinh_to_pv, p_ext_to_pv, GABA);
    
    
    //-- CONNECT FeedFORWARD
        // Pyr1 to pyr2 - STD
    float w_pyr1_to_pyr2 = 0.015; //0.017;
    float p_pyr1_to_pyr2 = 0.05;
    STPeTMConnection * pyr1_to_pyr2 = new STPeTMConnection(pyr1, pyr2, w_pyr1_to_pyr2, p_pyr1_to_pyr2, GLUT);
    set_Depressing_connection(pyr1_to_pyr2);
    pyr1_to_pyr2->set_target("g_ampa");
    
        // Pyr1 to PV - STD
    float w_pyr1_to_pv = 0.01;
    float p_pyr1_to_pv = 0.05;
    STPeTMConnection * pyr1_to_pv = new STPeTMConnection(pyr1, pv, w_pyr1_to_pv, p_pyr1_to_pv, GLUT);
    set_Depressing_connection(pyr1_to_pv);

        // PV to Pyr2
    float w_pv_to_pyr2 = 0.0; //0.04;
    float p_pv_to_pyr2 = 0.05;
    SparseConnection * pv_to_pyr2 = new SparseConnection(pv, pyr2, w_pv_to_pyr2, p_pv_to_pyr2, GABA);
    pv_to_pyr2->set_target("g_gaba");

    
    //-- CONNECT FeedBACK
        // Pyr2 to pyr1 - STF
    float w_pyr2_to_pyr1 = 0.2;//0.2;
    float p_pyr2_to_pyr1 = 0.05;
    STPeTMConnection * pyr2_to_pyr1 = new STPeTMConnection(pyr2, pyr1, w_pyr2_to_pyr1, p_pyr2_to_pyr1, GLUT);
    set_Facilitating_connection(pyr2_to_pyr1);
    pyr2_to_pyr1->set_target("g_ampa_dend");
    
        // Pyr2 to SOM - STD
    float w_pyr2_to_som = 0.02;
    float p_pyr2_to_som = 0.05;
    //SparseConnection * pyr2_to_som = new SparseConnection(pyr2, som, w_pyr2_to_som, p_pyr2_to_som, GLUT);
    
    STPeTMConnection * pyr2_to_som = new STPeTMConnection(pyr2, som, w_pyr2_to_som, p_pyr2_to_som, GLUT);
    set_Depressing_connection(pyr2_to_som);
    
        // SOM to pyr1
    float w_som_to_pyr1 = 0.03;//0.025;
    float p_som_to_pyr1 = 0.05;
    SparseConnection * som_to_pyr1 = new SparseConnection(som, pyr1, w_som_to_pyr1, p_som_to_pyr1, GABA);
    som_to_pyr1->set_target("g_gaba_dend");
    
    // -- OTHER CONNECTIONS
        // Pyr to Pyr - fake inh STF
    float w_pyr_to_pyr = 0.02; //0.02; //original: 0.2e-9/5.e-9;
    float p_pyr_to_pyr = 0.05;
    STPeTMConnection * pyr1_to_pyr1 = new STPeTMConnection(pyr1, pyr1, w_pyr_to_pyr, p_pyr_to_pyr, GABA);
    set_Facilitating_connection(pyr1_to_pyr1);
    pyr1_to_pyr1->set_target("g_gaba_dend");
    
    STPeTMConnection * pyr2_to_pyr2 = new STPeTMConnection(pyr2, pyr2, w_pyr_to_pyr, p_pyr_to_pyr, GABA);
    set_Facilitating_connection(pyr2_to_pyr2);
    pyr2_to_pyr2->set_target("g_gaba_dend");

    
    /**************************************************/
    /******         CURRENT INJECTORS         *********/
    /**************************************************/
    CurrentInjector * curr_inject_soma1 = new CurrentInjector(pyr1, "mem");
    CurrentInjector * curr_inject_dend1 = new CurrentInjector(pyr1, "Vd");
    CurrentInjector * curr_inject_soma2 = new CurrentInjector(pyr2, "mem");
    CurrentInjector * curr_inject_dend2 = new CurrentInjector(pyr2, "Vd");
    CurrentInjector * curr_inject_som   = new CurrentInjector(som, "mem");
    CurrentInjector * curr_inject_pv    = new CurrentInjector(pv, "mem");
    
    /**************************************************/
    /******              MONITORS             *********/
    /**************************************************/
    double binSize_rate = 20.e-3; // ms
    
    // Voltage monitors
    VoltageMonitor * vmon_som       = new VoltageMonitor( som, 0, sys->fn("memsom") );
    VoltageMonitor * vmon_pv        = new VoltageMonitor( pv, 0, sys->fn("mempv") );
    VoltageMonitor * vmon_pyrsoma1  = new VoltageMonitor( pyr1, 0, sys->fn("memsoma1") );
    VoltageMonitor * vmon_pyrsoma2  = new VoltageMonitor( pyr2, 0, sys->fn("memsoma2") );
    StateMonitor * smon_vd          = new StateMonitor( pyr1, 0, "Vd", sys->fn("Vd") );
    StateMonitor * smon_vd2         = new StateMonitor( pyr2, 0, "Vd", sys->fn("Vd2") );
    
    // Burst and event rate monitors
    // To plot burst rate in gnuplot, u 1:2; to plot event rate, u 1:3
    auto seed_str = std::to_string(seed);
    BurstRateMonitor * brmon1 = new BurstRateMonitor( pyr1, sys->fn("brate1_seed"+seed_str), binSize_rate);
    BurstRateMonitor * brmon2 = new BurstRateMonitor( pyr2, sys->fn("brate2_seed"+seed_str), binSize_rate);
    
    // Raster plots
    SpikeMonitor * smon1 = new SpikeMonitor( pyr1, sys->fn("ras1_seed"+seed_str));
    SpikeMonitor * smon2 = new SpikeMonitor( pyr2, sys->fn("ras2_seed"+seed_str));
    SpikeMonitor * smon_pv = new SpikeMonitor( pv, sys->fn("raspv_seed"+seed_str));
    SpikeMonitor * smon_som = new SpikeMonitor( som, sys->fn("rassom_seed"+seed_str));
    

    /**************************************************/
    /******             SIMULATION            *********/
    /**************************************************/
    logger->msg("Running ...",PROGRESS);
    
    // The alternating currents switched between a maximum and a minimum in both the dendrites and the somas.
    const double max_dendritic_current = 50e-12;//90e-12;
    const double min_dendritic_current = 0e-12;
    const double max_somatic_current = 100e-12;
    const double min_somatic_current = 0.; //50e-12;
    
    const double simtime = 1000e-3;
    const double period = 200e-3;
    const double segtime_maxsoma = period/2.;
    const double segtime_minsoma = period/2;
    const double segtime_maxdend = 0.7*period;
    const double segtime_mindend = 0.3*period;
    const double small_overlap = 0.1*period;
    
    // Burn-in period (i.e. relaxation) before alternating stimuation
    curr_inject_soma1->set_all_currents(min_somatic_current/pyr1[0].get_Cs());
    curr_inject_dend2->set_all_currents(min_dendritic_current/pyr2[0].get_Cd());
    curr_inject_som->set_all_currents(0e-12/som[0].get_c_mem());
    curr_inject_pv->set_all_currents(205e-12/pv[0].get_c_mem());
    sys->run(simtime);
    
    // Main simulation
    for (int i=0;i<10;i++)
    {
        curr_inject_soma1->set_all_currents(max_somatic_current/pyr1[0].get_Cs());
        curr_inject_dend2->set_all_currents(min_dendritic_current/pyr2[0].get_Cd());
        sys->run(segtime_mindend - small_overlap);
        
        curr_inject_soma1->set_all_currents(max_somatic_current/pyr1[0].get_Cs());
        curr_inject_dend2->set_all_currents(max_dendritic_current/pyr2[0].get_Cd());
        sys->run(segtime_maxsoma  - (segtime_mindend - small_overlap) );
        
        curr_inject_soma1->set_all_currents(min_somatic_current/pyr1[0].get_Cs());
        curr_inject_dend2->set_all_currents(max_dendritic_current/pyr2[0].get_Cd());
        sys->run(segtime_maxdend + segtime_mindend - small_overlap - segtime_maxsoma);
        
        curr_inject_soma1->set_all_currents(min_somatic_current/pyr1[0].get_Cs());
        curr_inject_dend2->set_all_currents(min_dendritic_current/pyr2[0].get_Cd());
        sys->run(small_overlap);
    }
    /*
     // alternating current to dendrites only (constant somatic current)
     double somatic_current = min_somatic_current;
     for (int i=0;i<10;i++)
     {
     curr_inject_soma->set_all_currents(somatic_current/pyr[0].get_Cs());
     curr_inject_dend->set_all_currents(min_dendritic_current/pyr[0].get_Cd());
     sys->run(segtime_mindend - small_overlap);
     
     curr_inject_soma->set_all_currents(somatic_current/pyr[0].get_Cs());
     curr_inject_dend->set_all_currents(max_dendritic_current/pyr[0].get_Cd());
     sys->run(segtime_maxsoma  - (segtime_mindend - small_overlap) );
     
     curr_inject_soma->set_all_currents(somatic_current/pyr[0].get_Cs());
     curr_inject_dend->set_all_currents(max_dendritic_current/pyr[0].get_Cd());
     sys->run(segtime_maxdend + segtime_mindend - small_overlap - segtime_maxsoma);
     
     curr_inject_soma->set_all_currents(somatic_current/pyr[0].get_Cs());
     curr_inject_dend->set_all_currents(min_dendritic_current/pyr[0].get_Cd());
     sys->run(small_overlap);
     }*/
    
    // alternating current to somata only (constant dendritic current):
    //double dendritic_current[1] = {min_dendritic_current, (min_dendritic_current + max_dendritic_current)/2., max_dendritic_current};
    /* double dendritic_current[1] = {max_dendritic_current};
     
     for (int j=0;j<1;j++){
     for (int i=0;i<10;i++)
     {
     //double somatic_current = min_somatic_current + (i + 1)*(max_somatic_current - min_somatic_current)/11;
     curr_inject_soma1->set_all_currents(max_somatic_current/pyr1[0].get_Cs());
     curr_inject_dend1->set_all_currents(dendritic_current[j]/pyr1[0].get_Cd());
     curr_inject_dend2->set_all_currents(dendritic_current[j]/pyr2[0].get_Cd());
     
     sys->run(segtime_mindend - small_overlap); //20 ms
     
     curr_inject_soma1->set_all_currents(max_somatic_current/pyr1[0].get_Cs());
     curr_inject_dend1->set_all_currents(dendritic_current[j]/pyr1[0].get_Cd());
     curr_inject_dend2->set_all_currents(dendritic_current[j]/pyr2[0].get_Cd());
     sys->run(segtime_maxsoma  - (segtime_mindend - small_overlap) ); //30 ms
     
     curr_inject_soma1->set_all_currents(min_somatic_current/pyr1[0].get_Cs());
     curr_inject_dend1->set_all_currents(dendritic_current[j]/pyr1[0].get_Cd());
     curr_inject_dend2->set_all_currents(dendritic_current[j]/pyr2[0].get_Cd());
     sys->run(segtime_maxdend + segtime_mindend - small_overlap - segtime_maxsoma); //40 ms
     
     curr_inject_soma1->set_all_currents(min_somatic_current/pyr1[0].get_Cs());
     curr_inject_dend1->set_all_currents(dendritic_current[j]/pyr1[0].get_Cd());
     curr_inject_dend2->set_all_currents(dendritic_current[j]/pyr2[0].get_Cd());
     sys->run(small_overlap); // 10 ms
     }
     }*/
    
    if (errcode)
        auryn_abort(errcode);
    
    logger->msg("Freeing ...",PROGRESS,true);
    auryn_free();
    return errcode;
}

void fix_parameters_pv_neurons(AdExGroup* pv) {
    AurynFloat rheobase_pv = 200e-12;
    AurynFloat gL_pv = 10e-9;
    AurynFloat taum_pv = 10e-3;
    AurynFloat EL_pv = -70.e-3;
    AurynFloat DeltaT_pv = 2.e-3;
    
    pv->set_refractory_period(1.e-3);
    pv->set_g_leak(gL_pv);
    pv->set_delta_t(DeltaT_pv);
    pv->set_c_mem(gL_pv*taum_pv);
    pv->set_a(0.);
    pv->set_b(0.);
    pv->set_e_rest(EL_pv);
    pv->set_e_reset(EL_pv + 15.e-3);
    pv->set_e_thr(rheobase_pv/gL_pv + EL_pv + DeltaT_pv);
    pv->set_tau_w(5.e-3); //no importance
    
    pv->random_mem(EL_pv, 5.e-3);
    pv->set_tau_ampa(5.e-3);
    pv->set_tau_gaba(10.e-3);
}

void fix_parameters_som_neurons(AdExGroup* som) {
    AurynFloat rheobase_som = 25e-12;
    AurynFloat gL_som = 5.e-9;
    AurynFloat taum_som = 20e-3;
    AurynFloat tauw_som = 500e-3;
    AurynFloat EL_som = -70.e-3;
    AurynFloat DeltaT_som = 4.e-3;
    AurynFloat a_som = 0.5e-9;
    AurynFloat b_som = 10.e-12;
    AurynFloat Vreset_som = EL_som + 5.e-3;
    AurynFloat tauref_som = 2.e-3;
    
    som->set_refractory_period(tauref_som);
    som->set_g_leak(gL_som);
    som->set_delta_t(DeltaT_som);
    som->set_c_mem(gL_som*taum_som);
    som->set_a(a_som);
    som->set_b(b_som);
    som->set_e_rest(EL_som);
    som->set_e_reset(Vreset_som);
    som->set_e_thr(rheobase_som/(gL_som+a_som) + EL_som + DeltaT_som - DeltaT_som*log(1. + a_som/gL_som));
    som->set_tau_w(tauw_som);
    
    som->random_mem(EL_som, 5.e-3);
    som->set_tau_ampa(5.e-3);
    som->set_tau_gaba(10.e-3);
}


void initialize_pyr_neurons(NaudGroup* pyr) {
    pyr->syn_exc_soma->set_nmda_ampa_current_ampl_ratio(0.);
    pyr->syn_exc_dend->set_nmda_ampa_current_ampl_ratio(0.);
    pyr->random_mem(-70e-3, 5.e-3);
}

void set_Facilitating_connection(STPeTMConnection * pre_to_post){
    pre_to_post->set_tau_d(100e-3);
    pre_to_post->set_tau_f(100e-3);
    pre_to_post->set_ujump(0.02);
    pre_to_post->set_f(0.1);
}

void set_Depressing_connection(STPeTMConnection * pre_to_post){
    pre_to_post->set_tau_d(20.e-3);
    pre_to_post->set_tau_f(1.);
    pre_to_post->set_ujump(0.9);
    pre_to_post->set_f(0.1);
}


