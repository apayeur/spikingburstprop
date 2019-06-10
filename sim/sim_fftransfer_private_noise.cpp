/***********************************************************
 Description: Compute the feedforward transfer function.
 No feedback propagation.
************************************************************/

#include "auryn.h"
#include <cmath>
#include <fstream>
#include <string>
#include "STPeTMConnection.h"
#include "BurstPoissonGroup.h"

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
    string dir = "../data/fftransfer/private-noise";
    string simname = "fftransfer";
    
    
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
    sys->set_master_seed(0);
    
    /**************************************************/
    /******          NEURON POPULATIONS       *********/
    /**************************************************/
    
    //-- First layer
    //  2-comp pyramidal neurons
    NeuronID number_of_neurons = 4000;
    NaudGroup* pyr1 = new NaudGroup(number_of_neurons);
    initialize_pyr_neurons(pyr1);
    
    //-- Second layer
    //  2-comp pyramidal neurons
    NaudGroup* pyr2 = new NaudGroup(number_of_neurons);
    initialize_pyr_neurons(pyr2);
    //  PV neurons
    NeuronID N_other_neuron = 1000;
    AdExGroup* pv2 = new AdExGroup(N_other_neuron);
    fix_parameters_pv_neurons(pv2);
    /*
    //-- Third layer
    //  2-comp pyramidal neurons
    BurstPoissonGroup* pyr3 = new BurstPoissonGroup(number_of_neurons);
    initialize_pyr_neurons(pyr3);
    //  PV neurons
    AdExGroup* pv3 = new AdExGroup(1000);
    fix_parameters_pv_neurons(pv3);
    */
    //-- Background
    // External populations of Poisson neurons (noise)
    // Pyr 1
    float poisson_rate = 5.;
    PoissonGroup* exc_Poisson1      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson1      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* exc_Poisson_dend1 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_dend1 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    
    // Pyr2
    PoissonGroup* exc_Poisson2      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson2      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* exc_Poisson_dend2 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_dend2 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    
    // PV
    PoissonGroup* exc_Poisson_pv    = new PoissonGroup(N_other_neuron, 100*poisson_rate);
    PoissonGroup* inh_Poisson_pv    = new PoissonGroup(N_other_neuron, 100*poisson_rate);
    
    /**************************************************/
    /******           CONNECTIVITY            *********/
    /**************************************************/
    
    //-- CONNECT BACKGROUND POISSON
    
    // External Poisson neurons -> Pyr1
    float w_epois_to_soma = 0.35;
    float ratio_ie_soma = 1.;
    IdentityConnection * epois_to_soma1 = new IdentityConnection(exc_Poisson1, pyr1, w_epois_to_soma, GLUT);
    epois_to_soma1->set_target("g_ampa");
    IdentityConnection * ipois_to_soma1 = new IdentityConnection(inh_Poisson1, pyr1, ratio_ie_soma*w_epois_to_soma, GABA);
    ipois_to_soma1->set_target("g_gaba");
    
    float w_epois_to_dend = 0.05;
    float ratio_ie_dend = 1.;
    IdentityConnection * epois_to_dend1 = new IdentityConnection(exc_Poisson_dend1, pyr1, w_epois_to_dend, GLUT);
    epois_to_dend1->set_target("g_ampa_dend");
    IdentityConnection * ipois_to_dend1 = new IdentityConnection(inh_Poisson_dend1, pyr1, ratio_ie_dend*w_epois_to_dend, GABA);
    ipois_to_dend1->set_target("g_gaba_dend");
    
    // External Poisson neurons -> Pyr2
    IdentityConnection * epois_to_soma2 = new IdentityConnection(exc_Poisson2, pyr2, w_epois_to_soma, GLUT);
    epois_to_soma2->set_target("g_ampa");
    IdentityConnection * ipois_to_soma2 = new IdentityConnection(inh_Poisson2, pyr2, ratio_ie_soma*w_epois_to_soma, GABA);
    ipois_to_soma2->set_target("g_gaba");
    IdentityConnection * epois_to_dend2 = new IdentityConnection(exc_Poisson_dend2, pyr2, w_epois_to_dend, GLUT);
    epois_to_dend2->set_target("g_ampa_dend");
    IdentityConnection * ipois_to_dend2 = new IdentityConnection(inh_Poisson_dend2, pyr2, ratio_ie_dend*w_epois_to_dend, GABA);
    ipois_to_dend2->set_target("g_gaba_dend");
    
    /*
    // External Poisson neurons -> Pyr3
    SparseConnection * con_ext_exc_soma3 = new SparseConnection(exc_Poisson, pyr3, w_exc, p_exc, GLUT);
    con_ext_exc_soma3->set_target("g_ampa");
    SparseConnection * con_ext_inh_soma3 = new SparseConnection(inh_Poisson, pyr3, w_inh, p_inh, GABA);
    con_ext_inh_soma3->set_target("g_gaba");
    SparseConnection * con_ext_exc_dend3 = new SparseConnection(exc_Poisson, pyr3, w_exc_dend, p_exc, GLUT);
    con_ext_exc_dend3->set_target("g_ampa_dend");
    SparseConnection * con_ext_inh_dend3 = new SparseConnection(inh_Poisson, pyr3, w_inh_dend, p_inh, GABA);
    con_ext_inh_dend3->set_target("g_gaba_dend");
    */
    
    // External Poisson neurons -> PV2
    float w_epois_to_pv = 0.75;
    float ratio_ie_pv = 1.5;
    IdentityConnection * epois_to_pv2 = new IdentityConnection(exc_Poisson_pv, pv2, w_epois_to_pv, GLUT);
    IdentityConnection * ipois_to_pv2 = new IdentityConnection(inh_Poisson_pv, pv2, w_epois_to_pv*ratio_ie_pv, GABA);
    
    /*
    // External Poisson neurons -> PV3
    SparseConnection * con_ext_exc_pv3 = new SparseConnection(exc_Poisson, pv3, w_extexc_to_pv, p_ext_to_pv, GLUT);
    SparseConnection * con_ext_inh_pv3 = new SparseConnection(inh_Poisson, pv3, w_extinh_to_pv, p_ext_to_pv, GABA);
    */
    
    //-- CONNECT FeedFORWARD
        // Pyr1 to pyr2 - STD
    float w_pyr1_to_pyr2 = 0.05; //0.014; 0.03 for EventBurstPoisson
    float p_pyr1_to_pyr2 = 0.05;
    STPeTMConnection * pyr1_to_pyr2 = new STPeTMConnection(pyr1, pyr2, w_pyr1_to_pyr2, p_pyr1_to_pyr2, GLUT);
    set_Depressing_connection(pyr1_to_pyr2);
    pyr1_to_pyr2->set_target("g_ampa");
    
        // Pyr1 to PV2 - STD
    float w_pyr1_to_pv2 = 0.05;//0.01
    float p_pyr1_to_pv2 = 0.05;
    STPeTMConnection * pyr1_to_pv2 = new STPeTMConnection(pyr1, pv2, w_pyr1_to_pv2, p_pyr1_to_pv2, GLUT);
    set_Depressing_connection(pyr1_to_pv2);

        // PV2 to Pyr2
    float w_pv2_to_pyr2 = 0.05; //0.05;
    float p_pv2_to_pyr2 = 0.05;
    STPeTMConnection * pv2_to_pyr2 = new STPeTMConnection(pv2, pyr2, w_pv2_to_pyr2, p_pv2_to_pyr2, GABA);
    pv2_to_pyr2->set_target("g_gaba");
    /*
    // Pyr2 to pyr3 - STD
    float w_pyr2_to_pyr3 = 0.05; //0.014;
    float p_pyr2_to_pyr3 = 0.05;
    STPeTMConnection * pyr2_to_pyr3 = new STPeTMConnection(pyr2, pyr3, w_pyr2_to_pyr3, p_pyr2_to_pyr3, GLUT);
    set_Depressing_connection(pyr2_to_pyr3);
    pyr2_to_pyr3->set_target("g_ampa");
    
    // Pyr2 to PV3 - STD
    float w_pyr2_to_pv3 = 0.01;
    float p_pyr2_to_pv3 = 0.05;
    STPeTMConnection * pyr2_to_pv3 = new STPeTMConnection(pyr2, pv3, w_pyr2_to_pv3, p_pyr2_to_pv3, GLUT);
    set_Depressing_connection(pyr2_to_pv3);
    
    // PV3 to Pyr3
    float w_pv3_to_pyr3 = 0.025; //0.05;
    float p_pv3_to_pyr3 = 0.05;
    STPeTMConnection * pv3_to_pyr3 = new STPeTMConnection(pv3, pyr3, w_pv3_to_pyr3, p_pv3_to_pyr3, GABA);
    pv3_to_pyr3->set_target("g_gaba");
    */


    /**************************************************/
    /******         CURRENT INJECTORS         *********/
    /**************************************************/
    CurrentInjector * curr_inject_soma1 = new CurrentInjector(pyr1, "mem");
    CurrentInjector * curr_inject_soma2 = new CurrentInjector(pyr2, "mem");
    CurrentInjector * curr_inject_pv2   = new CurrentInjector(pv2, "mem");
    //CurrentInjector * curr_inject_pv3    = new CurrentInjector(pv3, "mem");

    
    /**************************************************/
    /******              MONITORS             *********/
    /**************************************************/
    double binSize_rate = 20.e-3; // ms
    auto seed_str = std::to_string(seed);
    
    // Burst/event rate monitors
    BurstRateMonitor * brmon1 = new BurstRateMonitor( pyr1, sys->fn("brate1"), binSize_rate);
    BurstRateMonitor * brmon2 = new BurstRateMonitor( pyr2, sys->fn("brate2"), binSize_rate);
    //BurstRateMonitor * brmon3 = new BurstRateMonitor( pyr3, sys->fn("brate3"), binSize_rate);
    
    // Population activity monitors
    PopulationRateMonitor * pv2_activity   = new PopulationRateMonitor( pv2, sys->fn("pv2rate"), binSize_rate );
    //PopulationRateMonitor * pv3_activity   = new PopulationRateMonitor( pv3, sys->fn("pv3rate"), binSize_rate );
    
    // Voltage monitors
    VoltageMonitor *pyr1_mem = new VoltageMonitor(pyr1, 0,  sys->fn("pyr1_mem"), 1.e-3);
    VoltageMonitor *pyr2_mem = new VoltageMonitor(pyr2, 0,  sys->fn("pyr2_mem"), 1.e-3);
    //VoltageMonitor *pyr3_mem = new VoltageMonitor(pyr3, 0,  sys->fn("pyr3_mem"), 1.e-3);
    
    //Raster
    SpikeMonitor * smon1 = new SpikeMonitor( pyr1, sys->fn("ras1"));
    SpikeMonitor * smon2 = new SpikeMonitor( pyr2, sys->fn("ras2"));

    /**************************************************/
    /******             SIMULATION            *********/
    /**************************************************/
    logger->msg("Running ...",PROGRESS);
    
    const double max_somatic_current  = 500e-12;
    const double min_somatic_current  = -200e-12;
    const double somatic_current_incr = 100e-12;
    const double simtime = 1000e-3;
    
    double somatic_current = min_somatic_current;
    curr_inject_soma1->set_all_currents(somatic_current/pyr1[0].get_Cs());
    curr_inject_soma2->set_all_currents(-300e-12/pyr2[0].get_Cs());
    curr_inject_pv2->set_all_currents(205e-12/pv2[0].get_c_mem());
    //curr_inject_pv3->set_all_currents(205e-12/pv3[0].get_c_mem());

    while (somatic_current < max_somatic_current + somatic_current_incr/2)
    {
        curr_inject_soma1->set_all_currents(somatic_current/pyr1[0].get_Cd());
        somatic_current += somatic_current_incr;
        sys->run(simtime);
    }

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



