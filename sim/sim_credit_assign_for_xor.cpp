/***********************************************************
 Description: Example of "credit assignment" and error
 generation with BurstPoisson neurons and Transmit connections.
 ************************************************************/

#include "auryn.h"
#include <cmath>
#include <fstream>
#include <string>
#include "FileCurrentInjector.h"
#include "BurstPoissonGroup.h"
#include "TransmitBurstConnection.h"
#include "TransmitEventConnection.h"
#include "AdaptiveEBCPConnection.h"

using namespace auryn;

namespace po = boost::program_options;

void initialize_pyr_neurons(NaudGroup* pyr);


int main(int ac, char* av[])
{
    int errcode = 0;
    char strbuf [255];
    unsigned int seed = 1;
    string dir = "./";
    string simname = "credit_assign";
    const NeuronID number_of_neurons = 1000;
    const NeuronID N_other_neuron = number_of_neurons/4;
    
    float w_pyr1_to_pyr2_exc = 0.07*4000/number_of_neurons;
    float w_pyr1_to_pyr2_inh = 0.03*4000/number_of_neurons;

    float w_pyr2_to_pyr1_exc = 0.09*4000/number_of_neurons;
    float w_pyr2_to_pyr1_inh = 0.03*4000/number_of_neurons;

    /**************************************************/
    /********              OPTIONS            *********/
    /**************************************************/
    try {
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help", "produce help message")
        ("dir",  po::value<string>(), "output directory")
        ("seed", po::value<int>(),    "random seed")
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
    sys->set_master_seed(seed);
    
    /**************************************************/
    /******          NEURON POPULATIONS       *********/
    /**************************************************/
    
    // First layer of 2-comp pyramidal neurons
    BurstPoissonGroup* pyr1 = new BurstPoissonGroup(number_of_neurons);
    initialize_pyr_neurons(pyr1);
    
    // Second layer of 2-comp pyramidal neurons
    BurstPoissonGroup* pyr2 = new BurstPoissonGroup(number_of_neurons);
    initialize_pyr_neurons(pyr2);
    
    // External populations of Poisson neurons (noise)
    // Pyr 1
    float poisson_rate = 5.;
    PoissonGroup* exc_Poisson1      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson1      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* exc_Poisson_dend1 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_dend1 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    
    // Pyr2
    PoissonGroup* exc_Poisson2      = new PoissonGroup(number_of_neurons, 50*poisson_rate);
    PoissonGroup* inh_Poisson2      = new PoissonGroup(number_of_neurons, 50*poisson_rate);
    PoissonGroup* exc_Poisson_dend2 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_dend2 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    
    /**************************************************/
    /******           CONNECTIVITY            *********/
    /**************************************************/
    
    //----- CONNECT BACKGROUND POISSON -----//
    
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
    
    
    //------ CONNECT FeedFORWARD ------//
    // Pyr1 to pyr2 - Excitation
    float p_pyr1_to_pyr2 = 0.05; //0.05
    const float learning_rate = 2e-3;
    const float tau_pre = 20e-3;
    const float w_max = 1.;
    const float alpha = 3.; //time constant for moving average
    
    AdaptiveEBCPConnection * pyr1_to_pyr2_exc = new AdaptiveEBCPConnection(pyr1, pyr2, w_pyr1_to_pyr2_exc, p_pyr1_to_pyr2, learning_rate, w_max, tau_pre, GLUT);
    pyr1_to_pyr2_exc->set_target("g_ampa");
    pyr1_to_pyr2_exc->set_post_trace_event_tau(alpha);
    pyr1_to_pyr2_exc->set_post_trace_burst_tau(alpha);
    pyr1_to_pyr2_exc->max_rate = 15.;
    
    // Pyr1 to pyr2 - Inhibition
    TransmitEventConnection * pyr1_to_pyr2_inh = new TransmitEventConnection(pyr1, pyr2, w_pyr1_to_pyr2_inh, p_pyr1_to_pyr2, GABA);
    pyr1_to_pyr2_inh->set_target("g_gaba");
    
    
    //------ CONNECT FeedBACK ------//
    // Pyr2 to pyr1 - Excitation
    float p_pyr2_to_pyr1 = 0.05;
    TransmitBurstConnection * pyr2_to_pyr1_exc = new TransmitBurstConnection(pyr2, pyr1, w_pyr2_to_pyr1_exc, p_pyr2_to_pyr1, GLUT);
    pyr2_to_pyr1_exc->set_target("g_ampa_dend");
    
    // Pyr2 to pyr1 - Inhibition
    TransmitBurstConnection * pyr2_to_pyr1_inh = new TransmitBurstConnection(pyr2, pyr1, w_pyr2_to_pyr1_inh, p_pyr2_to_pyr1, GABA);
    pyr2_to_pyr1_inh->set_target("g_gaba_dend");
    
    
    /**************************************************/
    /******         CURRENT INJECTORS         *********/
    /**************************************************/
    // File-based current injections
    FileCurrentInjector * curr_inject_soma1 = new FileCurrentInjector(pyr1,"../data/credit-assign/current_soma.txt", "mem");
    FileCurrentInjector * curr_inject_dend2 = new FileCurrentInjector(pyr2,"../data/credit-assign/current_dend.txt", "Vd");
    
    // Standard current injections
    CurrentInjector * curr_inject_dend1 = new CurrentInjector(pyr1, "Vd");
    CurrentInjector * const_curr_inject_dend2 = new CurrentInjector(pyr2, "Vd");
    CurrentInjector * curr_inject_soma2 = new CurrentInjector(pyr2, "mem");
    
    
    /**************************************************/
    /******              MONITORS             *********/
    /**************************************************/
    double binSize_rate = 0.5; // s
    
    // Burst and event rate monitors
    auto seed_str = std::to_string(seed);
    BurstRateMonitor * brmon1 = new BurstRateMonitor( pyr1, sys->fn("brate1_seed"+seed_str), binSize_rate);
    BurstRateMonitor * brmon2 = new BurstRateMonitor( pyr2, sys->fn("brate2_seed"+seed_str), binSize_rate);
    
    // Voltage monitors
    VoltageMonitor *pyr1_mem = new VoltageMonitor(pyr1, 0,  sys->fn("mem1"), 1.e-3);
    VoltageMonitor *pyr2_mem = new VoltageMonitor(pyr2, 0,  sys->fn("mem2"), 1.e-3);
    StateMonitor *pyr1_Vd    = new StateMonitor(pyr1, 0, "Vd", sys->fn("Vd"), 1.e-3);
    StateMonitor *pyr2_Vd    = new StateMonitor(pyr2, 0, "Vd", sys->fn("Vd2"), 1.e-3);
    
    // Monitors for  estimating burst probability
    std::vector< StateMonitor* > smon_tr_post_ev;
    std::vector< StateMonitor* > smon_tr_post_b;
    
    for (int i=0;i<50;i++){
        StateMonitor * ev = new StateMonitor( pyr1_to_pyr2_exc->get_tr_event(), i, sys->fn(simname,i,"trevent"), binSize_rate);
        smon_tr_post_ev.push_back(ev);
        StateMonitor * b = new StateMonitor( pyr1_to_pyr2_exc->get_tr_burst(), i, sys->fn(simname,i,"trburst"), binSize_rate);
        smon_tr_post_b.push_back(b);
    }
    
    // Weight monitor
    WeightSumMonitor * wsmon = new WeightSumMonitor( pyr1_to_pyr2_exc, sys->fn("wsum"), binSize_rate );

    /**************************************************/
    /******             SIMULATION            *********/
    /**************************************************/
    curr_inject_soma2->set_all_currents(-200.e-12/pyr2[0].get_Cs());
    curr_inject_dend1->set_all_currents(00.e-12/pyr1[0].get_Cd());
    const_curr_inject_dend2->set_all_currents(325.e-12/pyr2[0].get_Cd());

    logger->msg("Running ...",PROGRESS);
    
    // relaxation period
    pyr1_to_pyr2_exc->stdp_active = false;
    sys->run(10.);
    
    const double durex = 15.;
    const double learning_duration = 0.5*durex;
    const double prediction_duration = 0.5*durex;
    
    // run first example
    sys->run(prediction_duration);
    pyr1_to_pyr2_exc->stdp_active = true;
    sys->run(1*learning_duration);
    pyr1_to_pyr2_exc->stdp_active = false;
    sys->run(0*learning_duration);

    // run second example
    sys->run(prediction_duration);
    pyr1_to_pyr2_exc->stdp_active = true;
    sys->run(1*learning_duration);
    pyr1_to_pyr2_exc->stdp_active = false;
    sys->run(0.*learning_duration);
    
    // run third example
    sys->run(prediction_duration);
    pyr1_to_pyr2_exc->stdp_active = true;
    sys->run(1*learning_duration);
    pyr1_to_pyr2_exc->stdp_active = false;
    sys->run(0.*learning_duration);
    
    if (errcode)
        auryn_abort(errcode);
    
    logger->msg("Freeing ...",PROGRESS,true);
    auryn_free();
    return errcode;
}


void initialize_pyr_neurons(NaudGroup* pyr) {
    pyr->syn_exc_soma->set_nmda_ampa_current_ampl_ratio(0.);
    pyr->syn_exc_dend->set_nmda_ampa_current_ampl_ratio(0.);
    pyr->random_mem(-70e-3, 5.e-3);
}

