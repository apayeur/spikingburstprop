/***********************************************************
 Description: Compute the feedforward transfer function.
************************************************************/

#include "auryn.h"
#include <cmath>
#include <fstream>
#include <string>
#include "STPeTMConnection.h"
#include "BurstPoissonGroup.h"
#include "TransmitBurstConnection.h"
#include "TransmitEventConnection.h"

using namespace auryn;

namespace po = boost::program_options;

void initialize_pyr_neurons(NaudGroup* pyr);


int main(int ac, char* av[])
{
    int errcode = 0;
    char strbuf [255];
    unsigned int seed = 1;
    string dir = "../data/fftransfer/for-xor";
    string simname = "fftransfer";
    
    const NeuronID number_of_neurons = 2000;
    const NeuronID N_other_neuron = number_of_neurons/4;
    
    float w_pyr1_to_pyr2_exc = 0.07*4000/number_of_neurons; //0.07
    float w_pyr1_to_pyr2_inh = 0.03*4000/number_of_neurons; //0.03
    
    float w_pyr2_to_pyr1_exc = 0.12*4000/number_of_neurons;
    float w_pyr2_to_pyr1_inh = 0.03*4000/number_of_neurons;
    
    float w_pyr2_to_pyr3_exc = 0.09*4000/number_of_neurons; //0.07
    float w_pyr2_to_pyr3_inh = 0.05*4000/number_of_neurons; //0.03
    
    float w_pyr3_to_pyr2_exc = 0.12*4000/number_of_neurons;
    float w_pyr3_to_pyr2_inh = 0.03*4000/number_of_neurons;


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
    sys->set_master_seed(seed);
    
    /**************************************************/
    /******          NEURON POPULATIONS       *********/
    /**************************************************/
    
    // First layer of 2-comp pyramidal neurons
    BurstPoissonGroup* pyr1 = new BurstPoissonGroup(number_of_neurons);
    initialize_pyr_neurons(pyr1);
    
    BurstPoissonGroup* pyr12 = new BurstPoissonGroup(number_of_neurons);
    initialize_pyr_neurons(pyr12);
    
    
    // Second layer of 2-comp pyramidal neurons
    BurstPoissonGroup* pyr2 = new BurstPoissonGroup(number_of_neurons);
    initialize_pyr_neurons(pyr2);
    
    BurstPoissonGroup* pyr22 = new BurstPoissonGroup(number_of_neurons);
    initialize_pyr_neurons(pyr22);
    
    
    // Third layer of 2-comp pyramidal neurons
    BurstPoissonGroup* pyr3 = new BurstPoissonGroup(number_of_neurons);
    initialize_pyr_neurons(pyr3);
    
    
    // External populations of Poisson neurons (noise)
        // Pyr 1
    float poisson_rate = 5.;
    PoissonGroup* exc_Poisson1      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson1      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* exc_Poisson_dend1 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_dend1 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    
        // Pyr2
    PoissonGroup* exc_Poisson2      = new PoissonGroup(number_of_neurons, 20*poisson_rate);
    PoissonGroup* inh_Poisson2      = new PoissonGroup(number_of_neurons, 50*poisson_rate);
    PoissonGroup* exc_Poisson_dend2 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_dend2 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    
        // Pyr3
    PoissonGroup* exc_Poisson3      = new PoissonGroup(number_of_neurons, 50*poisson_rate);
    PoissonGroup* inh_Poisson3      = new PoissonGroup(number_of_neurons, 50*poisson_rate);
    PoissonGroup* exc_Poisson_dend3 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_dend3 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    
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
    
    IdentityConnection * epois_to_soma12 = new IdentityConnection(exc_Poisson1, pyr12, w_epois_to_soma, GLUT);
    epois_to_soma12->set_target("g_ampa");
    IdentityConnection * ipois_to_soma12 = new IdentityConnection(inh_Poisson1, pyr12, ratio_ie_soma*w_epois_to_soma, GABA);
    ipois_to_soma12->set_target("g_gaba");
    
    float w_epois_to_dend = 0.05;
    float ratio_ie_dend = 1.;
    IdentityConnection * epois_to_dend1 = new IdentityConnection(exc_Poisson_dend1, pyr1, w_epois_to_dend, GLUT);
    epois_to_dend1->set_target("g_ampa_dend");
    IdentityConnection * ipois_to_dend1 = new IdentityConnection(inh_Poisson_dend1, pyr1, ratio_ie_dend*w_epois_to_dend, GABA);
    ipois_to_dend1->set_target("g_gaba_dend");
    
    IdentityConnection * epois_to_dend12 = new IdentityConnection(exc_Poisson_dend1, pyr12, w_epois_to_dend, GLUT);
    epois_to_dend12->set_target("g_ampa_dend");
    IdentityConnection * ipois_to_dend12 = new IdentityConnection(inh_Poisson_dend1, pyr12, ratio_ie_dend*w_epois_to_dend, GABA);
    ipois_to_dend12->set_target("g_gaba_dend");
    
    // External Poisson neurons -> Pyr2
    IdentityConnection * epois_to_soma2 = new IdentityConnection(exc_Poisson2, pyr2, w_epois_to_soma, GLUT);
    epois_to_soma2->set_target("g_ampa");
    IdentityConnection * ipois_to_soma2 = new IdentityConnection(inh_Poisson2, pyr2, ratio_ie_soma*w_epois_to_soma, GABA);
    ipois_to_soma2->set_target("g_gaba");
    
    IdentityConnection * epois_to_soma22 = new IdentityConnection(exc_Poisson2, pyr22, w_epois_to_soma, GLUT);
    epois_to_soma22->set_target("g_ampa");
    IdentityConnection * ipois_to_soma22 = new IdentityConnection(inh_Poisson2, pyr22, ratio_ie_soma*w_epois_to_soma, GABA);
    ipois_to_soma22->set_target("g_gaba");
    
    IdentityConnection * epois_to_dend2 = new IdentityConnection(exc_Poisson_dend2, pyr2, w_epois_to_dend, GLUT);
    epois_to_dend2->set_target("g_ampa_dend");
    IdentityConnection * ipois_to_dend2 = new IdentityConnection(inh_Poisson_dend2, pyr2, ratio_ie_dend*w_epois_to_dend, GABA);
    ipois_to_dend2->set_target("g_gaba_dend");
    
    IdentityConnection * epois_to_dend22 = new IdentityConnection(exc_Poisson_dend2, pyr22, w_epois_to_dend, GLUT);
    epois_to_dend22->set_target("g_ampa_dend");
    IdentityConnection * ipois_to_dend22 = new IdentityConnection(inh_Poisson_dend2, pyr22, ratio_ie_dend*w_epois_to_dend, GABA);
    ipois_to_dend22->set_target("g_gaba_dend");
    
    
    // External Poisson neurons -> Pyr3
    IdentityConnection * epois_to_soma3 = new IdentityConnection(exc_Poisson3, pyr3, w_epois_to_soma, GLUT);
    epois_to_soma3->set_target("g_ampa");
    IdentityConnection * ipois_to_soma3 = new IdentityConnection(inh_Poisson3, pyr3, ratio_ie_soma*w_epois_to_soma, GABA);
    ipois_to_soma3->set_target("g_gaba");
    IdentityConnection * epois_to_dend3 = new IdentityConnection(exc_Poisson_dend3, pyr3, w_epois_to_dend, GLUT);
    epois_to_dend3->set_target("g_ampa_dend");
    IdentityConnection * ipois_to_dend3 = new IdentityConnection(inh_Poisson_dend3, pyr3, ratio_ie_dend*w_epois_to_dend, GABA);
    ipois_to_dend3->set_target("g_gaba_dend");
    
    
    //------ CONNECT FeedFORWARD ------//
    // Pyr1 to pyr2 - Excitation
    float p_pyr1_to_pyr2 = 0.05; //0.05
    TransmitEventConnection * pyr1_to_pyr2_exc = new TransmitEventConnection(pyr1, pyr2, w_pyr1_to_pyr2_exc, p_pyr1_to_pyr2, GLUT);
    pyr1_to_pyr2_exc->set_target("g_ampa");
    
    // Pyr1 to pyr2 - Inhibition
    TransmitEventConnection * pyr1_to_pyr2_inh = new TransmitEventConnection(pyr1, pyr2, w_pyr1_to_pyr2_inh, p_pyr1_to_pyr2, GABA);
    pyr1_to_pyr2_inh->set_target("g_gaba");
    
    // Pyr1 to pyr22 - Excitation
    TransmitEventConnection * pyr1_to_pyr22_exc = new TransmitEventConnection(pyr1, pyr22, w_pyr1_to_pyr2_exc, p_pyr1_to_pyr2, GLUT);
    pyr1_to_pyr22_exc->set_target("g_ampa");
    
    // Pyr1 to pyr22 - Inhibition
    TransmitEventConnection * pyr1_to_pyr22_inh = new TransmitEventConnection(pyr1, pyr22, w_pyr1_to_pyr2_inh, p_pyr1_to_pyr2, GABA);
    pyr1_to_pyr22_inh->set_target("g_gaba");
    
    // Pyr12 to pyr2 - Excitation
    TransmitEventConnection * pyr12_to_pyr2_exc = new TransmitEventConnection(pyr12, pyr2, w_pyr1_to_pyr2_exc, p_pyr1_to_pyr2, GLUT);
    pyr12_to_pyr2_exc->set_target("g_ampa");
    
    // Pyr12 to pyr2 - Inhibition
    TransmitEventConnection * pyr12_to_pyr2_inh = new TransmitEventConnection(pyr12, pyr2, w_pyr1_to_pyr2_inh, p_pyr1_to_pyr2, GABA);
    pyr12_to_pyr2_inh->set_target("g_gaba");
    
    // Pyr12 to pyr22 - Excitation
    TransmitEventConnection * pyr12_to_pyr22_exc = new TransmitEventConnection(pyr12, pyr22, w_pyr1_to_pyr2_exc, p_pyr1_to_pyr2, GLUT);
    pyr12_to_pyr22_exc->set_target("g_ampa");
    
    // Pyr12 to pyr22 - Inhibition
    TransmitEventConnection * pyr12_to_pyr22_inh = new TransmitEventConnection(pyr12, pyr22, w_pyr1_to_pyr2_inh, p_pyr1_to_pyr2, GABA);
    pyr12_to_pyr22_inh->set_target("g_gaba");
    
    
    
    
    // Pyr2 to pyr3 - Excitation
    float p_pyr2_to_pyr3 = 0.05; //0.05
    TransmitEventConnection * pyr2_to_pyr3_exc = new TransmitEventConnection(pyr2, pyr3, w_pyr2_to_pyr3_exc, p_pyr2_to_pyr3, GLUT);
    pyr2_to_pyr3_exc->set_target("g_ampa");
    
    // Pyr2 to pyr3 - Inhibition
    TransmitEventConnection * pyr2_to_pyr3_inh = new TransmitEventConnection(pyr2, pyr3, w_pyr2_to_pyr3_inh, p_pyr2_to_pyr3, GABA);
    pyr2_to_pyr3_inh->set_target("g_gaba");
    
    // Pyr22 to pyr3 - Excitation
    TransmitEventConnection * pyr22_to_pyr3_exc = new TransmitEventConnection(pyr22, pyr3, w_pyr2_to_pyr3_exc, p_pyr2_to_pyr3, GLUT);
    pyr22_to_pyr3_exc->set_target("g_ampa");
    
    // Pyr22 to pyr3 - Inhibition
    TransmitEventConnection * pyr22_to_pyr3_inh = new TransmitEventConnection(pyr22, pyr3, w_pyr2_to_pyr3_inh, p_pyr2_to_pyr3, GABA);
    pyr22_to_pyr3_inh->set_target("g_gaba");
    
    
    //------ CONNECT FeedBACK ------//
    // Pyr3 to pyr2 - Excitation
    float p_pyr3_to_pyr2 = 0.05;
    TransmitBurstConnection * pyr3_to_pyr2_exc = new TransmitBurstConnection(pyr3, pyr2, w_pyr3_to_pyr2_exc, p_pyr3_to_pyr2, GLUT);
    pyr3_to_pyr2_exc->set_target("g_ampa_dend");
    
    // Pyr3 to pyr2 - Inhibition
    TransmitBurstConnection * pyr3_to_pyr2_inh = new TransmitBurstConnection(pyr3, pyr2, w_pyr3_to_pyr2_inh, p_pyr3_to_pyr2, GABA);
    pyr3_to_pyr2_inh->set_target("g_gaba_dend");
    
    // Pyr3 to pyr22 - Excitation
    TransmitBurstConnection * pyr3_to_pyr22_exc = new TransmitBurstConnection(pyr3, pyr22, w_pyr3_to_pyr2_exc, p_pyr3_to_pyr2, GLUT);
    pyr3_to_pyr22_exc->set_target("g_ampa_dend");
    
    // Pyr3 to pyr22 - Inhibition
    TransmitBurstConnection * pyr3_to_pyr22_inh = new TransmitBurstConnection(pyr3, pyr22, w_pyr3_to_pyr2_inh, p_pyr3_to_pyr2, GABA);
    pyr3_to_pyr22_inh->set_target("g_gaba_dend");


    /**************************************************/
    /******         CURRENT INJECTORS         *********/
    /**************************************************/
    CurrentInjector * curr_inject_soma1 = new CurrentInjector(pyr1, "mem");
    CurrentInjector * curr_inject_soma12 = new CurrentInjector(pyr12, "mem");
    CurrentInjector * curr_inject_soma2 = new CurrentInjector(pyr2, "mem");
    CurrentInjector * curr_inject_soma22 = new CurrentInjector(pyr22, "mem");
    CurrentInjector * curr_inject_soma3 = new CurrentInjector(pyr3, "mem");

    
    /**************************************************/
    /******              MONITORS             *********/
    /**************************************************/
    double binSize_rate = 20.e-3; // ms
    auto seed_str = std::to_string(seed);
    
    // Burst/event rate monitors
    BurstRateMonitor * brmon1 = new BurstRateMonitor( pyr1, sys->fn("brate1"), binSize_rate);
    BurstRateMonitor * brmon2 = new BurstRateMonitor( pyr2, sys->fn("brate2"), binSize_rate);
    BurstRateMonitor * brmon3 = new BurstRateMonitor( pyr3, sys->fn("brate3"), binSize_rate);
    
    // Voltage monitors
    VoltageMonitor *pyr1_mem = new VoltageMonitor(pyr1, 0,  sys->fn("pyr1_mem"), 1.e-3);
    VoltageMonitor *pyr2_mem = new VoltageMonitor(pyr2, 0,  sys->fn("pyr2_mem"), 1.e-3);
    VoltageMonitor *pyr3_mem = new VoltageMonitor(pyr3, 0,  sys->fn("pyr3_mem"), 1.e-3);
    
    //Raster
    SpikeMonitor * smon1 = new SpikeMonitor( pyr1, sys->fn("ras1"));
    SpikeMonitor * smon2 = new SpikeMonitor( pyr2, sys->fn("ras2"));

    /**************************************************/
    /******             SIMULATION            *********/
    /**************************************************/
    logger->msg("Running ...",PROGRESS);
    
    const double max_somatic_current  = 300e-12;
    const double min_somatic_current  = -200e-12;
    const double somatic_current_incr = 100e-12;
    const double simtime = 1000e-3;
    
    double somatic_current = min_somatic_current;
    curr_inject_soma1->set_all_currents(somatic_current/pyr1[0].get_Cs());
    curr_inject_soma12->set_all_currents(min_somatic_current/pyr12[0].get_Cs());
    
    curr_inject_soma2->set_all_currents(-400e-12/pyr2[0].get_Cs());
    curr_inject_soma22->set_all_currents(-400e-12/pyr22[0].get_Cs());
   
    curr_inject_soma3->set_all_currents(-200e-12/pyr3[0].get_Cs());

    while (somatic_current < max_somatic_current + somatic_current_incr/2)
    {
        curr_inject_soma1->set_all_currents(somatic_current/pyr1[0].get_Cs());
        somatic_current += somatic_current_incr;
        sys->run(simtime);
    }
    
    somatic_current = min_somatic_current;
    curr_inject_soma1->set_all_currents(somatic_current/pyr1[0].get_Cs());
    curr_inject_soma12->set_all_currents(max_somatic_current/pyr12[0].get_Cs());
    
    while (somatic_current < max_somatic_current + somatic_current_incr/2)
    {
        curr_inject_soma1->set_all_currents(somatic_current/pyr1[0].get_Cs());
        somatic_current += somatic_current_incr;
        sys->run(simtime);
    }
    

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



