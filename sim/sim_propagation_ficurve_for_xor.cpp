/***********************************************************
 Description: Compute ficurves of pre- et postsynaptic
 populations in the context of the FF and FB propagation.
 Helps fine-tune parameters.
************************************************************/

#include "auryn.h"
#include <cmath>
#include <fstream>
#include <string>
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
    string dir = "../data/propagation/ficurve/for-xor";
    string simname = "propagation_ficurve";
    double somatic_current = 100e-12;
    
    const NeuronID number_of_neurons = 4000;
    const NeuronID N_other_neuron = number_of_neurons/4;

    float w_pyr1_to_pyr2_exc = 0.07*4000/number_of_neurons;
    float w_pyr1_to_pyr2_inh = 0.03*4000/number_of_neurons; //0.03
    
    float w_pyr2_to_pyr1_exc = 0.12*4000/number_of_neurons;
    float w_pyr2_to_pyr1_inh = 0.03*4000/number_of_neurons;


    /**************************************************/
    /********              OPTIONS            *********/
    /**************************************************/
    try {
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help", "produce help message")
        ("dir", po::value<string>(), "output directory")
        ("seed", po::value<int>(), "random seed")
        ("somacurrent", po::value<double>(), "somatic current")
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
        
        if (vm.count("somacurrent")) {
            std::cout << "soma current set to "
            << vm["somacurrent"].as<double>() << ".\n";
            somatic_current = vm["somacurrent"].as<double>();
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
    TransmitEventConnection * pyr1_to_pyr2_exc = new TransmitEventConnection(pyr1, pyr2, w_pyr1_to_pyr2_exc, p_pyr1_to_pyr2, GLUT);
    pyr1_to_pyr2_exc->set_target("g_ampa");
    
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
    CurrentInjector * curr_inject_soma1 = new CurrentInjector(pyr1, "mem");
    CurrentInjector * curr_inject_dend1 = new CurrentInjector(pyr1, "Vd");
    CurrentInjector * curr_inject_soma2 = new CurrentInjector(pyr2, "mem");
    CurrentInjector * curr_inject_dend2 = new CurrentInjector(pyr2, "Vd");
    
    /**************************************************/
    /******              MONITORS             *********/
    /**************************************************/
    double binSize_rate = 20.e-3; // ms
    auto label_str = std::to_string(somatic_current*1.e12);
    auto seed_str = std::to_string(seed);
    
    // Burst/event rate monitors
    BurstRateMonitor * brmon1 = new BurstRateMonitor( pyr1, sys->fn("brate1_somacurrent"+label_str), binSize_rate);
    BurstRateMonitor * brmon2 = new BurstRateMonitor( pyr2, sys->fn("brate2_somacurrent"+label_str), binSize_rate);
    SpikeMonitor * smon1 = new SpikeMonitor( pyr1, sys->fn("ras1_somacurrent"+label_str) );
    SpikeMonitor * smon2 = new SpikeMonitor( pyr2, sys->fn("ras2_somacurrent"+label_str) );

    /**************************************************/
    /******             SIMULATION            *********/
    /**************************************************/
    logger->msg("Running ...",PROGRESS);
    
    const double max_dendritic_current  = 1600e-12;
    const double min_dendritic_current  = -400e-12;
    const double dendritic_current_incr = 200e-12;
    const double simtime = 1000e-3;
    
    double dendritic_current = min_dendritic_current;
    curr_inject_soma1->set_all_currents(somatic_current/pyr1[0].get_Cs());
    curr_inject_dend1->set_all_currents(0e-12/pyr1[0].get_Cd());

    while (dendritic_current < max_dendritic_current + dendritic_current_incr/2)
    {
        curr_inject_dend2->set_all_currents(dendritic_current/pyr2[0].get_Cd());
        dendritic_current += dendritic_current_incr;
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

