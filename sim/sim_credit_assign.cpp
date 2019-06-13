/***********************************************************
 Description: Example of "credit assignment" and error
 generation.
 ************************************************************/

#include "auryn.h"
#include <cmath>
#include <fstream>
#include <string>
#include "STPeTMConnection.h"
#include "FileCurrentInjector.h"

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
    string simname = "credit_assign";
    const NeuronID number_of_neurons = 8000;
    const NeuronID N_other_neuron = number_of_neurons/4;
    
    float w_pyr1_to_pyr2 = 0.05*4000/number_of_neurons; //0.012
    float w_pv_to_pyr2   = 0.05*1000/N_other_neuron; //0.
    float w_pyr2_to_pyr1 = 0.05;//0.05;
    float w_som_to_pyr1  = 0.0;//0.025;
    
    /**************************************************/
    /********              OPTIONS            *********/
    /**************************************************/
    try {
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help", "produce help message")
        ("dir", po::value<string>(), "output directory")
        ("seed", po::value<int>(), "random seed")
        ("w21", po::value<float>(), "weight from pyr1 to pyr2")
        ("w2pv", po::value<float>(), "weight from pv to pyr2")
        ("w12", po::value<float>(), "weight from pyr2 to pyr1")
        ("w1som", po::value<float>(), "weight from som to pyr1")
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
        
        if (vm.count("w21")) {
            std::cout << "weight from pyr1 to pyr2 set to "
            << vm["w21"].as<float>() << ".\n";
            w_pyr1_to_pyr2 = vm["w21"].as<float>();
        }
        
        if (vm.count("w2pv")) {
            std::cout << "weight from pv to pyr2 set to "
            << vm["w2pv"].as<float>() << ".\n";
            w_pv_to_pyr2 = vm["w2pv"].as<float>();
        }
        
        if (vm.count("w12")) {
            std::cout << "weight from pyr2 to pyr1 set to "
            << vm["w12"].as<float>() << ".\n";
            w_pyr2_to_pyr1 = vm["w12"].as<float>()*4000/number_of_neurons;
        }
        
        if (vm.count("w1som")) {
            std::cout << "weight from som to pyr1 set to "
            << vm["w1som"].as<float>() << ".\n";
            w_som_to_pyr1 = vm["w1som"].as<float>()*1000/N_other_neuron;
        }
    }
    catch(std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        std::cerr << "Exception of unknown type!\n";
    }
    
    std::cout<<"weight = "<<w_pyr2_to_pyr1<<std::endl;
    // INITIALIZE AURYN
    auryn_init( ac, av, dir, simname );
    sys->set_master_seed(seed);
    
    /**************************************************/
    /******          NEURON POPULATIONS       *********/
    /**************************************************/
    
    // First layer of 2-comp pyramidal neurons
    NaudGroup* pyr1 = new NaudGroup(number_of_neurons);
    initialize_pyr_neurons(pyr1);
    
    // Second layer of 2-comp pyramidal neurons
    NaudGroup* pyr2 = new NaudGroup(number_of_neurons);
    initialize_pyr_neurons(pyr2);
    
    // PV neurons
    AdExGroup* pv = new AdExGroup(N_other_neuron);
    fix_parameters_pv_neurons(pv);
    
    // SOM neurons
    AdExGroup* som = new AdExGroup(N_other_neuron);
    fix_parameters_som_neurons(som);
    
    // External populations of Poisson neurons (noise)
    // Pyr 1
    float poisson_rate = 5.;
    PoissonGroup* exc_Poisson1      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson1      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* exc_Poisson_dend1 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_dend1 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    
    // Pyr2
    PoissonGroup* exc_Poisson2      = new PoissonGroup(number_of_neurons, 75*poisson_rate);
    PoissonGroup* inh_Poisson2      = new PoissonGroup(number_of_neurons, 75*poisson_rate);
    PoissonGroup* exc_Poisson_dend2 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_dend2 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    
    // SOM
    PoissonGroup* exc_Poisson_som   = new PoissonGroup(N_other_neuron, 100*poisson_rate);
    PoissonGroup* inh_Poisson_som   = new PoissonGroup(N_other_neuron, 100*poisson_rate);
    
    // PV
    PoissonGroup* exc_Poisson_pv    = new PoissonGroup(N_other_neuron, 100*poisson_rate);
    PoissonGroup* inh_Poisson_pv    = new PoissonGroup(N_other_neuron, 100*poisson_rate);
    
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
    
    // External Poisson neurons -> SOM
    float w_epois_to_som = 0.35;
    float ratio_ie_som = 2.;
    IdentityConnection * epois_to_som = new IdentityConnection(exc_Poisson_som, som, w_epois_to_som, GLUT);
    IdentityConnection * ipois_to_som = new IdentityConnection(inh_Poisson_som, som, w_epois_to_som*ratio_ie_som, GABA);
    
    // External Poisson neurons -> PV
    float w_epois_to_pv = 0.75;
    float ratio_ie_pv = 1.5;
    IdentityConnection * epois_to_pv = new IdentityConnection(exc_Poisson_pv, pv, w_epois_to_pv, GLUT);
    IdentityConnection * ipois_to_pv = new IdentityConnection(inh_Poisson_pv, pv, w_epois_to_pv*ratio_ie_pv, GABA);
    
    
    //------ CONNECT FeedFORWARD ------//
    
    // Pyr1 to pyr2 - STD
    float p_pyr1_to_pyr2 = 0.05; //0.05
    STPeTMConnection * pyr1_to_pyr2 = new STPeTMConnection(pyr1, pyr2, w_pyr1_to_pyr2, p_pyr1_to_pyr2, GLUT);
    set_Depressing_connection(pyr1_to_pyr2);
    pyr1_to_pyr2->set_target("g_ampa");
    
    // Pyr1 to PV - STD
    float w_pyr1_to_pv = 0.05*4000/number_of_neurons; //0.04
    float p_pyr1_to_pv = 0.05;
    STPeTMConnection * pyr1_to_pv = new STPeTMConnection(pyr1, pv, w_pyr1_to_pv, p_pyr1_to_pv, GLUT);
    set_Depressing_connection(pyr1_to_pv);
    
    // PV to Pyr2
    float p_pv_to_pyr2 = 0.05;
    SparseConnection * pv_to_pyr2 = new SparseConnection(pv, pyr2, w_pv_to_pyr2, p_pv_to_pyr2, GABA);
    pv_to_pyr2->set_target("g_gaba");
    
    
    //------ CONNECT FeedBACK ------//
    // Pyr2 to pyr1 - STF
    float p_pyr2_to_pyr1 = 0.05;
    STPeTMConnection * pyr2_to_pyr1 = new STPeTMConnection(pyr2, pyr1, w_pyr2_to_pyr1, p_pyr2_to_pyr1, GLUT);
    set_Facilitating_connection(pyr2_to_pyr1);
    pyr2_to_pyr1->set_target("g_ampa_dend");
    
    // Pyr2 to SOM - STF
    float w_pyr2_to_som = 0.4*4000/number_of_neurons; //0.01
    float p_pyr2_to_som = 0.05;
    STPeTMConnection * pyr2_to_som = new STPeTMConnection(pyr2, som, w_pyr2_to_som, p_pyr2_to_som, GLUT);
    set_Facilitating_connection(pyr2_to_som);
    
    // SOM to pyr1
    float p_som_to_pyr1 = 0.05; //0.05
    SparseConnection * som_to_pyr1 = new SparseConnection(som, pyr1, w_som_to_pyr1, p_som_to_pyr1, GABA);
    som_to_pyr1->set_target("g_gaba_dend");
    

    // -- OTHER CONNECTIONS
    // Pyr to Pyr - fake inh STF
    float w_pyr_to_pyr = 0.1*4000/number_of_neurons; //0.1
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
    // File-based current injections
    FileCurrentInjector * curr_inject_soma1 = new FileCurrentInjector(pyr1,"../data/credit-assign/current_soma.txt", "mem");
    FileCurrentInjector * curr_inject_dend2 = new FileCurrentInjector(pyr2,"../data/credit-assign/current_dend.txt", "Vd");
    
    // Standard current injections
    CurrentInjector * curr_inject_som   = new CurrentInjector(som, "mem");
    CurrentInjector * curr_inject_pv    = new CurrentInjector(pv, "mem");
    CurrentInjector * curr_inject_dend1 = new CurrentInjector(pyr1, "Vd");
    CurrentInjector * curr_inject_soma2 = new CurrentInjector(pyr2, "mem");

    
    /**************************************************/
    /******              MONITORS             *********/
    /**************************************************/
    double binSize_rate = 20.e-3; // ms
    
    // Burst and event rate monitors
    auto seed_str = std::to_string(seed);
    BurstRateMonitor * brmon1 = new BurstRateMonitor( pyr1, sys->fn("brate1_seed"+seed_str), binSize_rate);
    BurstRateMonitor * brmon2 = new BurstRateMonitor( pyr2, sys->fn("brate2_seed"+seed_str), binSize_rate);
    
    // Population monitors
    PopulationRateMonitor * pv_activity  = new PopulationRateMonitor( pv,  sys->fn("pvrate_seed"+seed_str), binSize_rate );
    PopulationRateMonitor * som_activity = new PopulationRateMonitor( som, sys->fn("somrate_seed"+seed_str), binSize_rate );
    
    
    /**************************************************/
    /******             SIMULATION            *********/
    /**************************************************/
    curr_inject_soma2->set_all_currents(-100.e-12/pyr2[0].get_Cs());
    curr_inject_dend1->set_all_currents(-50.e-12/pyr1[0].get_Cd());
    curr_inject_som->set_all_currents(0e-12/som[0].get_c_mem());
    curr_inject_pv->set_all_currents(205e-12/pv[0].get_c_mem());
    
    logger->msg("Running ...",PROGRESS);
    sys->run(1.+3*1.5);
    
    
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



