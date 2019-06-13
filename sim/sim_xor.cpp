//  Training of a XOR gate using spiking burstprop.
//
//
//  Created by Alexandre Payeur on 2/19/19.
//

#include "auryn.h"
#include <cmath>
#include <fstream>
#include <string>
#include "BurstPoissonGroup.h"
#include "TransmitBurstConnection.h"
#include "TransmitEventConnection.h"

using namespace auryn;

namespace po = boost::program_options;

void fix_parameters_pv_neurons(AdExGroup* pv);
void initialize_pyr_neurons(NaudGroup* pyr);
void fix_parameters_som_neurons(AdExGroup* som);
void set_Facilitating_connection(STPConnectionETM * pre_to_post);
void set_Depressing_connection(STPConnectionETM * pre_to_post);

int main(int ac, char* av[])
{
    int errcode = 0;
    char strbuf [255];
    unsigned int seed = 1;
    string dir = "./";
    string simname = "xor";
    
    //************************************/
    //***          Options            ***//
    //************************************/
    
    try {
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help", "produce help message")
        ("dir", po::value<string>(), "output directory")
        ("seed", po::value<int>(), "random seed");
        
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
    
    
    
    //**********************************************/
    //***          Neural populations            ***/
    //**********************************************/
    NeuronID number_of_neurons = 4000;
    NeuronID number_inh_neurons = number_of_neurons/4;

    //-- PYRAMIDAL NEURONS
    // Output pyramidal population
    BurstPoissonGroup* output_pyr = new BurstPoissonGroup(number_of_neurons);
    initialize_pyr_neurons(output_pyr);

    // Hidden layer
    BurstPoissonGroup* hidden_pyr1 = new BurstPoissonGroup(number_of_neurons);
    initialize_pyr_neurons(hidden_pyr1);

    BurstPoissonGroup* hidden_pyr2 = new BurstPoissonGroup(number_of_neurons);
    initialize_pyr_neurons(hidden_pyr2);
    
    // Input layer
    BurstPoissonGroup* input_pyr1 = new BurstPoissonGroup(number_of_neurons);
    initialize_pyr_neurons(input_pyr1);
    
    BurstPoissonGroup* input_pyr2 = new BurstPoissonGroup(number_of_neurons);
    initialize_pyr_neurons(input_pyr2);
    
    //-- PARVALBUMIN-POSITIVE INTERNEURONS
    // Output layer
    //AdExGroup* output_pv = new AdExGroup(number_inh_neurons);
    //fix_parameters_pv_neurons(output_pv);
    
    // Hidden layer
    //AdExGroup* hidden_pv1 = new AdExGroup(number_inh_neurons);
    //fix_parameters_pv_neurons(hidden_pv1);
    
    //AdExGroup* hidden_pv2 = new AdExGroup(number_inh_neurons);
    //fix_parameters_pv_neurons(hidden_pv2);
    
    //-- SOMATOSTATIN-POSITIVE INTERNEURONS
    // Hidden layer
    //AdExGroup* hidden_som1 = new AdExGroup(number_inh_neurons);
    //fix_parameters_som_neurons(hidden_som1);
    
    //AdExGroup* hidden_som2 = new AdExGroup(number_inh_neurons);
    //fix_parameters_som_neurons(hidden_som2);
    
    //-- POISSON BACKGROUND NOISE
    float poisson_rate = 5.;
    //- Output layer
    // Pyr
    PoissonGroup* exc_Poisson_output_pyr       = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_output_pyr       = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* exc_Poisson_dend_output_pyr  = new PoissonGroup(number_of_neurons, 50*poisson_rate);
    PoissonGroup* inh_Poisson_dend_output_pyr  = new PoissonGroup(number_of_neurons, 50*poisson_rate);
    
    // PV
    //PoissonGroup* exc_Poisson_output_pv        = new PoissonGroup(number_inh_neurons, 100*poisson_rate);
    //PoissonGroup* inh_Poisson_output_pv        = new PoissonGroup(number_inh_neurons, 100*poisson_rate);
    
    //- Hidden layer
    // Pyr1
    PoissonGroup* exc_Poisson_hidden_pyr1      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_hidden_pyr1      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* exc_Poisson_dend_hidden_pyr1 = new PoissonGroup(number_of_neurons, 50*poisson_rate);
    PoissonGroup* inh_Poisson_dend_hidden_pyr1 = new PoissonGroup(number_of_neurons, 50*poisson_rate);

    // Pyr2
    PoissonGroup* exc_Poisson_hidden_pyr2      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_hidden_pyr2      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* exc_Poisson_dend_hidden_pyr2 = new PoissonGroup(number_of_neurons, 50*poisson_rate);
    PoissonGroup* inh_Poisson_dend_hidden_pyr2 = new PoissonGroup(number_of_neurons, 50*poisson_rate);

    // PV1
    //PoissonGroup* exc_Poisson_hidden_pv1       = new PoissonGroup(number_inh_neurons, 100*poisson_rate);
    //PoissonGroup* inh_Poisson_hidden_pv1       = new PoissonGroup(number_inh_neurons, 100*poisson_rate);
    
    // PV2
    //PoissonGroup* exc_Poisson_hidden_pv2       = new PoissonGroup(number_inh_neurons, 100*poisson_rate);
    //PoissonGroup* inh_Poisson_hidden_pv2       = new PoissonGroup(number_inh_neurons, 100*poisson_rate);
    
    // SOM1
    //PoissonGroup* exc_Poisson_hidden_som1       = new PoissonGroup(number_inh_neurons, 100*poisson_rate);
    //PoissonGroup* inh_Poisson_hidden_som1       = new PoissonGroup(number_inh_neurons, 100*poisson_rate);
    
    // SOM2
    //PoissonGroup* exc_Poisson_hidden_som2       = new PoissonGroup(number_inh_neurons, 100*poisson_rate);
    //PoissonGroup* inh_Poisson_hidden_som2       = new PoissonGroup(number_inh_neurons, 100*poisson_rate);
    
    //- Input layer
    // Pyr1
    PoissonGroup* exc_Poisson_input_pyr1      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_input_pyr1      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* exc_Poisson_dend_input_pyr1 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_dend_input_pyr1 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    
    // Pyr2
    PoissonGroup* exc_Poisson_input_pyr2      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_input_pyr2      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* exc_Poisson_dend_input_pyr2 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_dend_input_pyr2 = new PoissonGroup(number_of_neurons, 100*poisson_rate);

    
    //**********************************************************************/
    //***                         Connectivity                           ***/
    //**********************************************************************/
    
    //-- CONNECT BACKGROUND POISSON
    
    // (1) Connect Poisson noise to all pyramidal neuron populations
        // a - Input pyramidal neurons
    float w_epois_to_soma = 0.35; // somata
    float ratio_ie_soma = 1.;
    float w_epois_to_dend = 0.05; // dendrites
    float ratio_ie_dend = 1.;
    
            // population 1
    IdentityConnection * epois_to_soma_input_pyr1 = new IdentityConnection(exc_Poisson_input_pyr1, input_pyr1, w_epois_to_soma, GLUT);
    epois_to_soma_input_pyr1->set_target("g_ampa");
    IdentityConnection * ipois_to_soma_input_pyr1 = new IdentityConnection(inh_Poisson_input_pyr1, input_pyr1, ratio_ie_soma*w_epois_to_soma, GABA);
    ipois_to_soma_input_pyr1->set_target("g_gaba");
    
    IdentityConnection * epois_to_dend_input_pyr1 = new IdentityConnection(exc_Poisson_dend_input_pyr1, input_pyr1, w_epois_to_dend, GLUT);
    epois_to_dend_input_pyr1->set_target("g_ampa_dend");
    IdentityConnection * ipois_to_dend_input_pyr1 = new IdentityConnection(inh_Poisson_dend_input_pyr1, input_pyr1, ratio_ie_dend*w_epois_to_dend, GABA);
    ipois_to_dend_input_pyr1->set_target("g_gaba_dend");
    
            // population 2
    IdentityConnection * epois_to_soma_input_pyr2 = new IdentityConnection(exc_Poisson_input_pyr2, input_pyr2, w_epois_to_soma, GLUT);
    epois_to_soma_input_pyr2->set_target("g_ampa");
    IdentityConnection * ipois_to_soma_input_pyr2 = new IdentityConnection(inh_Poisson_input_pyr2, input_pyr2, ratio_ie_soma*w_epois_to_soma, GABA);
    ipois_to_soma_input_pyr2->set_target("g_gaba");
    
    IdentityConnection * epois_to_dend_input_pyr2 = new IdentityConnection(exc_Poisson_dend_input_pyr2, input_pyr2, w_epois_to_dend, GLUT);
    epois_to_dend_input_pyr2->set_target("g_ampa_dend");
    IdentityConnection * ipois_to_dend_input_pyr2 = new IdentityConnection(inh_Poisson_dend_input_pyr2, input_pyr2, ratio_ie_dend*w_epois_to_dend, GABA);
    ipois_to_dend_input_pyr2->set_target("g_gaba_dend");

        // b - Hidden layer pyramidal neurons
            // population 1
    IdentityConnection * epois_to_soma_hidden_pyr1 = new IdentityConnection(exc_Poisson_hidden_pyr1, hidden_pyr1, w_epois_to_soma, GLUT);
    epois_to_soma_hidden_pyr1->set_target("g_ampa");
    IdentityConnection * ipois_to_soma_hidden_pyr1 = new IdentityConnection(inh_Poisson_hidden_pyr1, hidden_pyr1, ratio_ie_soma*w_epois_to_soma, GABA);
    ipois_to_soma_hidden_pyr1->set_target("g_gaba");
    
    IdentityConnection * epois_to_dend_hidden_pyr1 = new IdentityConnection(exc_Poisson_dend_hidden_pyr1, hidden_pyr1, w_epois_to_dend, GLUT);
    epois_to_dend_hidden_pyr1->set_target("g_ampa_dend");
    IdentityConnection * ipois_to_dend_hidden_pyr1 = new IdentityConnection(inh_Poisson_dend_hidden_pyr1, hidden_pyr1, ratio_ie_dend*w_epois_to_dend, GABA);
    ipois_to_dend_hidden_pyr1->set_target("g_gaba_dend");
    
            // population 2
    IdentityConnection * epois_to_soma_hidden_pyr2 = new IdentityConnection(exc_Poisson_hidden_pyr2, hidden_pyr2, w_epois_to_soma, GLUT);
    epois_to_soma_hidden_pyr2->set_target("g_ampa");
    IdentityConnection * ipois_to_soma_hidden_pyr2 = new IdentityConnection(inh_Poisson_hidden_pyr2, hidden_pyr2, ratio_ie_soma*w_epois_to_soma, GABA);
    ipois_to_soma_hidden_pyr2->set_target("g_gaba");
    
    IdentityConnection * epois_to_dend_hidden_pyr2 = new IdentityConnection(exc_Poisson_dend_hidden_pyr2, hidden_pyr2, w_epois_to_dend, GLUT);
    epois_to_dend_hidden_pyr2->set_target("g_ampa_dend");
    IdentityConnection * ipois_to_dend_hidden_pyr2 = new IdentityConnection(inh_Poisson_dend_hidden_pyr2, hidden_pyr2, ratio_ie_dend*w_epois_to_dend, GABA);
    ipois_to_dend_hidden_pyr2->set_target("g_gaba_dend");

        // c - Output layer pyramidal neurons
    IdentityConnection * epois_to_soma_output_pyr = new IdentityConnection(exc_Poisson_output_pyr, output_pyr, w_epois_to_soma, GLUT);
    epois_to_soma_output_pyr->set_target("g_ampa");
    IdentityConnection * ipois_to_soma_output_pyr = new IdentityConnection(inh_Poisson_output_pyr, output_pyr, ratio_ie_soma*w_epois_to_soma, GABA);
    ipois_to_soma_output_pyr->set_target("g_gaba");
    
    IdentityConnection * epois_to_dend_output_pyr = new IdentityConnection(exc_Poisson_dend_output_pyr, output_pyr, w_epois_to_dend, GLUT);
    epois_to_dend_output_pyr->set_target("g_ampa_dend");
    IdentityConnection * ipois_to_dend_output_pyr = new IdentityConnection(inh_Poisson_dend_output_pyr, output_pyr, ratio_ie_dend*w_epois_to_dend, GABA);
    ipois_to_dend_output_pyr->set_target("g_gaba_dend");
    

    //-- CONNECT FeedFORWARD
        // (1) Connect input layer to hidden layer
            // a - input 1 to hidden 1
    
    EBCPConnection *  = new EBCPConnection(input_pyr1, hidden_pyr1, we_soma, 0.5, 1e-3); //eta
    con_ext_soma->set_target("g_ampa");
    con_ext_soma->set_max_weight(1.0);
    
    
    
    //-- CONNECT FeedBACK
    
    
    
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

void set_Facilitating_connection(STPConnectionETM * pre_to_post){
    pre_to_post->set_tau_d(100e-3);
    pre_to_post->set_tau_f(100e-3);
    pre_to_post->set_ujump(0.02);
    pre_to_post->set_f(0.1);
}

void set_Depressing_connection(STPConnectionETM * pre_to_post){
    pre_to_post->set_tau_d(20.e-3);
    pre_to_post->set_tau_f(1.);
    pre_to_post->set_ujump(0.9);
    pre_to_post->set_f(0.1);
}

void connect_noise_to_pyramids(NaudGroup* pyr, )
