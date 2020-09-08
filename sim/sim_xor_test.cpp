//  Training of a XOR gate using spiking burstprop.
//
//
//  Created by Alexandre Payeur on 2/19/19.
//

#include "auryn.h"
#include <cmath>
#include <fstream>
#include <string>
#include <cstdlib>
#include <time.h>
#include <iostream>
#include <algorithm>

#include "BurstPoissonGroup.h"
#include "TransmitBurstConnection.h"
#include "TransmitEventConnection.h"
#include "AdaptiveEBCPConnection.h"
#include "BiasConnection.h"
#include "BiasIdentityConnection.h"

using namespace auryn;

namespace po = boost::program_options;

// function declarations
void initialize_pyr_neurons(NaudGroup* pyr);
void generate_error(AdaptiveEBCPConnection * connect,  CurrentInjector * inject_dend, int output, float capacitance);
double compute_cost(AdaptiveEBCPConnection * connect, int output);
void fix_plastic_connections_parameters(AdaptiveEBCPConnection * connect, double alpha);
void generate_random_permutation(std::vector<int> &v);
void fix_bias_connections_parameters(BiasIdentityConnection * connect, double alpha);

// definition of global constants
NeuronID number_of_neurons = 2000;
const NeuronID number_inh_neurons = number_of_neurons/4;
const AurynFloat e_max = 10.; // maximum event rate
const AurynFloat e_min = 2.; // minimum event rate (for error gen and cost function)

// definition of global variables
int count_example = 0;

// structure for XOR gate input/output pairs
struct input_output_XOR {
    short input1;
    short input2;
    double input_current1;
    double input_current2;
    short output;
    double current_for_one;
    double current_for_zero;
    
    input_output_XOR(double curr_one, double curr_zero):current_for_one(curr_one), current_for_zero(curr_zero) {}
    input_output_XOR():current_for_one(200.e-12), current_for_zero(-200e-12) {}
    
    void select_XOR_example_random() {
        input1 = rand() % 2;
        input2 = rand() % 2;
        
        std::cout<<std::endl<<"input = ("<<input1<<", "<<input2<<")"<<std::endl;

        if ( input1==0 and input2==0) output = 0;
        else if ( input1==1 and input2==0) output = 1;
        else if ( input1==0 and input2==1) output = 1;
        else output = 0;
        
        input_current1 = current_for_zero + (current_for_one-current_for_zero)*input1;
        input_current2 = current_for_zero + (current_for_one-current_for_zero)*input2;

    }
    
    void select_XOR_example_deterministic() {
        int selected_input = count_example % 4;
        apply(selected_input);
    }
    
    void select_XOR_example_deterministic(int i) {
        int selected_input = i % 4;
        apply(selected_input);
    }
    
    void select_XOR_example_specific(short a1, short a2) {
        input1 = a1;
        input2 = a2;
        
        std::cout<<std::endl<<"input = ("<<input1<<", "<<input2<<")"<<std::endl;
        
        if ( input1==0 and input2==0) output = 0;
        else if ( input1==1 and input2==0) output = 1;
        else if ( input1==0 and input2==1) output = 1;
        else output = 0;
        
        input_current1 = current_for_zero + (current_for_one-current_for_zero)*input1;
        input_current2 = current_for_zero + (current_for_one-current_for_zero)*input2;
        
    }
    
    void apply(int selected_input){
        switch (selected_input) {
            case 0:
                input1=0;
                input2=0;
                output=0;
                break;
            case 1:
                input1=1;
                input2=0;
                output=1;
                break;
            case 2:
                input1=0;
                input2=1;
                output=1;
                break;
            case 3:
                input1=1;
                input2=1;
                output=0;
        }
        
        std::cout<<std::endl<<"input = ("<<input1<<", "<<input2<<")"<<std::endl;
        
        input_current1 = current_for_zero + (current_for_one-current_for_zero)*input1;
        input_current2 = current_for_zero + (current_for_one-current_for_zero)*input2;
    }
};



int main(int ac, char* av[])
{
    int errcode = 0;
    char strbuf [255];
    unsigned int seed = 1;
    double durex = 10.;
    string dir = "./";
    string simname = "xor-test";
    double moving_average_time_constant = 4.;
    int number_of_training_examples = 1;
    
    const double learning_rate_hid_to_out = 2.*2.e-3;
    const double learning_rate_in_to_hid  = 2.*4.e-3;
    const double learning_rate_bias       = 0.1e-3;
    const double max_weight    = 1.0;
    const double tau_pre       = 100.e-3;
    std::ofstream outfile_cost;
    std::ofstream outfile_cost_epoch;
    srand (1);


    //****************************************************//
    //***                    OPTIONS                   ***//
    //****************************************************//
    
    try {
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help", "produce help message")
        ("dir", po::value<string>(), "output directory")
        ("seed", po::value<int>(), "random seed")
        ("numex", po::value<int>(), "number of training examples")
        ("N", po::value<int>(), "number of neurons per population")
        ("durex", po::value<double>(), "duration of an example")
        ("alpha", po::value<double>(), "moving average time constant")
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
        
        if (vm.count("numex")) {
            std::cout << "number of examples set to "
            << vm["numex"].as<int>() << ".\n";
            number_of_training_examples = vm["numex"].as<int>();
        }

        if (vm.count("N")) {
            std::cout << "number of neurons set to "
            << vm["N"].as<int>() << ".\n";
            number_of_neurons = vm["N"].as<int>();
        }
        
        if (vm.count("durex")) {
            std::cout << "example duration set to "
            << vm["durex"].as<double>() << ".\n";
            durex = vm["durex"].as<double>();
        }
        
        if (vm.count("alpha")) {
            std::cout << "moving average time constant set to "
            << vm["alpha"].as<double>() << ".\n";
            moving_average_time_constant = vm["alpha"].as<double>();
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
    sys->set_master_seed(seed);
    outfile_cost.open(dir+"/cost2.dat");
    outfile_cost_epoch.open(dir+"/cost_epoch.dat");

    //*********************************************************//
    //***                NEURAL POPULATIONS                 ***//
    //*********************************************************//

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
    
    //-- POISSON BACKGROUND NOISE
    float poisson_rate = 5.;
    //- Output layer
    // Pyr
    PoissonGroup* exc_Poisson_output_pyr       = new PoissonGroup(number_of_neurons, 30*poisson_rate);
    PoissonGroup* inh_Poisson_output_pyr       = new PoissonGroup(number_of_neurons, 30*poisson_rate);
    PoissonGroup* exc_Poisson_dend_output_pyr  = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_dend_output_pyr  = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    
    //- Hidden layer
    // Pyr1
    PoissonGroup* exc_Poisson_hidden_pyr1      = new PoissonGroup(number_of_neurons, 30*poisson_rate);
    PoissonGroup* inh_Poisson_hidden_pyr1      = new PoissonGroup(number_of_neurons, 30*poisson_rate);
    PoissonGroup* exc_Poisson_dend_hidden_pyr1 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_dend_hidden_pyr1 = new PoissonGroup(number_of_neurons, 100*poisson_rate);

    // Pyr2
    PoissonGroup* exc_Poisson_hidden_pyr2      = new PoissonGroup(number_of_neurons, 30*poisson_rate);
    PoissonGroup* inh_Poisson_hidden_pyr2      = new PoissonGroup(number_of_neurons, 30*poisson_rate);
    PoissonGroup* exc_Poisson_dend_hidden_pyr2 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_dend_hidden_pyr2 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    
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
    //***                         CONNECTIVITY                           ***/
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
    double w_bias = 0.1*w_epois_to_soma;
    BiasIdentityConnection * epois_to_soma_hidden_pyr1 = new BiasIdentityConnection(exc_Poisson_hidden_pyr1, hidden_pyr1, 1.1*w_epois_to_soma, learning_rate_bias, GLUT);
    fix_bias_connections_parameters(epois_to_soma_hidden_pyr1, moving_average_time_constant);

    IdentityConnection * ipois_to_soma_hidden_pyr1 = new IdentityConnection(inh_Poisson_hidden_pyr1, hidden_pyr1, 1.4*ratio_ie_soma*w_epois_to_soma, GABA);
    ipois_to_soma_hidden_pyr1->set_target("g_gaba");
    
    IdentityConnection * epois_to_dend_hidden_pyr1 = new IdentityConnection(exc_Poisson_dend_hidden_pyr1, hidden_pyr1, w_epois_to_dend, GLUT);
    epois_to_dend_hidden_pyr1->set_target("g_ampa_dend");
    IdentityConnection * ipois_to_dend_hidden_pyr1 = new IdentityConnection(inh_Poisson_dend_hidden_pyr1, hidden_pyr1, ratio_ie_dend*w_epois_to_dend, GABA);
    ipois_to_dend_hidden_pyr1->set_target("g_gaba_dend");
    
            // population 2
    BiasIdentityConnection * epois_to_soma_hidden_pyr2 = new BiasIdentityConnection(exc_Poisson_hidden_pyr2, hidden_pyr2, 1.1*w_epois_to_soma, learning_rate_bias, GLUT);
    fix_bias_connections_parameters(epois_to_soma_hidden_pyr2, moving_average_time_constant);

    IdentityConnection * ipois_to_soma_hidden_pyr2 = new IdentityConnection(inh_Poisson_hidden_pyr2, hidden_pyr2, 1.4*ratio_ie_soma*w_epois_to_soma, GABA);
    ipois_to_soma_hidden_pyr2->set_target("g_gaba");
    
    IdentityConnection * epois_to_dend_hidden_pyr2 = new IdentityConnection(exc_Poisson_dend_hidden_pyr2, hidden_pyr2, w_epois_to_dend, GLUT);
    epois_to_dend_hidden_pyr2->set_target("g_ampa_dend");
    IdentityConnection * ipois_to_dend_hidden_pyr2 = new IdentityConnection(inh_Poisson_dend_hidden_pyr2, hidden_pyr2, ratio_ie_dend*w_epois_to_dend, GABA);
    ipois_to_dend_hidden_pyr2->set_target("g_gaba_dend");

        // c - Output layer pyramidal neurons
    BiasIdentityConnection * epois_to_soma_output_pyr = new BiasIdentityConnection(exc_Poisson_output_pyr, output_pyr, w_epois_to_soma, learning_rate_bias, GLUT);
    fix_bias_connections_parameters(epois_to_soma_output_pyr, moving_average_time_constant);

    IdentityConnection * ipois_to_soma_output_pyr = new IdentityConnection(inh_Poisson_output_pyr, output_pyr, ratio_ie_soma*w_epois_to_soma, GABA);
    ipois_to_soma_output_pyr->set_target("g_gaba");
    
    IdentityConnection * epois_to_dend_output_pyr = new IdentityConnection(exc_Poisson_dend_output_pyr, output_pyr, w_epois_to_dend, GLUT);
    epois_to_dend_output_pyr->set_target("g_ampa_dend");
    IdentityConnection * ipois_to_dend_output_pyr = new IdentityConnection(inh_Poisson_dend_output_pyr, output_pyr, ratio_ie_dend*w_epois_to_dend, GABA);
    ipois_to_dend_output_pyr->set_target("g_gaba_dend");
    

    //-- CONNECT FeedFORWARD
    const float p_ff = 0.05*2000/number_of_neurons;
    float w_in_to_hid_exc  = 0.18;
    float w_in_to_hid_inh  = 0.1;
    float w_hid_to_out_exc = 0.18;
    float w_hid_to_out_inh = 0.1;
    
        // (1) Connect input layer to hidden layer
            // Excitation
    AdaptiveEBCPConnection * con_input1_to_hidden1_exc = new AdaptiveEBCPConnection(input_pyr1, hidden_pyr1, 1.3*w_in_to_hid_exc, p_ff, learning_rate_in_to_hid, max_weight, tau_pre);
    fix_plastic_connections_parameters(con_input1_to_hidden1_exc, moving_average_time_constant);
    
    AdaptiveEBCPConnection * con_input2_to_hidden1_exc = new AdaptiveEBCPConnection(input_pyr2, hidden_pyr1, 1.3*w_in_to_hid_exc, p_ff, learning_rate_in_to_hid, max_weight, tau_pre);
    fix_plastic_connections_parameters(con_input2_to_hidden1_exc, moving_average_time_constant);
    
    AdaptiveEBCPConnection * con_input1_to_hidden2_exc = new AdaptiveEBCPConnection(input_pyr1, hidden_pyr2, 1.3*w_in_to_hid_exc, p_ff, learning_rate_in_to_hid, max_weight, tau_pre);
    fix_plastic_connections_parameters(con_input1_to_hidden2_exc, moving_average_time_constant);
    
    AdaptiveEBCPConnection * con_input2_to_hidden2_exc = new AdaptiveEBCPConnection(input_pyr2, hidden_pyr2, 1.3*w_in_to_hid_exc, p_ff, learning_rate_in_to_hid, max_weight, tau_pre);
    fix_plastic_connections_parameters(con_input2_to_hidden2_exc, moving_average_time_constant);

            // Inhibition
    TransmitEventConnection * con_input1_to_hidden1_inh = new TransmitEventConnection(input_pyr1, hidden_pyr1, 2.5*w_in_to_hid_inh, p_ff, GABA);
    con_input1_to_hidden1_inh->set_target("g_gaba");
    
    TransmitEventConnection * con_input2_to_hidden1_inh = new TransmitEventConnection(input_pyr2, hidden_pyr1, 2.5*w_in_to_hid_inh, p_ff, GABA);
    con_input2_to_hidden1_inh->set_target("g_gaba");
    
    TransmitEventConnection * con_input1_to_hidden2_inh = new TransmitEventConnection(input_pyr1, hidden_pyr2, 2.5*w_in_to_hid_inh, p_ff, GABA);
    con_input1_to_hidden2_inh->set_target("g_gaba");
    
    TransmitEventConnection * con_input2_to_hidden2_inh = new TransmitEventConnection(input_pyr2, hidden_pyr2, 2.5*w_in_to_hid_inh, p_ff, GABA);
    con_input2_to_hidden2_inh->set_target("g_gaba");
    
    
    // (2) Connect hidden layer to output layer
    // Excitation
    AdaptiveEBCPConnection * con_hidden1_to_output_exc = new AdaptiveEBCPConnection(hidden_pyr1, output_pyr, 1.3*w_hid_to_out_exc, p_ff, learning_rate_hid_to_out, max_weight, tau_pre);
    fix_plastic_connections_parameters(con_hidden1_to_output_exc, moving_average_time_constant);
    
    AdaptiveEBCPConnection * con_hidden2_to_output_exc = new AdaptiveEBCPConnection(hidden_pyr2, output_pyr, 1.3*w_hid_to_out_exc, p_ff, learning_rate_hid_to_out, max_weight, tau_pre);
    fix_plastic_connections_parameters(con_hidden2_to_output_exc, moving_average_time_constant);

    // Inhibition
    TransmitEventConnection * con_hidden1_to_output_inh = new TransmitEventConnection(hidden_pyr1, output_pyr, 2.5*w_hid_to_out_inh, p_ff, GABA);
    con_hidden1_to_output_inh->set_target("g_gaba");
    
    TransmitEventConnection * con_hidden2_to_output_inh = new TransmitEventConnection(hidden_pyr2, output_pyr, 2.5*w_hid_to_out_inh, p_ff, GABA);
    con_hidden2_to_output_inh->set_target("g_gaba");
    

    //-- CONNECT FeedBACK
    float w_fb_exc = 0.26;
    float w_fb_inh = 0.06;
    float p_fb = p_ff;
    
    // Excitation
    TransmitBurstConnection * con_output_to_hidden1_exc = new TransmitBurstConnection(output_pyr, hidden_pyr1, w_fb_exc, p_fb, GLUT);
    con_output_to_hidden1_exc->set_target("g_ampa_dend");
    TransmitBurstConnection * con_output_to_hidden2_exc = new TransmitBurstConnection(output_pyr, hidden_pyr2, 0.1*w_fb_exc, p_fb, GLUT);
    con_output_to_hidden2_exc->set_target("g_ampa_dend");
    
    // Inhibition
    TransmitBurstConnection * con_output_to_hidden1_inh = new TransmitBurstConnection(output_pyr, hidden_pyr1, w_fb_inh, p_fb, GABA);
    con_output_to_hidden1_inh->set_target("g_gaba_dend");
    TransmitBurstConnection * con_output_to_hidden2_inh = new TransmitBurstConnection(output_pyr, hidden_pyr2, 4.5*w_fb_inh, p_fb, GABA);
    con_output_to_hidden2_inh->set_target("g_gaba_dend");
    
    
    /**************************************************/
    /******         CURRENT INJECTORS         *********/
    /**************************************************/
    CurrentInjector * baseline_curr_output_pyr_dend   = new CurrentInjector(output_pyr, "Vd");
    CurrentInjector * baseline_curr_hidden_pyr2_dend  = new CurrentInjector(hidden_pyr2, "Vd");
    CurrentInjector * curr_inject_output_pyr_dend     = new CurrentInjector(output_pyr, "Vd");
    CurrentInjector * curr_inject_hidden_pyr1_dend    = new CurrentInjector(hidden_pyr1, "Vd");
    CurrentInjector * curr_inject_hidden_pyr2_dend    = new CurrentInjector(hidden_pyr2, "Vd");

    CurrentInjector * curr_inject_output_pyr_soma  = new CurrentInjector(output_pyr, "mem");
    CurrentInjector * curr_inject_hidden_pyr1_soma = new CurrentInjector(hidden_pyr1, "mem");
    CurrentInjector * curr_inject_hidden_pyr2_soma = new CurrentInjector(hidden_pyr2, "mem");
    CurrentInjector * curr_inject_input_pyr1_soma  = new CurrentInjector(input_pyr1, "mem");
    CurrentInjector * curr_inject_input_pyr2_soma  = new CurrentInjector(input_pyr2, "mem");
    
    // constant currents injected in dendrites (mostly to shift dendritic transfer function)
    baseline_curr_output_pyr_dend->set_all_currents(400.e-12/output_pyr[0].get_Cd());
    baseline_curr_hidden_pyr2_dend->set_all_currents(500.e-12/hidden_pyr2[0].get_Cd());
    
    // constant currents injected in somata (adjust gain; equivalent to changing leak reversal)
    curr_inject_hidden_pyr1_soma->set_all_currents(-00.e-12/hidden_pyr1[0].get_Cs());
    curr_inject_hidden_pyr2_soma->set_all_currents(-00.e-12/hidden_pyr2[0].get_Cs());
    curr_inject_output_pyr_soma->set_all_currents(0.e-12/output_pyr[0].get_Cs());
    
    // error-encoding current
    curr_inject_output_pyr_dend->set_all_currents(0.);
    
    
    /**************************************************/
    /******              MONITORS             *********/
    /**************************************************/
    double binSize_rate = 1;
    double binSize_wsum = binSize_rate;
    auto seed_str = std::to_string(seed);

    // Online monitor
    //sys->set_online_rate_monitor_target(hidden_pyr1);
    //sys->set_online_rate_monitor_tau(binSize_rate);

    // Burst/event rate monitors
    BurstRateMonitor * brmon_input1  = new BurstRateMonitor( input_pyr1, sys->fn("brate_input1_seed"+seed_str), binSize_rate);
    BurstRateMonitor * brmon_input2  = new BurstRateMonitor( input_pyr2, sys->fn("brate_input2_seed"+seed_str), binSize_rate);
    BurstRateMonitor * brmon_hidden1 = new BurstRateMonitor( hidden_pyr1, sys->fn("brate_hidden1_seed"+seed_str), binSize_rate);
    BurstRateMonitor * brmon_hidden2 = new BurstRateMonitor( hidden_pyr2, sys->fn("brate_hidden2_seed"+seed_str), binSize_rate);
    BurstRateMonitor * brmon_output  = new BurstRateMonitor( output_pyr, sys->fn("brate_output_seed"+seed_str), binSize_rate);

    // Weight sum monitors
    WeightSumMonitor * wsmon_in1_to_hid1 = new WeightSumMonitor( con_input1_to_hidden1_exc, sys->fn("wsum_in1_to_hid1_seed"+seed_str), binSize_wsum );
    WeightSumMonitor * wsmon_in1_to_hid2 = new WeightSumMonitor( con_input1_to_hidden2_exc, sys->fn("wsum_in1_to_hid2_seed"+seed_str), binSize_wsum );
    WeightSumMonitor * wsmon_in2_to_hid1 = new WeightSumMonitor( con_input2_to_hidden1_exc, sys->fn("wsum_in2_to_hid1_seed"+seed_str), binSize_wsum );
    WeightSumMonitor * wsmon_in2_to_hid2 = new WeightSumMonitor( con_input2_to_hidden2_exc, sys->fn("wsum_in2_to_hid2_seed"+seed_str), binSize_wsum );
    WeightSumMonitor * wsmon_hid1_to_out = new WeightSumMonitor( con_hidden1_to_output_exc, sys->fn("wsum_hid1_to_out_seed"+seed_str), binSize_wsum );
    WeightSumMonitor * wsmon_hid2_to_out = new WeightSumMonitor( con_hidden2_to_output_exc, sys->fn("wsum_hid2_to_out_seed"+seed_str), binSize_wsum );
    
    WeightSumMonitor * wsmon_noise_to_out = new WeightSumMonitor( epois_to_soma_output_pyr, sys->fn("wsum_noise_to_out_seed"+seed_str), binSize_wsum );
    WeightSumMonitor * wsmon_noise_to_hid1 = new WeightSumMonitor( epois_to_soma_hidden_pyr1, sys->fn("wsum_noise_to_hid1_seed"+seed_str), binSize_wsum );
    WeightSumMonitor * wsmon_noise_to_hid2 = new WeightSumMonitor( epois_to_soma_hidden_pyr2, sys->fn("wsum_noise_to_hid2_seed"+seed_str), binSize_wsum );

    // Monitors for  estimating burst probability (a sample of 50 neurons)
    std::vector< StateMonitor* > smon_tr_post_ev;
    std::vector< StateMonitor* > smon_tr_post_b;
    
    for (int i=0;i<50;i++){
        StateMonitor * ev = new StateMonitor( con_hidden1_to_output_exc->get_tr_event(), i, sys->fn(simname,i,"trevent_seed"+seed_str), binSize_rate);
        smon_tr_post_ev.push_back(ev);
        StateMonitor * b = new StateMonitor( con_hidden1_to_output_exc->get_tr_burst(), i, sys->fn(simname,i,"trburst_seed"+seed_str), binSize_rate);
        smon_tr_post_b.push_back(b);
    }
     
    //**********************************************************************/
    //***                          SIMULATIONS                           ***/
    //**********************************************************************/
    const double max_input_current = 150e-12;
    const double min_input_current = -350e-12;
    
    input_output_XOR example(max_input_current, min_input_current);
    const int number_of_epochs = number_of_training_examples/4;
    
    double Cs_input_pyr1 = input_pyr1[0].get_Cs();
    double Cs_input_pyr2 = input_pyr2[0].get_Cs();
    double Cd_output_pyr = output_pyr[0].get_Cd();
    
    double cost = 0;
    double cost_epoch = 0;
    std::vector<int> permutation; // when the examples are presented in random order
    
    // burn-in period (plasticity is off)
    con_input1_to_hidden1_exc->stdp_active = false;
    con_input1_to_hidden2_exc->stdp_active = false;
    con_input2_to_hidden1_exc->stdp_active = false;
    con_input2_to_hidden2_exc->stdp_active = false;
    
    con_hidden1_to_output_exc->stdp_active = false;
    con_hidden2_to_output_exc->stdp_active = false;

    epois_to_soma_hidden_pyr1->stdp_active = false;
    epois_to_soma_hidden_pyr2->stdp_active = false;
    epois_to_soma_output_pyr->stdp_active =  false;
    
    const float prediction_time = 0.9*durex;
    const float learning_time = durex - prediction_time;
    const int number_of_error_recomputing = 1;

    sys->run(3*moving_average_time_constant);

    // test set before learning
    // input (0,0)
    curr_inject_input_pyr1_soma->set_all_currents(min_input_current/Cs_input_pyr1);
    curr_inject_input_pyr2_soma->set_all_currents(min_input_current/Cs_input_pyr2);
    sys->run(durex);
    
    // input (1,0)
    curr_inject_input_pyr1_soma->set_all_currents(max_input_current/Cs_input_pyr1);
    curr_inject_input_pyr2_soma->set_all_currents(min_input_current/Cs_input_pyr2);
    sys->run(durex);
    
    // input (0,1)
    curr_inject_input_pyr1_soma->set_all_currents(min_input_current/Cs_input_pyr1);
    curr_inject_input_pyr2_soma->set_all_currents(max_input_current/Cs_input_pyr2);
    sys->run(durex);
    
    // input (1,1)
    curr_inject_input_pyr1_soma->set_all_currents(max_input_current/Cs_input_pyr1);
    curr_inject_input_pyr2_soma->set_all_currents(max_input_current/Cs_input_pyr2);
    sys->run(durex);
    
    
    // loop over set of training examples
    count_example = 1;
    for (int epoch=1;epoch<=number_of_epochs;epoch++){
        // Uncomment the following line if the examples are randomly selected.
        //generate_random_permutation(permutation);
        cost_epoch = 0;
        for (int i=0;i<4;i++){
            // select example (random or deterministic)
            example.select_XOR_example_deterministic(i);
            // Uncomment the following line and comment the previous one if the examples are randomly selected.
            //example.apply(permutation[i]);
            
            // inject input currents according to the selected example
            curr_inject_input_pyr1_soma->set_all_currents(example.input_current1/Cs_input_pyr1);
            curr_inject_input_pyr2_soma->set_all_currents(example.input_current2/Cs_input_pyr2);
            
            // first, plasticity is off to build up the moving average and formulate prediction
            con_input1_to_hidden1_exc->stdp_active = false;
            con_input1_to_hidden2_exc->stdp_active = false;
            con_input2_to_hidden1_exc->stdp_active = false;
            con_input2_to_hidden2_exc->stdp_active = false;

            con_hidden1_to_output_exc->stdp_active = false;
            con_hidden2_to_output_exc->stdp_active = false;
            
            epois_to_soma_hidden_pyr1->stdp_active = false;
            epois_to_soma_hidden_pyr2->stdp_active = false;
            epois_to_soma_output_pyr->stdp_active  = false;
            
            // run for the prediction time
            sys->run(prediction_time);
            
            // compute cost at the end of the prediction interval
            cost = compute_cost(con_hidden1_to_output_exc, example.output);
            std::cout<<"example #"<< count_example <<": cost (pre) = "<<cost<<std::endl;
            outfile_cost<<count_example<<"\t"<<cost;
            
            cost_epoch += cost;
            
            // turn plasticity on
            con_input1_to_hidden1_exc->stdp_active = true;
            con_input1_to_hidden2_exc->stdp_active = true;
            con_input2_to_hidden1_exc->stdp_active = true;
            con_input2_to_hidden2_exc->stdp_active = true;
            
            con_hidden1_to_output_exc->stdp_active = true;
            con_hidden2_to_output_exc->stdp_active = true;
            
            // change "false" to "true" to turn ON the plasticity of the noise/bias (default : false)
            epois_to_soma_hidden_pyr1->stdp_active = false;
            epois_to_soma_hidden_pyr2->stdp_active = false;
            epois_to_soma_output_pyr->stdp_active  = false;

            for (int k=0;k<number_of_error_recomputing;k++){
                // inject current in output pyramids' dendrites to generate error
                generate_error(con_hidden1_to_output_exc, curr_inject_output_pyr_dend, example.output, Cd_output_pyr);
                
                // propagate credit and learn
                sys->run(learning_time/number_of_error_recomputing);
            }
            
            // compute cost after learning interval for reference
            cost = compute_cost(con_hidden1_to_output_exc, example.output);
            std::cout<<"example #"<< count_example <<": cost (post) = "<<cost<<std::endl;
            outfile_cost<<"\t"<<cost<<std::endl;

            // back to baseline current injection in output dendrites
            curr_inject_output_pyr_dend->set_all_currents(0.);
            
            count_example++;
        }
        std::cout << "**** Cost epoch " << epoch << " = " << cost_epoch << std::endl;
        outfile_cost_epoch<<epoch<<"\t"<<cost_epoch<<std::endl;
    }
    // test set post
    con_input1_to_hidden1_exc->stdp_active = false;
    con_input1_to_hidden2_exc->stdp_active = false;
    con_input2_to_hidden1_exc->stdp_active = false;
    con_input2_to_hidden2_exc->stdp_active = false;
    
    con_hidden1_to_output_exc->stdp_active = false;
    con_hidden2_to_output_exc->stdp_active = false;
    
    epois_to_soma_hidden_pyr1->stdp_active = false;
    epois_to_soma_hidden_pyr2->stdp_active = false;
    epois_to_soma_output_pyr->stdp_active =  false;
    
    // input (0,0)
    curr_inject_input_pyr1_soma->set_all_currents(min_input_current/Cs_input_pyr1);
    curr_inject_input_pyr2_soma->set_all_currents(min_input_current/Cs_input_pyr2);
        sys->run(durex);
    
    // input (1,0)
    curr_inject_input_pyr1_soma->set_all_currents(max_input_current/Cs_input_pyr1);
    curr_inject_input_pyr2_soma->set_all_currents(min_input_current/Cs_input_pyr2);
    sys->run(durex);
    
    // input (0,1)
    curr_inject_input_pyr1_soma->set_all_currents(min_input_current/Cs_input_pyr1);
    curr_inject_input_pyr2_soma->set_all_currents(max_input_current/Cs_input_pyr2);
    sys->run(durex);
    
    // input (1,1)
    curr_inject_input_pyr1_soma->set_all_currents(max_input_current/Cs_input_pyr1);
    curr_inject_input_pyr2_soma->set_all_currents(max_input_current/Cs_input_pyr2);
    sys->run(durex);
    
    if (errcode)
        auryn_abort(errcode);
    
    logger->msg("Freeing ...",PROGRESS,true);
    auryn_free();
    outfile_cost.close();
    return errcode;
}

void initialize_pyr_neurons(NaudGroup* pyr) {
    pyr->syn_exc_soma->set_nmda_ampa_current_ampl_ratio(0.);
    pyr->syn_exc_dend->set_nmda_ampa_current_ampl_ratio(0.);
    pyr->random_mem(-70e-3, 5.e-3);
}

void generate_error(AdaptiveEBCPConnection * connect,  CurrentInjector * inject_dend, int output, float capacitance){
    Trace * tmp_event_trace = connect->get_tr_event();
    AurynFloat injected_current;
    AurynFloat er_estimate;
    double base_deflection = 1000e-12;
    double desired_deltap = 0.;
    double tolerance = 1.; // how far from e_min (for output=0) and e_max (output=1) must the prediction be to be correct
    
    for (NeuronID i=0; i<number_of_neurons;i++){
        er_estimate = tmp_event_trace->normalized_get(i);
        if ( output == 0 ){
            desired_deltap += -1/(e_max - er_estimate);
            if (er_estimate < e_max - tolerance) {
                if (er_estimate <= e_min + tolerance) injected_current = 0.;
                else injected_current = -base_deflection/(e_max - er_estimate);
            }
            else injected_current = -base_deflection;
        }
        else {
            desired_deltap += 1/(er_estimate - e_min);
            if (er_estimate > e_min + tolerance) {
                if (er_estimate >= e_max - tolerance) injected_current = 0.;
                else injected_current = base_deflection/(er_estimate - e_min);
            }
            else injected_current = base_deflection;
        }
        inject_dend->set_current(i, injected_current/capacitance);
    }
    std::cout<<" desired delta p = "<<desired_deltap/number_of_neurons<<std::endl;

}

double compute_cost(AdaptiveEBCPConnection * connect, int output){
    // compute mean event rate over population
    double mean_er = 0.;
    Trace * tmp_event_trace = connect->get_tr_event();

    for (NeuronID i=0; i<number_of_neurons;i++){
        mean_er += tmp_event_trace->normalized_get(i);
    }
    
    mean_er /= number_of_neurons;
    std::cout<<" Event rate = "<<mean_er<<std::endl;

    mean_er = (mean_er - e_min)/(e_max - e_min);
    
    // compute cost
    if (mean_er > 0.99) {
        mean_er = 0.99;
    }
    else if (mean_er < 0.01) mean_er = 0.01;
    double cost = -output*log(mean_er) - (1 - output)*log(1. - mean_er);
    return cost;
    
}

void fix_plastic_connections_parameters(AdaptiveEBCPConnection * connect, double alpha){
    connect->set_post_trace_event_tau(alpha);
    connect->set_post_trace_burst_tau(alpha);
    connect->set_target("g_ampa");
    connect->max_rate = e_max;
    connect->min_rate = e_min;
}

void fix_bias_connections_parameters(BiasIdentityConnection * connect, double alpha){
    connect->set_target("g_ampa");
    connect->set_post_trace_event_tau(alpha);
    connect->set_post_trace_burst_tau(alpha);
    connect->max_rate = e_max;
    connect->min_rate = e_min;
    connect->stdp_active = false;
}


void generate_random_permutation(std::vector<int> &v){
    v.clear();
    for (int i=0;i<4;i++) v.push_back(i);
    
    random_shuffle ( v.begin(), v.end() );
}
