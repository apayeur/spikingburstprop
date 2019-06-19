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

#include "BurstPoissonGroup.h"
#include "TransmitBurstConnection.h"
#include "TransmitEventConnection.h"
#include "AdaptiveEBCPConnection.h"


using namespace auryn;

namespace po = boost::program_options;

// function declarations
void fix_parameters_pv_neurons(AdExGroup* pv);
void initialize_pyr_neurons(NaudGroup* pyr);
void fix_parameters_som_neurons(AdExGroup* som);
//void set_Facilitating_connection(STPConnectionETM * pre_to_post);
//void set_Depressing_connection(STPConnectionETM * pre_to_post);
void generate_error(AdaptiveEBCPConnection * connect,  CurrentInjector * inject_dend, int output, float capacitance);
double compute_cost(AdaptiveEBCPConnection * connect, int output);
void fix_plastic_connections_parameters(AdaptiveEBCPConnection * connect, double alpha);

// definitions of global constants
const NeuronID number_of_neurons = 4000;
const NeuronID number_inh_neurons = number_of_neurons/4;
const AurynFloat e_max = 10.; // maximum event rate

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
    
    void select_XOR_example() {
        input1 = rand() % 2;
        input2 = rand() % 2;
        
        std::cout<<"input = ("<<input1<<", "<<input2<<")"<<std::endl;

        if ( input1==0 and input2==0) output = 0;
        else if ( input1==1 and input2==0) output = 1;
        else if ( input1==0 and input2==1) output = 1;
        else output = 0;
        
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
    string simname = "xor";
    double moving_average_time_constant = 4.;
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
    PoissonGroup* exc_Poisson_output_pyr       = new PoissonGroup(number_of_neurons, 50*poisson_rate);
    PoissonGroup* inh_Poisson_output_pyr       = new PoissonGroup(number_of_neurons, 50*poisson_rate);
    PoissonGroup* exc_Poisson_dend_output_pyr  = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_dend_output_pyr  = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    
    //- Hidden layer
    // Pyr1
    PoissonGroup* exc_Poisson_hidden_pyr1      = new PoissonGroup(number_of_neurons, 50*poisson_rate);
    PoissonGroup* inh_Poisson_hidden_pyr1      = new PoissonGroup(number_of_neurons, 50*poisson_rate);
    PoissonGroup* exc_Poisson_dend_hidden_pyr1 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_dend_hidden_pyr1 = new PoissonGroup(number_of_neurons, 100*poisson_rate);

    // Pyr2
    PoissonGroup* exc_Poisson_hidden_pyr2      = new PoissonGroup(number_of_neurons, 50*poisson_rate);
    PoissonGroup* inh_Poisson_hidden_pyr2      = new PoissonGroup(number_of_neurons, 50*poisson_rate);
    PoissonGroup* exc_Poisson_dend_hidden_pyr2 = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_dend_hidden_pyr2 = new PoissonGroup(number_of_neurons, 100*poisson_rate);

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
    const float p_ff = 0.05; //0.05
    float w_ff_exc = 0.035*4000/number_of_neurons;
    float w_ff_inh = 0.015*4000/number_of_neurons;
    const double learning_rate = 1e-3;
    const double max_weight    = 1.0;
    const double tau_pre       = 20.e-3;
    
        // (1) Connect input layer to hidden layer
            // Excitation
    AdaptiveEBCPConnection * con_input1_to_hidden1_exc = new AdaptiveEBCPConnection(input_pyr1, hidden_pyr1, 1.3*w_ff_exc, p_ff, learning_rate, max_weight, tau_pre);
    fix_plastic_connections_parameters(con_input1_to_hidden1_exc, moving_average_time_constant);
    
    AdaptiveEBCPConnection * con_input1_to_hidden2_exc = new AdaptiveEBCPConnection(input_pyr1, hidden_pyr2, 0.8*w_ff_exc, p_ff, learning_rate, max_weight, tau_pre);
    fix_plastic_connections_parameters(con_input1_to_hidden2_exc, moving_average_time_constant);

    AdaptiveEBCPConnection * con_input2_to_hidden1_exc = new AdaptiveEBCPConnection(input_pyr2, hidden_pyr1, 1.2*w_ff_exc, p_ff, learning_rate, max_weight, tau_pre);
    fix_plastic_connections_parameters(con_input2_to_hidden1_exc, moving_average_time_constant);
    
    AdaptiveEBCPConnection * con_input2_to_hidden2_exc = new AdaptiveEBCPConnection(input_pyr2, hidden_pyr2, 0.7*w_ff_exc, p_ff, learning_rate, max_weight, tau_pre);
    fix_plastic_connections_parameters(con_input2_to_hidden2_exc, moving_average_time_constant);

            // Inhibition
    TransmitEventConnection * con_input1_to_hidden1_inh = new TransmitEventConnection(input_pyr1, hidden_pyr1, w_ff_inh, p_ff, GABA);
    con_input1_to_hidden1_inh->set_target("g_gaba");
    
    TransmitEventConnection * con_input1_to_hidden2_inh = new TransmitEventConnection(input_pyr1, hidden_pyr2, w_ff_inh, p_ff, GABA);
    con_input1_to_hidden2_inh->set_target("g_gaba");
    
    TransmitEventConnection * con_input2_to_hidden1_inh = new TransmitEventConnection(input_pyr2, hidden_pyr1, w_ff_inh, p_ff, GABA);
    con_input2_to_hidden1_inh->set_target("g_gaba");
    
    TransmitEventConnection * con_input2_to_hidden2_inh = new TransmitEventConnection(input_pyr2, hidden_pyr2, w_ff_inh, p_ff, GABA);
    con_input2_to_hidden2_inh->set_target("g_gaba");
    
    
    // (2) Connect hidden layer to output layer
    // Excitation
    AdaptiveEBCPConnection * con_hidden1_to_output_exc = new AdaptiveEBCPConnection(hidden_pyr1, output_pyr, w_ff_exc, p_ff, learning_rate, max_weight, tau_pre);
    fix_plastic_connections_parameters(con_hidden1_to_output_exc, moving_average_time_constant);
    
    AdaptiveEBCPConnection * con_hidden2_to_output_exc = new AdaptiveEBCPConnection(hidden_pyr2, output_pyr, w_ff_exc, p_ff, learning_rate, max_weight, tau_pre);
    fix_plastic_connections_parameters(con_hidden2_to_output_exc, moving_average_time_constant);
    
    // Inhibition
    TransmitEventConnection * con_hidden1_to_output_inh = new TransmitEventConnection(hidden_pyr1, output_pyr, w_ff_inh, p_ff, GABA);
    con_hidden1_to_output_inh->set_target("g_gaba");
    
    TransmitEventConnection * con_hidden2_to_output_inh = new TransmitEventConnection(hidden_pyr2, output_pyr, w_ff_inh, p_ff, GABA);
    con_hidden2_to_output_inh->set_target("g_gaba");
    

    //-- CONNECT FeedBACK
    float w_fb_exc = 0.12*4000/number_of_neurons;
    float w_fb_inh = 0.03*4000/number_of_neurons;
    float p_fb = p_ff;
    
    // Excitation
    TransmitBurstConnection * con_output_to_hidden1_exc = new TransmitBurstConnection(output_pyr, hidden_pyr1, w_fb_exc, p_fb, GLUT);
    con_output_to_hidden1_exc->set_target("g_ampa_dend");
    TransmitBurstConnection * con_output_to_hidden2_exc = new TransmitBurstConnection(output_pyr, hidden_pyr2, w_fb_exc, p_fb, GLUT);
    con_output_to_hidden2_exc->set_target("g_ampa_dend");
    
    // Inhibition
    TransmitBurstConnection * con_output_to_hidden1_inh = new TransmitBurstConnection(output_pyr, hidden_pyr1, w_fb_inh, p_fb, GABA);
    con_output_to_hidden1_inh->set_target("g_gaba_dend");
    TransmitBurstConnection * con_output_to_hidden2_inh = new TransmitBurstConnection(output_pyr, hidden_pyr2, w_fb_inh, p_fb, GABA);
    con_output_to_hidden2_inh->set_target("g_gaba_dend");
    
    
    /**************************************************/
    /******         CURRENT INJECTORS         *********/
    /**************************************************/
    CurrentInjector * baseline_curr_output_pyr_dend  = new CurrentInjector(output_pyr, "Vd");
    CurrentInjector * curr_inject_output_pyr_dend  = new CurrentInjector(output_pyr, "Vd");
    CurrentInjector * curr_inject_hidden_pyr1_dend = new CurrentInjector(hidden_pyr1, "Vd");
    CurrentInjector * curr_inject_hidden_pyr2_dend = new CurrentInjector(hidden_pyr2, "Vd");

    CurrentInjector * curr_inject_output_pyr_soma  = new CurrentInjector(output_pyr, "mem");
    CurrentInjector * curr_inject_hidden_pyr1_soma = new CurrentInjector(hidden_pyr1, "mem");
    CurrentInjector * curr_inject_hidden_pyr2_soma = new CurrentInjector(hidden_pyr2, "mem");
    CurrentInjector * curr_inject_input_pyr1_soma  = new CurrentInjector(input_pyr1, "mem");
    CurrentInjector * curr_inject_input_pyr2_soma  = new CurrentInjector(input_pyr2, "mem");
    
    baseline_curr_output_pyr_dend->set_all_currents(400.e-12/output_pyr[0].get_Cd());
    curr_inject_hidden_pyr1_soma->set_all_currents(-100.e-12/hidden_pyr1[0].get_Cs());
    curr_inject_hidden_pyr2_soma->set_all_currents(-100.e-12/hidden_pyr2[0].get_Cs());
    curr_inject_output_pyr_soma->set_all_currents(-200.e-12/output_pyr[0].get_Cs());

    /**************************************************/
    /******              MONITORS             *********/
    /**************************************************/
    double binSize_rate = 1.;
    double binSize_wsum = binSize_rate;
    
    // Burst/event rate monitors
    BurstRateMonitor * brmon_input1  = new BurstRateMonitor( input_pyr1, sys->fn("brate_input1"), binSize_rate);
    BurstRateMonitor * brmon_input2  = new BurstRateMonitor( input_pyr2, sys->fn("brate_input2"), binSize_rate);
    BurstRateMonitor * brmon_hidden1 = new BurstRateMonitor( hidden_pyr1, sys->fn("brate_hidden1"), binSize_rate);
    BurstRateMonitor * brmon_hidden2 = new BurstRateMonitor( hidden_pyr2, sys->fn("brate_hidden2"), binSize_rate);
    BurstRateMonitor * brmon_output  = new BurstRateMonitor( output_pyr, sys->fn("brate_output"), binSize_rate);

    // Weight sum monitors
    WeightSumMonitor * wsmon_in1_to_hid1 = new WeightSumMonitor( con_input1_to_hidden1_exc, sys->fn("wsum_in1_to_hid1"), binSize_wsum );
    WeightSumMonitor * wsmon_in1_to_hid2 = new WeightSumMonitor( con_input1_to_hidden2_exc, sys->fn("wsum_in1_to_hid2"), binSize_wsum );
    WeightSumMonitor * wsmon_in2_to_hid1 = new WeightSumMonitor( con_input2_to_hidden1_exc, sys->fn("wsum_in2_to_hid1"), binSize_wsum );
    WeightSumMonitor * wsmon_in2_to_hid2 = new WeightSumMonitor( con_input2_to_hidden2_exc, sys->fn("wsum_in2_to_hid2"), binSize_wsum );
    WeightSumMonitor * wsmon_hid1_to_out = new WeightSumMonitor( con_hidden1_to_output_exc, sys->fn("wsum_hid1_to_out"), binSize_wsum );
    WeightSumMonitor * wsmon_hid2_to_out = new WeightSumMonitor( con_hidden2_to_output_exc, sys->fn("wsum_hid2_to_out"), binSize_wsum );

    //**********************************************************************/
    //***                          SIMULATIONS                           ***/
    //**********************************************************************/
    const double max_input_current = 200e-12;
    const double min_input_current = -200e-12;
    
    input_output_XOR example(max_input_current, min_input_current);
    const int number_of_training_examples = 10;
    const int number_of_examples_before_test = 25;
    const int number_of_examples_before_cost = 5;

    double Cs_input_pyr1 = input_pyr1[0].get_Cs();
    double Cs_input_pyr2 = input_pyr2[0].get_Cs();
    double Cd_output_pyr = output_pyr[0].get_Cd();
    
    // burn-in period
    sys->run(3*moving_average_time_constant);
    
    for (int i=1;i<=number_of_training_examples;i++){
        // select example at random
        example.select_XOR_example();
        
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
        
        // run for half the duration of the example
        sys->run(durex/2.);
        
        // inject current in output pyramids' dendrites to generate error
        generate_error(con_hidden1_to_output_exc, curr_inject_output_pyr_dend, example.output, Cd_output_pyr);
        
        // turn plasticity on
        con_input1_to_hidden1_exc->stdp_active = true;
        con_input1_to_hidden2_exc->stdp_active = true;
        con_input2_to_hidden1_exc->stdp_active = true;
        con_input2_to_hidden2_exc->stdp_active = true;
        
        con_hidden1_to_output_exc->stdp_active = true;
        con_hidden2_to_output_exc->stdp_active = true;
        
        // propagate credit and learn
        sys->run(durex/2.);
        
        // compute cost
        if ( i % number_of_examples_before_cost == 0 ){
            double cost = 0;
            con_input1_to_hidden1_exc->stdp_active = false;
            con_input1_to_hidden2_exc->stdp_active = false;
            con_input2_to_hidden1_exc->stdp_active = false;
            con_input2_to_hidden2_exc->stdp_active = false;
            
            con_hidden1_to_output_exc->stdp_active = false;
            con_hidden2_to_output_exc->stdp_active = false;
            
            // input (0,0)
            curr_inject_input_pyr1_soma->set_all_currents(min_input_current/Cs_input_pyr1);
            curr_inject_input_pyr2_soma->set_all_currents(min_input_current/Cs_input_pyr2);
            
            cost += compute_cost(con_hidden1_to_output_exc, 0);
            sys->run(durex/2.);
            
            // input (0,1)
            curr_inject_input_pyr1_soma->set_all_currents(min_input_current/Cs_input_pyr1);
            curr_inject_input_pyr2_soma->set_all_currents(max_input_current/Cs_input_pyr2);
            
            cost += compute_cost(con_hidden1_to_output_exc, 1);
            sys->run(durex/2.);
            
            // input (1,0)
            curr_inject_input_pyr1_soma->set_all_currents(max_input_current/Cs_input_pyr1);
            curr_inject_input_pyr2_soma->set_all_currents(min_input_current/Cs_input_pyr2);
            
            cost += compute_cost(con_hidden1_to_output_exc, 1);
            sys->run(durex/2.);
            
            // input (1,1)
            curr_inject_input_pyr1_soma->set_all_currents(max_input_current/Cs_input_pyr1);
            curr_inject_input_pyr2_soma->set_all_currents(max_input_current/Cs_input_pyr2);
            
            cost += compute_cost(con_hidden1_to_output_exc, 0);
            sys->run(durex/2.);
            
            std::cout<<"example #"<<i <<": cost = "<<cost<<std::endl;

        }
        
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
/*
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
}*/

void generate_error(AdaptiveEBCPConnection * connect,  CurrentInjector * inject_dend, int output, float capacitance){
    Trace * tmp_event_trace = connect->get_tr_event();
    AurynFloat injected_current;
    AurynFloat er_estimate;
    
    for (NeuronID i=0; i<number_of_neurons;i++){
        er_estimate = tmp_event_trace->normalized_get(i);
        if ( output == 0 ){
            if (er_estimate < e_max - 0.5) injected_current = (-500e-12/0.4)/(e_max - er_estimate) ;
            else injected_current = (-500e-12/0.4)/0.5;
            inject_dend->set_current(i, injected_current/capacitance);
        }
        else {
            if (abs(er_estimate - e_max) < 0.1) {injected_current = 0;}
            else {injected_current = (500e-12/0.4)/er_estimate;}
            inject_dend->set_current(i, injected_current/capacitance);
        }
    }
    
}


double compute_cost(AdaptiveEBCPConnection * connect, int output){
    // compute mean event rate over population
    double mean_er = 0.;
    Trace * tmp_event_trace = connect->get_tr_event();

    for (NeuronID i=0; i<number_of_neurons;i++){
        mean_er += tmp_event_trace->normalized_get(i);
    }
    mean_er /= number_of_neurons;
    
    // compute cost
    if (mean_er > e_max) {
        mean_er = 0.99*e_max;
    }
    double cost = -output*log(mean_er/e_max) - (1 - output)*log(1. - mean_er/e_max);

    return cost;
    
}

void fix_plastic_connections_parameters(AdaptiveEBCPConnection * connect, double alpha){
    connect->set_post_trace_event_tau(alpha);
    connect->set_post_trace_burst_tau(alpha);
    connect->set_target("g_ampa");
    connect->max_rate = e_max;
}
