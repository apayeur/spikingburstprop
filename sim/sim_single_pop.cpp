/************************************************************************
 Description: Reproduces figure 1 in Naud & Sprekeler (2018),
 with Poisson input noise.
 
 The noise parameters (average number of connections, weights, etc.) are
 chosen so that the average event rate is ~4-5 Hz without injected current,
 with subthreshold voltage fluctuations around 4 mV (Destexhe et al.).
 
 The burst probability at zero injected dendritic and somatic currents
 is around 12%, and Vdend fluctuations are about 2.5 mV when disregarding
 bAP and Ca spikes ( thresholding the distribution std(Vd[Vd<-60mV]) ).
 ************************************************************************/

#include "auryn.h"
#include "STPeTMConnection.h"

using namespace auryn;

namespace po = boost::program_options;

void initialize_pyr_neurons(NaudGroup* pyr);
void set_Facilitating_connection(STPeTMConnection * pre_to_post);

int main(int ac, char* av[])
{
    int errcode = 0;
    char strbuf [255];
    unsigned int seed = 1;
    string dir = "../data/single-pop";
    string simname = "single_pop";
    float w_epois_to_soma = 0.35;
    float w_epois_to_dend = 0.05;

    /**************************************************/
    /********              OPTIONS            *********/
    /**************************************************/
    try {
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help", "produce help message")
        ("dir", po::value<string>(), "output directory")
        ("seed", po::value<int>(), "random seed")
        ("w_epois_to_soma", po::value<float>(), "weight from exc poisson to soma")
        ("w_epois_to_dend", po::value<float>(), "weight from exc poisson to dendrite")
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
        
        if (vm.count("w_epois_to_soma")) {
            std::cout << "w_epois_to_soma set to "
            << vm["w_epois_to_soma"].as<float>() << ".\n";
            w_epois_to_soma = vm["w_epois_to_soma"].as<float>();
        }
        
        if (vm.count("w_epois_to_dend")) {
            std::cout << "w_epois_to_dend set to "
            << vm["w_epois_to_dend"].as<float>() << ".\n";
            w_epois_to_dend = vm["w_epois_to_dend"].as<float>();
        }
    }
    catch(std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        std::cerr << "Exception of unknown type!\n";
    }
    // Initialize auryn
    auryn_init( ac, av, dir, simname );

    
    /**************************************************/
    /******          NEURON POPULATIONS       *********/
    /**************************************************/
    
    // define 2-comp neuron group
    NeuronID number_of_neurons = 4000;
    NaudGroup* pyr = new NaudGroup(number_of_neurons);
    initialize_pyr_neurons(pyr);

    // define external population of Poisson neurons
    float poisson_rate = 5.;
    PoissonGroup* exc_Poisson      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson      = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* exc_Poisson_dend = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    PoissonGroup* inh_Poisson_dend = new PoissonGroup(number_of_neurons, 100*poisson_rate);
    exc_Poisson->seed(seed);
    inh_Poisson->seed(seed);
    exc_Poisson_dend->seed(seed);
    inh_Poisson_dend->seed(seed);
    
    
    /**************************************************/
    /******           CONNECTIVITY            *********/
    /**************************************************/
    float ratio_ie_soma = 1.;
    IdentityConnection * epois_to_soma = new IdentityConnection(exc_Poisson, pyr, w_epois_to_soma, GLUT);
    epois_to_soma->set_target("g_ampa");
    IdentityConnection * ipois_to_soma = new IdentityConnection(inh_Poisson, pyr, ratio_ie_soma*w_epois_to_soma, GABA);
    ipois_to_soma->set_target("g_gaba");
    
    float ratio_ie_dend = 1.;
    IdentityConnection * epois_to_dend = new IdentityConnection(exc_Poisson_dend, pyr, w_epois_to_dend, GLUT);
    epois_to_dend->set_target("g_ampa_dend");
    IdentityConnection * ipois_to_dend = new IdentityConnection(inh_Poisson_dend, pyr, ratio_ie_dend*w_epois_to_dend, GABA);
    ipois_to_dend->set_target("g_gaba_dend");

    float w_pyr_to_pyr = 0.1; //0.05;
    float p_pyr_to_pyr = 0.05;
    STPeTMConnection * pyr_to_pyr = new STPeTMConnection(pyr, pyr, w_pyr_to_pyr, p_pyr_to_pyr, GABA);
    set_Facilitating_connection(pyr_to_pyr);
    pyr_to_pyr->set_target("g_gaba_dend");
    
    
    /**************************************************/
    /******         CURRENT INJECTORS         *********/
    /**************************************************/
    CurrentInjector * curr_inject_soma = new CurrentInjector(pyr, "mem");
    CurrentInjector * curr_inject_dend = new CurrentInjector(pyr, "Vd");


    /**************************************************/
    /******              MONITORS             *********/
    /**************************************************/
    double binSize_rate = 20.e-3; // ms
    auto seed_str = std::to_string(seed);
    
    SpikeMonitor * smon          = new SpikeMonitor( pyr, sys->fn("ras_seed"+seed_str) );
    VoltageMonitor * vmon        = new VoltageMonitor( pyr, 0, sys->fn("mem") );
    VoltageMonitor * vmon1       = new VoltageMonitor( pyr, 1, sys->fn("mem1") );
    StateMonitor * smon_vd       = new StateMonitor( pyr, 0, "Vd", sys->fn("Vd") );
    StateMonitor * smon_vd1      = new StateMonitor( pyr, 1, "Vd", sys->fn("Vd1") );
    PopulationRateMonitor * pmon = new PopulationRateMonitor( pyr, sys->fn("prate"), binSize_rate );
    BurstRateMonitor * brmon     = new BurstRateMonitor( pyr, sys->fn("brate_seed"+seed_str), binSize_rate);
    //StateMonitor * smon_ws  = new StateMonitor( pyr, 0, "wsoma", sys->fn("wsoma") );
    
    
    /**************************************************/
    /******             SIMULATION            *********/
    /**************************************************/
    logger->msg("Running ...",PROGRESS);
    
    const double max_dendritic_current = 100e-12;
    const double min_dendritic_current = 0e-12;
    const double max_somatic_current   = 150e-12;
    const double min_somatic_current   = 0e-12;

    const double period = 200e-3;
    const double simtime = 1000e-3;
    const double segtime_maxsoma = period/2.;
    const double segtime_minsoma = period/2;
    const double segtime_maxdend = 0.7*period;
    const double segtime_mindend = 0.3*period;
    const double small_overlap = 0.1*period;

    // Relaxation period
    curr_inject_soma->set_all_currents(min_somatic_current/pyr[0].get_Cs());
    curr_inject_dend->set_all_currents(min_dendritic_current/pyr[0].get_Cd());
    sys->run(simtime);
    
    // Alternating currents
    for (int i=0;i<10;i++)
    {
        curr_inject_soma->set_all_currents(max_somatic_current/pyr[0].get_Cs());
        curr_inject_dend->set_all_currents(min_dendritic_current/pyr[0].get_Cd());
        sys->run(segtime_mindend - small_overlap);

        curr_inject_soma->set_all_currents(max_somatic_current/pyr[0].get_Cs());
        curr_inject_dend->set_all_currents(max_dendritic_current/pyr[0].get_Cd());
        sys->run(segtime_maxsoma  - (segtime_mindend - small_overlap) );
        
        curr_inject_soma->set_all_currents(min_somatic_current/pyr[0].get_Cs());
        curr_inject_dend->set_all_currents(max_dendritic_current/pyr[0].get_Cd());
        sys->run(segtime_maxdend + segtime_mindend - small_overlap - segtime_maxsoma);
        
        curr_inject_soma->set_all_currents(min_somatic_current/pyr[0].get_Cs());
        curr_inject_dend->set_all_currents(min_dendritic_current/pyr[0].get_Cd());
        sys->run(small_overlap);
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

void set_Facilitating_connection(STPeTMConnection * pre_to_post){
    pre_to_post->set_tau_d(100e-3);
    pre_to_post->set_tau_f(100e-3);
    pre_to_post->set_ujump(0.02);
    pre_to_post->set_f(0.1);
}
