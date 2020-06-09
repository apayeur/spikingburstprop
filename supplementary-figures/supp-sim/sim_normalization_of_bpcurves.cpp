/*
 * Figure S1: Test the effect of recurrent STF feedback on
 * the dendritic transfer function.
 */

#include "auryn.h"
#include <string>
#include "STPeTMConnection.h"

using namespace auryn;

namespace po = boost::program_options;

void fix_parameters_som_neurons(AdExGroup* som);
void set_Facilitating_connection(STPeTMConnection * pre_to_post);

int main(int ac, char* av[])
{
    int errcode = 0;
    char strbuf [255];
    unsigned int seed = 1;
    double somatic_current = 0e-12;
    AurynFloat w_pyr_to_pyr(0.02);
    string dir = "./";
    string simname = "normalization";
    
    //***          Options            ***//
    try {
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help", "produce help message")
        ("dir", po::value<string>(), "output directory")
        ("seed", po::value<int>(), "random seed")
        ("somacurrent", po::value<double>(), "somatic current")
        ("w_pyr_to_pyr", po::value<double>(), "fake feedback inhibition")
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
        
        if (vm.count("somacurrent")) {
            std::cout << "somatic current set to "
            << vm["somacurrent"].as<double>() << ".\n";
             somatic_current = vm["somacurrent"].as<double>();
        }
        
        if (vm.count("w_pyr_to_pyr")) {
            std::cout << "w_pyr_to_pyr (fake FB inh) set to "
            << vm["w_pyr_to_pyr"].as<double>() << ".\n";
             w_pyr_to_pyr = vm["w_pyr_to_pyr"].as<double>();
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
    auryn_init( ac, av, ".", simname );
    sys->set_output_dir(dir);

    //***          Neuronal populations            ***//
    // 2-comp neuron group
    NeuronID number_of_neurons = 4000;
    NaudGroup* pyr = new NaudGroup(number_of_neurons);
    pyr->syn_exc_soma->set_nmda_ampa_current_ampl_ratio(0.);
    pyr->syn_exc_dend->set_nmda_ampa_current_ampl_ratio(0.);
    pyr->random_mem(-70e-3, 5.e-3);
    
    // SOM neurons
    //AdExGroup* som = new AdExGroup(1000);
    //fix_parameters_som_neurons(som);
    
    // define external population of Poisson neurons
    NeuronID nb_exc_Poisson_neurons = 25000;
    NeuronID nb_inh_Poisson_neurons = nb_exc_Poisson_neurons;
    float poisson_rate_exc = 5.;
    float poisson_rate_inh = 5.;
    PoissonGroup* exc_Poisson = new PoissonGroup(nb_exc_Poisson_neurons, poisson_rate_exc);
    PoissonGroup* inh_Poisson = new PoissonGroup(nb_inh_Poisson_neurons, poisson_rate_inh);
    exc_Poisson->seed(1);
    inh_Poisson->seed(1);
    
    
    //***          Current injections            ***//
    CurrentInjector * curr_inject_soma = new CurrentInjector(pyr, "mem");
    CurrentInjector * curr_inject_dend = new CurrentInjector(pyr, "Vd");
    
    //***          Connectivity         ***//
    // These connections yield an avg ER = 4-5 Hz without injected current, with subthreshold
    // voltage fluctuations around 4 mV (Destexhe et al.). The BP at zero injected dendritic and somatic
    // currents is around 12%, and Vd fluctuations are about 2.5 mV when disregarding bAP and Ca spikes (thresholding the distribution std(Vd[Vd<-60mV])).
    int number_of_ext_conn[2] = {100, 30}; //exc, inh
    
    float p_exc = float(number_of_ext_conn[0])/float(nb_exc_Poisson_neurons);
    float p_inh = float(number_of_ext_conn[1])/float(nb_inh_Poisson_neurons);
    
    const float w_exc = 0.17;   // conductance amplitude in units of leak conductance
    const float w_inh = 0.26;
    
    SparseConnection * con_ext_exc_soma = new SparseConnection(exc_Poisson, pyr, w_exc, p_exc, GLUT);
    con_ext_exc_soma->set_target("g_ampa");
    SparseConnection * con_ext_inh_soma = new SparseConnection(inh_Poisson, pyr, w_inh, p_inh, GABA);
    con_ext_inh_soma->set_target("g_gaba");
    
    const float w_exc_dend = 0.043;        // conductance amplitude in units of leak conductance
    const float w_inh_dend = 0.07;
    
    SparseConnection * con_ext_exc_dend = new SparseConnection(exc_Poisson, pyr, w_exc_dend, p_exc, GLUT);
    con_ext_exc_dend->set_target("g_ampa_dend");
    SparseConnection * con_ext_inh_dend = new SparseConnection(inh_Poisson, pyr, w_inh_dend, p_inh, GABA);
    con_ext_inh_dend->set_target("g_gaba_dend");
     
    // Pyr to Pyr - fake inh STF
    float p_pyr_to_pyr = 0.05;
    STPeTMConnection * pyr_to_pyr = new STPeTMConnection(pyr, pyr, w_pyr_to_pyr, p_pyr_to_pyr, GABA);
    set_Facilitating_connection(pyr_to_pyr);
    pyr_to_pyr->set_target("g_gaba_dend");
    
    
    // Uncomment the next two blocks and comment the previous one (i.e. Pyr to Pyr - fake inh STF)
    // to test the SOM-mediated feeback
    
    // External Poisson neurons -> SOM
    /*float p_exc_to_som = 100./float(nb_exc_Poisson_neurons);
    SparseConnection * con_ext_exc_som = new SparseConnection(exc_Poisson, som, 0.2, p_exc_to_som, GLUT);
    SparseConnection * con_ext_inh_som = new SparseConnection(inh_Poisson, som, 0.3, p_exc_to_som, GABA);
    */
    
    // SOM to Pyr
    //AllToAllConnection * som_to_pyr = new AllToAllConnection(som, pyr, 3*0.001, GABA);
    //som_to_pyr->set_target("g_gaba_dend");

    
    //***          Monitors            ***//
    double binSize_rate = 20.e-3; // ms
    auto label_str = std::to_string(somatic_current*1.e12);
    auto label_w = std::to_string(w_pyr_to_pyr);

    //PopulationRateMonitor * pmon = new PopulationRateMonitor( pyr, sys->fn("prate_somacurrent"+label_str), binSize_rate );
    //PopulationRateMonitor * somrate_mon = new PopulationRateMonitor( som, sys->fn("somrate"), binSize_rate );
    BurstRateMonitor * brmon = new BurstRateMonitor( pyr, sys->fn("somacurrent"+label_str+"_wpyrtopyr"+label_w,"brate"), binSize_rate);
    //SpikeMonitor * smon = new SpikeMonitor( pyr, sys->fn("ras_somacurrent"+label_str));
    
    //***          Simulation            ***//
    logger->msg("Running ...",PROGRESS);

    // current intensities
    const double max_dendritic_current  = 350e-12;
    const double min_dendritic_current  = -100e-12;
    const double dendritic_current_incr = 50e-12;
    const double simtime = 1000e-3;
    
    
    double dendritic_current = min_dendritic_current;
    curr_inject_soma->set_all_currents(somatic_current/pyr[0].get_Cs());
    
    while (dendritic_current < max_dendritic_current + dendritic_current_incr/2)
    {
        curr_inject_dend->set_all_currents(dendritic_current/pyr[0].get_Cd());
        dendritic_current += dendritic_current_incr;
        sys->run(simtime);
    }

    
    if (errcode)
        auryn_abort(errcode);
    
    logger->msg("Freeing ...",PROGRESS,true);
    auryn_free();
    return errcode;
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

void set_Facilitating_connection(STPeTMConnection * pre_to_post){
    pre_to_post->set_tau_d(100e-3);
    pre_to_post->set_tau_f(100e-3);
    pre_to_post->set_ujump(0.02);
    pre_to_post->set_f(0.1);
}
