/*******************************************
 * This test compute the burst fraction as
 * a function of the dendritic current for
 * the BurstPoisson model.
 ******************************************/
#include "auryn.h"
#include <string>
#include "BurstPoissonGroup.h"

using namespace auryn;

namespace po = boost::program_options;

int main(int ac, char* av[])
{
    int errcode = 0;
    char strbuf [255];
    unsigned int seed = 1;
    string dir = "./test-output/";
    string simname = "ficurve_burstpoisson";
    double somatic_current = 100e-12;

    //***          Options            ***//
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
    auryn_init( ac, av, dir, simname );
    
    //***          Neuronal populations            ***//
    // 2-comp neuron group
    NeuronID number_of_neurons = 4000;
    BurstPoissonGroup* pyr = new BurstPoissonGroup(number_of_neurons);
    pyr->syn_exc_soma->set_nmda_ampa_current_ampl_ratio(0.);
    pyr->syn_exc_dend->set_nmda_ampa_current_ampl_ratio(0.);
    pyr->random_mem(-70e-3, 5.e-3);
    
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
    
    // specify connectivity
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
    
    
    //***          Monitors            ***//
    double binSize_rate = 20.e-3; // ms
    auto label_str = std::to_string(somatic_current*1.e12);

    PopulationRateMonitor * pmon = new PopulationRateMonitor( pyr, sys->fn("prate_somacurrent"+label_str), binSize_rate );
    BurstRateMonitor * brmon = new BurstRateMonitor( pyr, sys->fn("brate_somacurrent"+label_str), binSize_rate);
    
    
    //***          Simulation            ***//
    logger->msg("Running ...",PROGRESS);

    // current intensities
    const double max_dendritic_current  = 1500e-12;
    const double min_dendritic_current  = -200e-12;
    const double dendritic_current_incr = 100e-12;
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





