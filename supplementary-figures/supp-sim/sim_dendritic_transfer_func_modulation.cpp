/* sim_dendritic_transfer_func_modulation.cpp
 *
 * Figure S4: Modulation of dendritic transfer function
 * Created by Alexandre Payeur.
 */

#include "auryn.h"
#include <cmath>
#include <fstream>
#include "STPeTMConnection.h"

using namespace auryn;

void fix_parameters_som_neurons(AdExGroup* som);
void fix_parameters_pv_neurons(AdExGroup* pv);
void set_Facilitating_connection(STPeTMConnection * pre_to_post);
void set_Depressing_connection(STPeTMConnection * pre_to_post);

namespace po = boost::program_options;

int main(int ac, char* av[])
{
    int errcode = 0;
    char strbuf [255];
    string dir = "./";
    string filename = "";
    string simname = "dendritictransfer";
    AurynFloat wsomtopyr(0.001);
    AurynFloat vip_current(0.e-12);
    AurynFloat releaseprob(0.02);
    AurynFloat burstiness(1200e-12);
    
    //***          Options            ***//
       try {
           po::options_description desc("Allowed options");
           desc.add_options()
           ("help", "produce help message")
           ("dir", po::value<string>(), "output directory")
           ("filename", po::value<string>(), "output file name")
           ("wsomtopyr", po::value<float>(), "weight from som to pyr")
           ("vipdisinhibition", po::value<float>(), "vip-like disinhibitory current")
           ("releaseprob", po::value<float>(), "release probability for connection PC to SOM")
           ("burstiness", po::value<float>(), "control propensity to burst - dendritic conductance")
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
           
           if (vm.count("filename")) {
               std::cout << "file name set to "
               << vm["filename"].as<string>() << ".\n";
               filename = vm["filename"].as<string>();
           }
           
           if (vm.count("wsomtopyr")) {
               std::cout << "weight from som to pyr set to "
               << vm["wsomtopyr"].as<float>() << ".\n";
               wsomtopyr = vm["wsomtopyr"].as<float>();
           }
           
           if (vm.count("vipdisinhibition")) {
               std::cout << "vip-like disinhibitory current set to "
               << vm["vipdisinhibition"].as<float>() << ".\n";
               vip_current = vm["vipdisinhibition"].as<float>();
           }
                
           if (vm.count("releaseprob")) {
               std::cout << "release probability set to "
               << vm["releaseprob"].as<float>() << ".\n";
               releaseprob = vm["releaseprob"].as<float>();
           }
           
           if (vm.count("burstiness")) {
               std::cout << "burstiness - dendritic conductance - set to "
               << vm["burstiness"].as<float>() << ".\n";
               burstiness = vm["burstiness"].as<float>();
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
    //**********************************************//
    //***          Neural populations            ***//
    //**********************************************//

    // 2-comp pyramidal neurons
    NeuronID number_of_neurons = 10000;
    NaudGroup* pyr = new NaudGroup(number_of_neurons);
    pyr->syn_exc_soma->set_nmda_ampa_current_ampl_ratio(0.);
    pyr->syn_exc_dend->set_nmda_ampa_current_ampl_ratio(0.);
    pyr->set_gd(burstiness);
    
    // PV neurons
    AdExGroup* pv = new AdExGroup(1000);
    fix_parameters_pv_neurons(pv);
    AurynFloat gL_pv = 10e-9;
    AurynFloat taum_pv = 10e-3;
    
    // SOM neurons
    AdExGroup* som = new AdExGroup(1000);
    fix_parameters_som_neurons(som);
    AurynFloat gL_som = 5.e-9;
    AurynFloat taum_som = 20e-3;
    
    // External population of Poisson neurons
    NeuronID nb_exc_Poisson_neurons = 25000;
    NeuronID nb_inh_Poisson_neurons = nb_exc_Poisson_neurons;
    float poisson_rate = 5.;
    PoissonGroup* exc_Poisson = new PoissonGroup(nb_exc_Poisson_neurons, poisson_rate);
    PoissonGroup* inh_Poisson = new PoissonGroup(nb_inh_Poisson_neurons, poisson_rate);
    exc_Poisson->seed(5);
    inh_Poisson->seed(5);
    
    //****************************************//
    //***          Connectivity            ***//
    //****************************************//
    // External Poisson neurons -> Pyr
    float p_ext_to_soma = 100./float(nb_exc_Poisson_neurons); // probability of connection
    float w_ext_to_soma = 0.05;    // conductance amplitude in units of leak conductance
    
    float p_ext_to_dend = 100./float(nb_exc_Poisson_neurons);
    float w_ext_to_dend = 0.05;
    
    SparseConnection * con_ext_exc_soma = new SparseConnection(exc_Poisson, pyr, w_ext_to_soma, p_ext_to_soma, GLUT);
    con_ext_exc_soma->set_target("g_ampa");
    SparseConnection * con_ext_inh_soma = new SparseConnection(inh_Poisson, pyr, w_ext_to_soma, p_ext_to_soma, GABA);
    con_ext_inh_soma->set_target("g_gaba");
    SparseConnection * con_ext_exc_dend = new SparseConnection(exc_Poisson, pyr, w_ext_to_dend, p_ext_to_dend, GLUT);
    con_ext_exc_dend->set_target("g_ampa_dend");
    SparseConnection * con_ext_inh_dend = new SparseConnection(inh_Poisson, pyr, w_ext_to_dend, p_ext_to_dend, GABA);
    con_ext_inh_dend->set_target("g_gaba_dend");
    
    // External Poisson neurons -> PV
    float p_ext_to_pv = 100./float(nb_exc_Poisson_neurons);
    float w_extexc_to_pv = 0.12;
    float w_extinh_to_pv = 0.15;
    
    SparseConnection * con_ext_exc_pv = new SparseConnection(exc_Poisson, pv, w_extexc_to_pv, p_ext_to_pv, GLUT);
    SparseConnection * con_ext_inh_pv = new SparseConnection(inh_Poisson, pv, w_extinh_to_pv, p_ext_to_pv, GABA);
    
    // External Poisson neurons -> SOM
    SparseConnection * con_ext_exc_som = new SparseConnection(exc_Poisson, som, 0.2, p_ext_to_soma, GLUT);
    SparseConnection * con_ext_inh_som = new SparseConnection(inh_Poisson, som, 0.3, p_ext_to_soma, GABA);
    
    // Pyr to PV - STD
    float w_pyr_to_pv = 0.01e-9/10.e-9;
    float p_pyr_to_pv = 0.2;
    STPeTMConnection * pyr_to_pv = new STPeTMConnection(pyr, pv, w_pyr_to_pv, p_pyr_to_pv, GLUT);
    set_Depressing_connection(pyr_to_pv);
    
    // Pyr to SOM - STF
    float w_pyr_to_som = 0.2e-9/5.e-9; //original: 0.2e-9/5.e-9;
    float p_pyr_to_som = 0.2;
    STPeTMConnection * pyr_to_som = new STPeTMConnection(pyr, som, w_pyr_to_som, p_pyr_to_som, GLUT);
    set_Facilitating_connection(pyr_to_som);
    pyr_to_som->set_ujump(releaseprob);
    
    // PV to SOM
    float w_pv_to_som = 0.002;
    AllToAllConnection * pv_to_som = new AllToAllConnection(pv, som, 0.002, GABA);
    
    // SOM to Pyr
    AllToAllConnection * som_to_pyr = new AllToAllConnection(som, pyr, wsomtopyr, GABA);
    som_to_pyr->set_target("g_gaba_dend");
    
    //******************************************//
    //***          Input currents            ***//
    //******************************************//
    CurrentInjector * curr_inject_soma = new CurrentInjector(pyr, "mem");
    CurrentInjector * curr_inject_dend = new CurrentInjector(pyr, "Vd");
    CurrentInjector * curr_inject_pv   = new CurrentInjector(pv, "mem");
    CurrentInjector * curr_inject_som  = new CurrentInjector(som, "mem");
    
    //************************************//
    //***          Monitors            ***//
    //************************************//
    //PopulationRateMonitor * pyr_activity = new PopulationRateMonitor( pyr, sys->fn("prate"), 25.e-3 );
    BurstRateMonitor * brmon = new BurstRateMonitor( pyr, sys->fn(filename, "brate"), 25.e-3); //to plot burst rate, u 1:2; to plot event rate, u 1:3
    //PopulationRateMonitor * pv_activity = new PopulationRateMonitor( pv, sys->fn("pvrate"), 25.e-3 );
    //SpikeMonitor * pv_raster = new SpikeMonitor( pv, sys->fn("pvras"));
    //PopulationRateMonitor * som_activity = new PopulationRateMonitor( som, sys->fn("somrate"), 25.e-3 );
    
    
    //**************************************//
    //***          Simulation            ***//
    //**************************************//
    logger->msg("Running ...",PROGRESS);
    
    // current intensities
    const double max_dendritic_current  = 500e-12;
    const double min_dendritic_current  = -50e-12;
    const double dendritic_current_incr = 50e-12;
    const double max_somatic_current    = 610e-12;
    const double min_somatic_current    = 560e-12;
    const double simtime = 1000e-3;
    
    double dendritic_current = min_dendritic_current;
    curr_inject_soma->set_all_currents(min_somatic_current/pyr[0].get_Cs());
    curr_inject_pv->set_all_currents(205e-12/(taum_pv*gL_pv));
    curr_inject_som->set_all_currents((-5e-12-vip_current)/(taum_som*gL_som));

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

void set_Facilitating_connection(STPeTMConnection * pre_to_post){
    pre_to_post->set_tau_d(100e-3);
    pre_to_post->set_tau_f(100e-3);
    pre_to_post->set_f(0.1);
}

void set_Depressing_connection(STPeTMConnection * pre_to_post){
    pre_to_post->set_tau_d(20.e-3);
    pre_to_post->set_tau_f(1.);
    pre_to_post->set_ujump(0.9);
    pre_to_post->set_f(0.1);
}
