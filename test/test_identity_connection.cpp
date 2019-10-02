/***********************************************************
 Description: This program tests the identity connection with
 BiasConnection
************************************************************/

#include "auryn.h"
#include "BiasIdentityConnection.h"

using namespace auryn;

namespace po = boost::program_options;

int main(int ac, char* av[])
{
    int errcode = 0;
    int seed = 123;
    double w0 = 0.06;
    double d0 = 0.2;
    double kappa = 100.0;
    double tau_pre = 20.e-3;
    double simtime = 1800;
    double moving_average_time_constant = 10.;
    string simname = "test_identity_connection";
    std::string dir = ".";
    
    
    
    
    // INITIALIZE AURYN
    auryn_init( ac, av, dir, simname);
    sys->set_master_seed(seed);
    
    // Main neuron group
    NaudGroup * neurons_exc = new NaudGroup(1000);
    PoissonGroup * poisson = new PoissonGroup(1000, kappa);
    
    /*BCPConnection * con_ext_soma = new BCPConnection(poisson,
                                                       neurons_exc,
                                                       GLUT,
                                                       learning_rate,
                                                       max_weight, 20.e-3);
    */
    /*
    BiasConnection * con_ext_soma = new BiasConnection(input,
                                                       neurons_exc,
                                                       we_soma,
                                                       sparseness,
                                                       learning_rate,
                                                       max_weight);
    
     */
    BiasIdentityConnection * con_ext_soma = new BiasIdentityConnection(poisson, neurons_exc, w0, 0.1e-3, GLUT);
    
    //con_ext_soma->load_from_complete_file("ident_conn.wmat");
    //con_ext_soma->write_to_file("test.wmat");
    /*
    for (int j = 0; j < 1000; j++){
        for (int i=0;i<clist.size();++i){
            std::cout << clist[j][i].i << " " << clist[j][i].j << std::endl;
        }
    }*/
    //con_ext_soma->set_target("g_ampa");
    //con_ext_soma->set_post_trace_event_tau(moving_average_time_constant);
    //con_ext_soma->set_post_trace_burst_tau(moving_average_time_constant);
    //con_ext_soma->max_rate = 20.;
    //con_ext_soma->min_rate = 2.;
    

    //sys->run(0);
    if (errcode)
        auryn_abort(errcode);
    
    logger->msg("Freeing ...",PROGRESS,true);
    auryn_free();
    return errcode;
}
