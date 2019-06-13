//
//  BurstPoissonGroup.cpp
//  
//
//  Created by Alexandre Payeur on 4/7/19.
//

#include "BurstPoissonGroup.h"

using namespace auryn;

boost::mt19937 BurstPoissonGroup::gen = boost::mt19937();

BurstPoissonGroup::BurstPoissonGroup( NeuronID size, NodeDistributionMode distmode ) : NaudGroup(size, distmode)
{
    set_name("BurstPoissonGroup");
    if ( evolve_locally() ) {
        dist = new boost::random::uniform_real_distribution<>(0., 1.);
        die  = new boost::variate_generator<boost::mt19937&, boost::random::uniform_real_distribution<> > ( gen, *dist );
        salt = sys->get_seed();
        seed(sys->get_seed());
        
        refractory_state = new AurynVector<unsigned int>(get_vector_size());
        refractory_state->set_all(0);
        burst_state = new AurynVector<unsigned int>(get_vector_size());
        burst_state->set_all(0);
        
        abs_ref_period = int(20.e-3/auryn_timestep);
        burst_duration = int(10.e-3/auryn_timestep);
        this->e_dend = -57e-3;
    }
}

BurstPoissonGroup::~BurstPoissonGroup()
{
    if ( evolve_locally() ) {
        delete dist;
        delete die;
        delete refractory_state;
        delete burst_state;
    }
}

/*! \brief This method applies the Euler integration step to the membrane dynamics. */
void BurstPoissonGroup::integrate_membrane()
{
    // somatic dynamics
    t_leak->diff(e_rest, state_soma); // leak current
    state_soma->saxpy(mul_soma, t_leak);
    state_soma->saxpy(mul_wsoma, state_wsoma);
    
    // dendritic dynamics
    t_leak->diff(e_rest, state_dend);  // dendritic leak
    state_dend->saxpy(mul_dend, t_leak);
    state_dend->saxpy(mul_wdend, state_wdend);
    
    // adjust refractory and burst state
    for ( NeuronID i = 0 ; i < get_post_size() ; ++i ) {
        if ( refractory_state->get(i) ) {
            // decrease refractory state counter
            refractory_state->add_specific(i,-1);
        }
        if ( burst_state->get(i) ) {
            // decrease burst state
            burst_state->add_specific(i,-1);
        }
    }
    
    // synaptic input
    state_soma->saxpy(mul_soma, syn_current_exc_soma);
    state_soma->saxpy(mul_soma, syn_current_inh_soma);
    state_dend->saxpy(mul_dend, syn_current_exc_dend);
    state_dend->saxpy(mul_dend, syn_current_inh_dend);
    
    // update somatic adapation variable (S42)
    state_wsoma->scale(scale_wsoma);
    
    // update dendritic adaptation variable (S44)
    temp->diff(state_dend,e_rest);  // dendritic leak
    temp->saxpy(-1.0,state_wdend);
    state_wdend->saxpy(aux_mul_wdend,temp);
    
    // decay moving threshold
    thr->follow_scalar(e_thr, mul_thr);
}


void BurstPoissonGroup::check_thresholds()
{
    temp->sigmoid(state_dend, xi, e_dend ); // burst probability
    
    for ( AurynState * i = mem->data ; i != mem->data+get_rank_size() ; ++i ) { // it's important to use rank_size here otherwise there might be spikes from units that do not exist
        NeuronID unit = i-mem->data;
        if ( *i > thr->get(unit) and !refractory_state->get(unit) ) {
            
            push_spike(unit);
            
            refractory_state->set(unit, abs_ref_period);
            mem->set( unit, e_reset); // reset
            thr->add_specific( unit, e_spk_thr); // not needed
            state_wsoma->add_specific(unit, 1.0); // increments somatic adaptation variable
            
            // check whether a burst occurs
            AurynDouble r = (*die)();
            if ( r<temp->get(unit) ){
                burst_state->set(unit, burst_duration);
            }
        }
        else if (burst_state->get(unit)==1){
            push_spike(unit);
            refractory_state->set(unit, abs_ref_period);
            // commenting the next three lines mean that there is no soma-dendrite coupling
            //mem->set( unit, e_reset);
            //thr->add_specific( unit, e_spk_thr);
            //state_wsoma->add_specific(unit, 1.0);
        }
    }
    
}

void BurstPoissonGroup::evolve()
{
    syn_exc_soma->evolve(); //!< integrate_linear_nmda_synapses
    syn_inh_soma->evolve();
    syn_exc_dend->evolve(); //!< integrate_linear_nmda_synapses
    syn_inh_dend->evolve();
    integrate_membrane();
    check_thresholds();
}

void BurstPoissonGroup::seed(unsigned int s)
{
    std::stringstream oss;
    oss << "BurstPoissonGroup:: Seeding with " << s
    << " and " << salt << " salt";
    auryn::logger->msg(oss.str(),NOTIFICATION);
    
    gen.seed( s + salt );
}

