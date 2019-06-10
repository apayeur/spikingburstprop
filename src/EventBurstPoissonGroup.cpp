//
//  EventBurstPoissonGroup.cpp
//  
//
//  Created by Alexandre Payeur on 4/7/19.
//

#include "EventBurstPoissonGroup.h"

using namespace auryn;

//boost::mt19937 EventBurstPoissonGroup::gen = boost::mt19937();

EventBurstPoissonGroup::EventBurstPoissonGroup( NeuronID size, NodeDistributionMode distmode ) : BurstPoissonGroup(size, distmode)
{
    set_name("EventBurstPoissonGroup");
    
    this->e_thr = -52e-3;
    this->e_spk_thr = 2e-3;

    // declare temp vector
    t_prob_dend = get_state_vector("_prob_dend");
    t_prob_soma = get_state_vector("_prob_soma");

}

EventBurstPoissonGroup::~EventBurstPoissonGroup()
{
    if ( evolve_locally() ) {
        delete t_prob_dend;
        delete t_prob_soma;
    }
}

/*! \brief This method applies the Euler integration step to the membrane dynamics. */
void EventBurstPoissonGroup::integrate_membrane()
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


void EventBurstPoissonGroup::check_thresholds()
{
    // burst probability
    t_prob_dend->sigmoid(state_dend, xi, e_dend ); // burst probability
    
    // somatic spike probability = 1 - exp(-auryn_timestep*exp(slope*(V-thr)))
    t_prob_soma->diff(state_soma,e_thr);
    t_prob_soma->mul(1./e_spk_thr);
    t_prob_soma->exp();
    t_prob_soma->mul(-auryn_timestep);
    t_prob_soma->exp();
    t_prob_soma->mul(-1.);
    t_prob_soma->add(1.0);
    
    
    for ( AurynState * i = mem->data ; i != mem->data+get_rank_size() ; ++i ) { // it's important to use rank_size here otherwise there might be spikes from units that do not exist
        NeuronID unit = i-mem->data;
        // check whether a somatic spike occurs
        AurynDouble r = (*die)();
        if ( r<t_prob_soma->get(unit) and !refractory_state->get(unit) ) {
            
            push_spike(unit);
            
            refractory_state->set(unit, abs_ref_period);
            mem->set( unit, e_reset); // reset
            thr->add_specific( unit, e_spk_thr); // not needed
            state_wsoma->add_specific(unit, 1.0); // increments somatic adaptation variable
            
            // check whether a burst occurs
            r = (*die)();
            if ( r<t_prob_dend->get(unit) ){
                burst_state->set(unit, burst_duration);
            }
        }
        else if (burst_state->get(unit)==1){
            push_spike(unit);
            refractory_state->set(unit, abs_ref_period);
            mem->set( unit, e_reset);
            thr->add_specific( unit, e_spk_thr);
            state_wsoma->add_specific(unit, 1.0);
        }
    }
    
}

void EventBurstPoissonGroup::evolve()
{
    syn_exc_soma->evolve(); //!< integrate_linear_nmda_synapses
    syn_inh_soma->evolve();
    syn_exc_dend->evolve(); //!< integrate_linear_nmda_synapses
    syn_inh_dend->evolve();
    integrate_membrane();
    check_thresholds();
}
