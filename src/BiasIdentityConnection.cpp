//
//  BiasIdentityConnection.cpp
//  
//
//  Created by Alexandre Payeur on 9/9/19.
//

#include "BiasIdentityConnection.h"

using namespace auryn;

void BiasIdentityConnection::init(AurynFloat eta, AurynWeight w_init)
{
    if ( dst->get_post_size() == 0 ) return;
    
    // learning rate
    eta_ = eta;
    auryn::logger->parameter("eta",eta);

    // homeostasis constants (min and max event rate)
    min_rate = 2.;
    max_rate = 12.;
    
    // burst detector
    const double default_tau = 16e-3;
    tr_post = dst->get_post_trace(default_tau);
    burst_state = new AurynVector<unsigned short>(dst->get_vector_size() );
    burst_state->set_all(0);
    burst_thr = std::exp(-1.0);
    
    // weight vector
    weights = new AurynVectorFloat( dst->get_vector_size() );
    weights->set_all(w_init);
    
    // postsynaptic event and burst traces
    tr_event = new EulerTrace(dst->get_vector_size(), 5.); //!< event rate average
    tr_burst = new EulerTrace(dst->get_vector_size(), 5.); //!< burst rate average
    
    // min/max weight
    set_min_weight(0.);
    set_max_weight(1.);
    
    offset = 0 ;
    every = 1;
    
    stdp_active = true;
}

BiasIdentityConnection::BiasIdentityConnection(SpikingGroup * source,
                                               NeuronGroup * destination,
                                               AurynWeight weight,
                                               AurynFloat eta,
                                               TransmitterType transmitter,
                                               string name) :
IdentityConnection(source, destination, weight, transmitter, name) {
    init(eta, weight);
}

void BiasIdentityConnection::finalize()
{
    IdentityConnection::finalize();
}

void BiasIdentityConnection::free()
{
    delete burst_state;
    delete tr_event;
    delete tr_burst;
}

BiasIdentityConnection::~BiasIdentityConnection()
{
    free();
}


void BiasIdentityConnection::propagate()
{
    propagate_forward();
    compute_burst_rate();
}

AurynWeight BiasIdentityConnection::on_pre(NeuronID post)
{
    const NeuronID ts = dst->global2rank(post); // only to be used for post traces
    // compute hom part of local error signal
    double hom = 0.0;
    const double rate = tr_event->normalized_get(ts);
    if (rate<min_rate) {
        hom = (min_rate-rate);
        // std::cout << "hom event" << std::endl;
    }
    AurynDouble dw = eta_*hom;
    return dw;
}


void BiasIdentityConnection::propagate_forward()
{
    SpikeContainer::const_iterator spikes_end = src->get_spikes()->end();
    SpikeContainer::const_iterator spikes_begin = src->get_spikes()->begin();
    AurynWeight * weight;
    
    for (SpikeContainer::const_iterator spike = spikes_begin; spike != spikes_end; ++spike ) {
        
        if ( *spike%every == 0 )
        {
            // IMPORTANT the use of safe_transmit
            // is important here since there is
            // no weight matrix with only the
            // corresponding columns
            weight = weights->ptr((*spike)/every + offset);
            safe_transmit( (*spike)/every + offset , *weight );
            
            // handle plasticity
            if ( stdp_active ) {
                // performs weight update
                *weight += on_pre(( *spike/every ) + offset);
                
                // clips weights
                if ( *weight > get_max_weight() ) *weight = get_max_weight();
                else
                    if ( *weight < get_min_weight() ) *weight = get_min_weight();
            }
        }
    }
}

void BiasIdentityConnection::propagate_backward(const NeuronID translated_post, const AurynState valence)
{
    if (stdp_active) {
        // compute hom part of local error signal
        double hom_neg = 0.0;
        const double estimated_event_rate = tr_event->normalized_get(translated_post);
        if (estimated_event_rate>max_rate) {
            // std::cout << "hom dep" << std::endl;
            hom_neg = -(estimated_event_rate - max_rate)/estimated_event_rate;
        }
        const NeuronID post = dst->rank2global(translated_post);
        
        // computes plasticity update
        AurynWeight * weight = weights->ptr(post);
        *weight += eta_*(valence+hom_neg);
        //std::cout<<std::setprecision(16)<<hom_neg<<std::endl;
        
        // clips weights
        if ( *weight > get_max_weight() ) *weight = get_max_weight();
        else
            if ( *weight < get_min_weight() ) *weight = get_min_weight();
    }
}

void BiasIdentityConnection::compute_burst_rate()
{
    // loop over all spikes to separate events from bursts
    SpikeContainer::const_iterator spk;
    for ( spk = dst->get_spikes_immediate()->begin() ;
         spk < dst->get_spikes_immediate()->end() ;
         ++spk ) {
        const NeuronID s = dst->global2rank(*spk);
        
        // detect first spike in bursts and non-burst spikes
        if ( tr_post->get(s) < burst_thr ) {
            
            burst_state->set(s,0);
            AurynFloat BP_estimate = tr_burst->normalized_get(s)/tr_event->normalized_get(s);
            if (BP_estimate < 1.)
            {
                propagate_backward(s, -BP_estimate);
            }
            else propagate_backward(s, -1.);
            // it is preferable to update tr_event after propagate backward
            tr_event->inc(s);
        }
        
        // detect second spike in burst
        if ( tr_post->get(s) > burst_thr && burst_state->get(s) == 0 )
        {
            tr_burst->inc(s);
            burst_state->set(s,1);
            propagate_backward(s, 1);
        }
    }
}

void BiasIdentityConnection::evolve()
{
    tr_event->evolve();
    tr_burst->evolve();
}

void BiasIdentityConnection::set_post_trace_event_tau(AurynFloat x){
    tr_event->set_timeconstant(x);
}

void BiasIdentityConnection::set_post_trace_burst_tau(AurynFloat x){
    tr_burst->set_timeconstant(x);
}

void BiasIdentityConnection::set_min_weight(AurynWeight value){
    wmin = value;
}

void BiasIdentityConnection::set_max_weight(AurynWeight value){
    wmax = value;
}

AurynWeight BiasIdentityConnection::get_min_weight(){
    return wmin;
}

AurynWeight BiasIdentityConnection::get_max_weight(){
    return wmax;
}

void BiasIdentityConnection::stats(AurynDouble &mean, AurynDouble &std, StateID zid){
    //Warning: This compute means on a given rank
    mean = weights->mean();
    std  = weights->std();
}

Trace * BiasIdentityConnection::get_tr_event()
{
    return tr_event;
}

Trace * BiasIdentityConnection::get_tr_burst()
{
    return tr_burst;
}
