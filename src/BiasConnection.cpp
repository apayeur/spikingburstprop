//
//  BiasConnection.cpp
//  
//
//  Created by Alexandre Payeur on 9/4/19.
//

#include "BiasConnection.h"

using namespace auryn;

BiasConnection::BiasConnection(SpikingGroup * source,
                               NeuronGroup * destination,
                               AurynWeight weight,
                               AurynFloat sparseness,
                               AurynFloat eta,
                               AurynFloat maxweight,
                               TransmitterType transmitter,
                               std::string name)
: AdaptiveEBCPConnection(source, destination, weight, sparseness, eta, maxweight, 20.e-3, transmitter, name)
{
}

BiasConnection::BiasConnection(SpikingGroup * source,
                               NeuronGroup * destination,
                               TransmitterType transmitter,
                               AurynFloat eta,
                               AurynFloat maxweight)
: AdaptiveEBCPConnection(source, destination, transmitter, eta, maxweight, 20.e-3)
{
}

BiasConnection::BiasConnection(SpikingGroup * source,
                               NeuronGroup * destination,
                               TransmitterType transmitter)
: AdaptiveEBCPConnection(source, destination, transmitter)
{
}

void BiasConnection::finalize() {
    AdaptiveEBCPConnection::finalize();
}

void BiasConnection::free()
{
}

BiasConnection::~BiasConnection()
{
}


void BiasConnection::propagate_forward()
{
    // loop over all spikes
    for (SpikeContainer::const_iterator spike = src->get_spikes()->begin() ; // spike = pre_spike
         spike != src->get_spikes()->end() ; ++spike ) {
        // loop over all postsynaptic partners
        for (const NeuronID * c = w->get_row_begin(*spike) ;
             c != w->get_row_end(*spike) ;
             ++c ) { // c = post index
            
            // transmit signal to target at postsynaptic neuron
            AurynWeight * weight = w->get_data_ptr(c);
            transmit( *c , *weight );
            
            /*
            // handle plasticity
            if ( stdp_active ) {
                // performs weight update
                *weight += on_pre(*c);
                
                // clips weights
                if ( *weight > get_max_weight() ) *weight = get_max_weight();
                else
                    if ( *weight < get_min_weight() ) *weight = get_min_weight();
            }*/
        }
    }
}



void BiasConnection::propagate_backward(const NeuronID translated_post, const AurynState valence)
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
        
        // loop over all presynaptic partners
        for (const NeuronID * c = bkw->get_row_begin(post) ; c != bkw->get_row_end(post) ; ++c ) {
            // computes plasticity update
            AurynWeight * weight = bkw->get_data(c);
            *weight += eta_*(valence+hom_neg);
            
            // clips weights
            if ( *weight > get_max_weight() ) *weight = get_max_weight();
            else
                if ( *weight < get_min_weight() ) *weight = get_min_weight();
        }
    }
}

void BiasConnection::propagate()
{
    propagate_forward();
    compute_burst_rate();
}

void BiasConnection::compute_burst_rate()
{
    // loop over all spikes to separate events from bursts
    SpikeContainer::const_iterator spk;
    for ( spk = dst->get_spikes_immediate()->begin() ;
         spk < dst->get_spikes_immediate()->end() ;
         ++spk ) {
        const NeuronID s = dst->global2rank(*spk);
        
        // detect first spike in bursts and non-burst spikes
        if ( tr_post->get(s) < burst_thr ) {
            
            tr_event_aux->set(s, tr_event->get(s));
            burst_state->set(s,0);
            if (tr_event->get(s) > 0) propagate_backward(s, -tr_burst->normalized_get(s)/tr_event->normalized_get(s));
            else propagate_backward(s, -1.);
            // it is preferable to update tr_event after propagate backward
            tr_event->inc(s);
        }
        
        // detect second spike in burst
        if ( tr_post->get(s) > burst_thr && burst_state->get(s) == 0 ) {
            tr_burst->inc(s);
            burst_state->set(s,1);
            propagate_backward(s, 1);
        }
    }
}

void BiasConnection::construct_identity_connection(AurynWeight weight_value){
    if ( dst->get_post_size() != src->get_post_size()) {
        std::cout<<"Warning: Unequal pre and post neuron numbers..."<<std::endl;
    }
    else {
        for (NeuronID i=0; i<dst->get_post_size(); ++i){
            if(!push_back(i, i, weight_value)){
                std::cout<<"Weight allocation failed!"<<std::endl;
            }
        }
    }
}


void BiasConnection::verify_identity_connection(){
    for (NeuronID i=1; i<=dst->get_post_size(); ++i){
        std::cout<<get(i,i)<<std::endl;
    }
}
