/*
 * Copyright 2014-2017 Friedemann Zenke
 *
 * This file is part of Auryn, a simulation package for plastic
 * spiking neural networks.
 *
 * Auryn is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Auryn is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Auryn.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "EBCPConnection.h"

using namespace auryn;

EBCPConnection::EBCPConnection(SpikingGroup * source,
                               NeuronGroup * destination,
                               AurynWeight weight,
                               AurynFloat sparseness,
                               AurynFloat eta,
                               AurynFloat maxweight,
                               AurynFloat tau_pre,
                               TransmitterType transmitter,
                               std::string name)
: BCPConnection(source, destination, weight, sparseness, eta, maxweight, tau_pre, transmitter, name)
{
    burst_state_pre = new AurynVector<unsigned short>(src->get_vector_size() );
    burst_state_pre->set_all(0);

    tr_event_pre = new EulerTrace(src->get_vector_size(), tau_pre);
    
    if ( !name.empty() )
        set_name("EBCPConnection");
}

EBCPConnection::EBCPConnection(SpikingGroup * source,
                               NeuronGroup * destination,
                               TransmitterType transmitter,
                               AurynFloat eta,
                               AurynFloat maxweight,
                               AurynFloat tau_pre)
: BCPConnection(source, destination, transmitter, eta, maxweight, tau_pre)
{
    burst_state_pre = new AurynVector<unsigned short>(src->get_vector_size() );
    burst_state_pre->set_all(0);
    
    tr_event_pre = new EulerTrace(src->get_vector_size(), tau_pre);
}

EBCPConnection::EBCPConnection(SpikingGroup * source,
                               NeuronGroup * destination,
                               TransmitterType transmitter)
: BCPConnection(source, destination, transmitter)
{
}


void EBCPConnection::finalize() {
    BCPConnection::finalize();
}


void EBCPConnection::free()
{
    //delete burst_state;
    delete burst_state_pre;
    //delete tr_event;
    //delete tr_burst;
    delete tr_event_pre;
}

EBCPConnection::~EBCPConnection()
{
    if ( dst->get_post_size() > 0 )
        free();
}


void EBCPConnection::compute_presyn_event_rate()
{
    // loop over all presynaptic spikes to separate events from bursts
    SpikeContainer::const_iterator spk;
    for ( spk = src->get_spikes_immediate()->begin() ;
         spk < src->get_spikes_immediate()->end() ;
         ++spk ) {
        
        // detect first spike in bursts and non-burst spikes
        if ( tr_pre->get(*spk) < burst_thr ) {
            burst_state_pre->set(*spk,0);
            propagate_forward(*spk);
            tr_event_pre->inc(*spk);
        }
    }
}

void EBCPConnection::propagate_forward(const NeuronID pre)
{
    // loop over all postsynaptic partners
    for (const NeuronID * c = w->get_row_begin(pre) ;
         c != w->get_row_end(pre) ;
         ++c ) { // c = post index
        
        // transmit signal to target at postsynaptic neuron
        AurynWeight * weight = w->get_data_ptr(c);
        transmit( *c , *weight );
        
        // handle plasticity
        if ( stdp_active ) {
            // performs weight update
            *weight += on_pre(*c);
            
            // clips weights
            if ( *weight > get_max_weight() ) *weight = get_max_weight();
            else
                if ( *weight < get_min_weight() ) *weight = get_min_weight();
        }
    }
}

void EBCPConnection::propagate_backward(const NeuronID translated_post, const AurynState valence)
{
    if (stdp_active) {
        // compute hom part of local error signal
        double hom_neg = 0.0;
        const double rate = tr_hom->normalized_get(translated_post);
        if (rate>max_rate) {
            // std::cout << "hom dep" << std::endl;
            hom_neg = -(rate-max_rate);
        }
        const NeuronID post = dst->rank2global(translated_post);
        
        // loop over all presynaptic partners
        for (const NeuronID * c = bkw->get_row_begin(post) ; c != bkw->get_row_end(post) ; ++c ) {
            // computes plasticity update
            AurynWeight * weight = bkw->get_data(c);
            
            const NeuronID pre = *c;
            *weight += eta_*((valence+hom_neg)*tr_event_pre->get(pre));
            
            // clips weights
            if ( *weight > get_max_weight() ) *weight = get_max_weight();
            else
                if ( *weight < get_min_weight() ) *weight = get_min_weight();
        }
    }
}


void EBCPConnection::propagate()
{
    compute_presyn_event_rate();
    compute_burst_rate();
}


void EBCPConnection::evolve()
{
    tr_event->evolve();
    tr_burst->evolve();
    tr_event_pre->evolve();
}
