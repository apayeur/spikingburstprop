//
//  TransmitBurstConnection.cpp
//  
//
//  Created by Alexandre Payeur on 3/4/19.
//

#include "TransmitBurstConnection.h"

using namespace auryn;

void TransmitBurstConnection::init()
{
    if ( dst->get_post_size() == 0 ) return;
    
    // burst detector
    const double default_tau = 16e-3;
    tr_pre = src->get_pre_trace(default_tau);
    burst_state = new AurynVector<unsigned short>(src->get_vector_size() );
    burst_state->set_all(0);
    burst_thr = std::exp(-1.0);
}

TransmitBurstConnection::TransmitBurstConnection(
                            SpikingGroup * source,
                            NeuronGroup * destination,
                            AurynWeight weight,
                            AurynDouble sparseness,
                            TransmitterType transmitter,
                            std::string name)
: SparseConnection(source,
                   destination,
                   weight,
                   sparseness,
                   transmitter,
                   name)
{
    init();
    if ( name.empty() )
        set_name("TransmitBurstConnection");

}


void TransmitBurstConnection::free()
{
    delete burst_state;
}


TransmitBurstConnection::~TransmitBurstConnection()
{
    if ( dst->get_post_size() > 0 )
        free();
}

void TransmitBurstConnection::propagate()
{
    // loop over all spikes to separate bursts from isolated spikes
    SpikeContainer::const_iterator spk;
    for ( spk = src->get_spikes_immediate()->begin() ;
         spk < src->get_spikes_immediate()->end() ;
         ++spk ) {
        
        // detect first spike in bursts and non-burst spikes
        if ( tr_pre->get(*spk) < burst_thr ){
            burst_state->set(*spk,0);
        }
        
        // detect second spike in burst
        if ( tr_pre->get(*spk) > burst_thr && burst_state->get(*spk) == 0 ){
            burst_state->set(*spk,1);
            propagate_forward(*spk);
        }
    }
}

void TransmitBurstConnection::propagate_forward(const NeuronID pre)
{
    // loop over all postsynaptic partners
    for (const NeuronID * c = w->get_row_begin(pre) ;
         c != w->get_row_end(pre) ;
         ++c ) { // c = post index
        
        // transmit signal to target at postsynaptic neuron
        AurynWeight * weight = w->get_data_ptr(c);
        transmit( *c , *weight );
    }
}

