//
//  TransmitEventConnection.cpp
//  
//
//  Created by Alexandre Payeur on 3/5/19.
//

#include "TransmitEventConnection.h"

using namespace auryn;

TransmitEventConnection::TransmitEventConnection(SpikingGroup * source,
                                                 NeuronGroup * destination,
                                                 AurynWeight weight,
                                                 AurynDouble sparseness,
                                                 TransmitterType transmitter,
                                                 std::string name):
TransmitBurstConnection(source, destination, weight, sparseness, transmitter, name) {
    if ( !name.empty() )
        // name will be unempty if no name is given to connection at creation because parent class constructor is called first
        set_name("TransmitEventConnection");
}

void TransmitEventConnection::propagate()
{
    // loop over all spikes to separate bursts from isolated spikes
    SpikeContainer::const_iterator spk;
    for ( spk = src->get_spikes_immediate()->begin() ;
         spk < src->get_spikes_immediate()->end() ;
         ++spk ) {
        
        // detect first spike in bursts and non-burst spikes
        if ( tr_pre->get(*spk) < burst_thr ){
            burst_state->set(*spk,0);
            propagate_forward(*spk);
        }
        
        // detect second spike in burst
        if ( tr_pre->get(*spk) > burst_thr && burst_state->get(*spk) == 0 ){
            burst_state->set(*spk,1);
        }
    }
}
