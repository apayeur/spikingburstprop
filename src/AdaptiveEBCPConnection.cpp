#include "AdaptiveEBCPConnection.h"

using namespace auryn;

AdaptiveEBCPConnection::AdaptiveEBCPConnection(
                             SpikingGroup * source,
                             NeuronGroup * destination,
                             AurynWeight weight,
                             AurynFloat sparseness,
                             AurynFloat eta,
                             AurynFloat maxweight,
                             AurynFloat tau_pre,
                             TransmitterType transmitter,
                               std::string name)
: EBCPConnection(source, destination, weight, sparseness, eta, maxweight, tau_pre, transmitter, name)
{
    if ( !name.empty() )
        set_name("AdaptiveEBCPConnection");
}

void AdaptiveEBCPConnection::finalize() {
    EBCPConnection::finalize();
}

void AdaptiveEBCPConnection::compute_burst_rate()
{
    // loop over all spikes to separate events from bursts
    SpikeContainer::const_iterator spk;
    for ( spk = dst->get_spikes_immediate()->begin() ;
         spk < dst->get_spikes_immediate()->end() ;
         ++spk ) {
        const NeuronID s = dst->global2rank(*spk);
        
        // detect first spike in bursts and non-burst spikes
        if ( tr_post->get(s) < burst_thr ) {
            tr_event->inc(s);
            burst_state->set(s,0);
            propagate_backward(s, -tr_burst->normalized_get(s));
        }
        
        // detect second spike in burst
        if ( tr_post->get(s) > burst_thr && burst_state->get(s) == 0 ) {
            tr_burst->inc(s);
            burst_state->set(s,1);
            propagate_backward(s, tr_event->normalized_get(s));
            //propagate_backward(s, tr_event->normalized_get(s)*exp(10e-3/tr_event->get_tau()) - 1./tr_event->get_tau()); //DEBUG
        }
    }
}

void AdaptiveEBCPConnection::propagate()
{
    compute_presyn_event_rate();
    compute_burst_rate();
}

void AdaptiveEBCPConnection::set_post_trace_event_tau(AurynFloat x){
    tr_event->set_timeconstant(x);
}

void AdaptiveEBCPConnection::set_post_trace_burst_tau(AurynFloat x){
    tr_burst->set_timeconstant(x);
}

