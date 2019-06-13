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
    // homeostasis constants (min and max event rate)
    min_rate = 0.1;
    max_rate = 10.0;
}

void AdaptiveEBCPConnection::finalize() {
    EBCPConnection::finalize();
}

AurynWeight AdaptiveEBCPConnection::on_pre(NeuronID post)
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

void AdaptiveEBCPConnection::propagate_backward(const NeuronID translated_post, const AurynState valence)
{
    if (stdp_active) {
        // compute hom part of local error signal
        double hom_neg = 0.0;
        const double rate = tr_event->normalized_get(translated_post);
        if (rate>max_rate) {
            // std::cout << "hom dep" << std::endl;
            hom_neg = -(rate - max_rate);
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

