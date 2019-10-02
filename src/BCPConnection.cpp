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

#include "BCPConnection.h"

using namespace auryn;

void BCPConnection::init(AurynFloat eta, AurynFloat maxweight, AurynFloat tau_pre)
{
	if ( dst->get_post_size() == 0 ) return;

	eta_ = eta; // learning rate

	Aplus = 1.0; // TODO set to something more meaningful
	Aminus = -0.15;

	// homeostasis constants
	min_rate = 0.1;
	max_rate = 20.0;

	auryn::logger->parameter("eta",eta);

	// presynaptic trace
	tr_pre  = src->get_pre_trace(tau_pre);

	// burst detector
	const double default_tau = 16e-3;
	tr_post = dst->get_post_trace(default_tau);
	burst_state = new AurynVector<unsigned short>( dst->get_vector_size() );
	burst_state->set_all(0);
	burst_thr = std::exp(-1.0);

	// other traces
	tr_event = new EulerTrace(dst->get_vector_size(), 5.5); //!< event rate average
	tr_burst = new EulerTrace(dst->get_vector_size(), 5.5); //!< burst rate average (currently not needed)
	tr_hom = dst->get_post_trace(30.0); //!< homeostatic firing rate average
	tr_hom->set_all(30.0*min_rate);

	// min/max weight
	set_min_weight(0.0);
	set_max_weight(maxweight);

	stdp_active = true;
}


void BCPConnection::finalize() {
	DuplexConnection::finalize();
}

void BCPConnection::free()
{
	delete burst_state;
	delete tr_event;
	delete tr_burst;
}

BCPConnection::BCPConnection(SpikingGroup * source,
                             NeuronGroup * destination,
                             TransmitterType transmitter)
: DuplexConnection(source,
                   destination,
                   transmitter)
{
}

BCPConnection::BCPConnection(SpikingGroup * source,
                             NeuronGroup * destination,
                             TransmitterType transmitter,
                             AurynFloat eta,
                             AurynFloat maxweight,
                             AurynFloat tau_pre)
: DuplexConnection(source,
                   destination,
                   transmitter)
{
    init(eta, maxweight, tau_pre);
}


BCPConnection::BCPConnection(SpikingGroup * source,
                             NeuronGroup * destination,
                             AurynWeight weight,
                             AurynFloat sparseness,
                             AurynFloat eta,
                             AurynFloat maxweight,
                             AurynFloat tau_pre,
                             TransmitterType transmitter,
                             std::string name)
: DuplexConnection(source, 
                   destination,
                   weight,
                   sparseness,
                   transmitter,
                   name)
{
	init(eta, maxweight, tau_pre);
	if ( name.empty() )
		set_name("BCPConnection");
}

BCPConnection::~BCPConnection()
{
	if ( dst->get_post_size() > 0 ) 
		free();
}


AurynWeight BCPConnection::on_pre(NeuronID post)
{
	const NeuronID ts = dst->global2rank(post); // only to be used for post traces
	// compute hom part of local error signal
	double hom = 0.0;
	const double rate = tr_hom->normalized_get(ts);
	if (rate<min_rate) {
		hom = (min_rate-rate);
		// std::cout << "hom event" << std::endl;
	} 
	AurynDouble dw = eta_*hom;
	return dw;
}

AurynWeight BCPConnection::on_post(NeuronID pre)
{
	AurynDouble dw = tr_pre->get(pre);
	return dw;
}


void BCPConnection::propagate_forward()
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
}

void BCPConnection::propagate_backward(const NeuronID translated_post, const AurynState valence)
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
			*weight += eta_*((valence+hom_neg)*tr_pre->get(pre));

			// clips weights
			if ( *weight > get_max_weight() ) *weight = get_max_weight(); 
			else
			if ( *weight < get_min_weight() ) *weight = get_min_weight();
		}
	}
}

void BCPConnection::compute_burst_rate()
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
			propagate_backward(s,Aminus);
		}

		// detect second spike in burst
		if ( tr_post->get(s) > burst_thr && burst_state->get(s) == 0 ) { 
			tr_burst->inc(s);
			burst_state->set(s,1);
			propagate_backward(s,Aplus);
		}
	}
}

void BCPConnection::propagate()
{
	propagate_forward();
	compute_burst_rate();
}

void BCPConnection::evolve()
{
	tr_event->evolve();
	tr_burst->evolve();
}

Trace * BCPConnection::get_tr_event()
{
    return tr_event;
}

Trace * BCPConnection::get_tr_burst()
{
    return tr_burst;
}
