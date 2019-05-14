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

#ifndef BCPCONNECTION_H_
#define BCPCONNECTION_H_

#include "auryn_definitions.h"
#include "AurynVector.h"
#include "DuplexConnection.h"
#include "Trace.h"
#include "LinearTrace.h"
#include "SpikeDelay.h"

//#include "auryn.h"


namespace auryn {


    /*! \brief BCP connection is an event-based synaptic plasticity rule which updates weights only upon the detection of a postsynaptic event and burst.
     *
     */
class BCPConnection : public DuplexConnection
{

private:
	void init(AurynFloat eta, AurynFloat maxweight, AurynFloat tau_pre);

protected:
	Trace * tr_pre;
	Trace * tr_post;

	Trace * tr_event;
	Trace * tr_hom;
	Trace * tr_burst;
	AurynVector<unsigned short> * burst_state;

	AurynState burst_thr;

	void propagate_forward();
	void propagate_backward(const NeuronID translated_post, const AurynState valence);
	void compute_burst_rate();

	AurynWeight on_pre(NeuronID post);
	AurynWeight on_post(NeuronID pre);
    
public:
	AurynFloat eta_; /*!< learning rate */

	// homeostasis constants
	double min_rate;
	double max_rate;


	AurynFloat Aplus, Aminus;

	bool stdp_active;

	BCPConnection(SpikingGroup * source, NeuronGroup * destination, 
			TransmitterType transmitter=GLUT);

	BCPConnection(SpikingGroup * source, NeuronGroup * destination, 
			AurynWeight weight, AurynFloat sparseness=0.05, 
			AurynFloat eta=1e-3, 
			AurynFloat maxweight=1.,
            AurynFloat tau_pre=20.e-3,
			TransmitterType transmitter=GLUT,
			string name = "BCPConnection" );


	virtual ~BCPConnection();
	virtual void finalize();
	void free();

	virtual void propagate();
	virtual void evolve();
    
    Trace * get_tr_event();
    Trace * get_tr_burst();

};

}

#endif /*BCPCONNECTION_H_*/
