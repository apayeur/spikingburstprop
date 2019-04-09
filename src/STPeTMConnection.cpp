/*
 * Copyright 2014-2018 Friedemann Zenke
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
 *
 * If you are using Auryn or parts of it for your work please cite:
 * Zenke, F. and Gerstner, W., 2014. Limits to high-speed simulations
 * of spiking neural networks using general-purpose computers.
 * Front Neuroinform 8, 76. doi: 10.3389/fninf.2014.00076
 */

#include "STPeTMConnection.h"

using namespace auryn;


STPeTMConnection::STPeTMConnection(const char * filename) : STPConnection(filename)
{
    if ( src->get_rank_size() > 0 ) {
        f = 0.1; // hard-coded at construction time, but can be changed afterwards with setter
    }
}


STPeTMConnection::STPeTMConnection(NeuronID rows, NeuronID cols) : STPConnection(rows, cols)
{
    if ( src->get_rank_size() > 0 ) {
        f = 0.1;
    }
}


STPeTMConnection::STPeTMConnection(SpikingGroup * source, NeuronGroup * destination, TransmitterType transmitter) : STPConnection(source, destination, transmitter)
{
    if ( src->get_rank_size() > 0 ) {
        f = 0.1;
    }
}


STPeTMConnection::STPeTMConnection(SpikingGroup * source, NeuronGroup * destination, const char * filename , TransmitterType transmitter) : STPConnection(source, destination, filename, transmitter)
{
    if ( src->get_rank_size() > 0 ) {
        f = 0.1;
    }
}

STPeTMConnection::STPeTMConnection(SpikingGroup * source, NeuronGroup * destination, AurynWeight weight, AurynFloat sparseness, TransmitterType transmitter, string name) : STPConnection(source, destination, weight, sparseness, transmitter, name)
{
    if ( src->get_rank_size() > 0 ) {
        f = 0.1;
    }
}


void STPeTMConnection::push_attributes()
{
    // need to push one attribute for each spike
    SpikeContainer * spikes = src->get_spikes_immediate();
    for (SpikeContainer::const_iterator spike = spikes->begin() ;
         spike != spikes->end() ; ++spike ) {
        // dynamics
        NeuronID spk = src->global2rank(*spike);
        double x = auryn_vector_float_get( state_x, spk );
        double u = auryn_vector_float_get( state_u, spk );
        auryn_vector_float_set( state_x, spk, x-u*x );
        auryn_vector_float_set( state_u, spk, u+f*(1-u) );
        
        // TODO spike translation or introduce local_spikes
        // function in SpikingGroup and implement this there ... (better option)
        src->push_attribute( x*u );
        
    }
    
    // If we had two spike attributes in this connection we push
    // the second attribute for each spike here:
    //
    // SpikeContainer * spikes = src->get_spikes_immediate();
    // for (SpikeContainer::const_iterator spike = spikes->begin() ;
    //         spike != spikes->end() ; ++spike ) {
    //     AurynFloat other_attribute = foo+bar;
    //     src->push_attribute( other_attribute );
    // }
}


void STPeTMConnection::set_f(AurynFloat f_) {
    f = f_;
}


void STPeTMConnection::clone_parameters(STPeTMConnection * con)
{
    tau_d = con->tau_d;
    tau_f = con->tau_f;
    Ujump = con->Ujump;
    f = con->f;
    clear();
}
