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

#include "SineCurrentInjector.h"

using namespace auryn;

SineCurrentInjector::SineCurrentInjector(NeuronGroup * target, AurynFloat amplitude, AurynFloat frequency, AurynFloat shift, std::string neuron_state_name) : CurrentInjector(target, neuron_state_name, 0.0 ), A(amplitude), f(frequency), t_0(shift)
{
    target_neuron_ids = new std::vector<NeuronID>;
    mode = ALL;
}


void SineCurrentInjector::free( )
{
    delete target_neuron_ids;
}


SineCurrentInjector::~SineCurrentInjector()
{
    free();
}

void SineCurrentInjector::execute()
{
    if ( dst->evolve_locally() ) {
        AurynState cur = A*sin(2*M_PI*f*(sys->get_time()-t_0));
        switch ( mode ) {
            case LIST:
                currents->set_all(0.0);
                for ( int i = 0 ; i < target_neuron_ids->size() ; ++i ) {
                    currents->set(target_neuron_ids->at(i),cur);
                }
                break;
            case ALL:
            default:
                currents->set_all(cur);
        }
        CurrentInjector::execute();
    }
}

void SineCurrentInjector::set_amplitude(AurynFloat ampl)
{
    A = ampl;
}

void SineCurrentInjector::set_frequency(AurynFloat freq)
{
    f = freq;
}

void SineCurrentInjector::set_time_shift(AurynFloat s)
{
    t_0 = s;
}
