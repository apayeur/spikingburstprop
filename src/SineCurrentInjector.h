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

#ifndef SineCurrentInjector_hpp
#define SineCurrentInjector_hpp

#include "auryn_definitions.h"
#include "AurynVector.h"
#include "System.h"
#include "Logger.h"
#include "Device.h"
#include "NeuronGroup.h"
#include "CurrentInjector.h"
#include <cmath>


namespace auryn {
    
    /*! \brief Stimulator class to inject a sinusoidal current to a list of neurons within a
     group or all neurons within that group: A*sin(2*pi*f*(t-t_0))
     *
     * Used to inject sinusoidal "currents" to arbitraty neuronal states. Maintains an internal vector with
     * numbers which are added (times auryn_timestep) in each timestep to the neuronal target vector
     * (per default that is the membrane voltage and hence the operation corresponds to injecting a current).
     * Note that because of this current units of SineCurrentInjector are in a sense arbitrary because they depend
     * on the neuron model.
     *
     */
    
    
    
    class SineCurrentInjector : public CurrentInjector
    {
    private:
        AurynFloat A;   // amplitude
        AurynFloat f;   // frequency
        AurynFloat t_0; // time shift
        void free();
        
        enum SineCurrentInjectorMode {
            LIST,
            ALL,
        };
        
    public:        
        /*! \brief Default Constructor
         * @param[target] The target group
         * @param[amplitude] Amplitude of the sine wave (see note on units of current)
         * @param[frequency] Frequency of the sine wave (in Hz)
         * @param[neuron_state_name] The state to manipulate
         */
        SineCurrentInjector(NeuronGroup * target, AurynFloat amplitude=0.01, AurynFloat frequency=1., AurynFloat shift=0., std::string neuron_state_name="mem");
        
        /*! \brief Default Destructor */
        virtual ~SineCurrentInjector();
        
        /*! \brief Mode of operation
         *
         * Determines whether the current should be injected to a list of neurons or all neurons. */
        SineCurrentInjectorMode mode;
        
        /*! \brief The array holding the target neuron ids */
        std::vector<NeuronID> * target_neuron_ids;
        
        /*! Implementation of necessary propagate() function. */
        void execute();
        
        /*! Set amplitude of the sine function. */
        void set_amplitude(AurynFloat ampl);
        
        /*! Set frequency of the sine function. */
        void set_frequency(AurynFloat freq);
        
        /*! Set time shift of the sine function. */
        void set_time_shift(AurynFloat s);
        
    };
    
}

#endif /* SineCurrentInjector_hpp */
