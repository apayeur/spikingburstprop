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

#ifndef STPeTMConnection_h
#define STPeTMConnection_h

#include "auryn_definitions.h"
#include "AurynVector.h"
#include "SparseConnection.h"
#include "STPConnection.h"


namespace auryn {
    
    
    
    /*! \brief This class implements short term plasticity according to the extended Tsodyks-Markram synapse. It is a simple modification of the class STPConnection.
     *
     * This class implements the short-term plasticity model following
     * Costa, Sjostrom, Van Rossum. Front. Comp. Neurosci. 2013
     *
     */
    
    class STPeTMConnection : public STPConnection
    {
    private:
        double f;
        
    public:
        
        /*! Minimal constructor for from file init -- deprecated
         * \param filename Filename to load
         */
        STPeTMConnection(const char * filename);
        
        /*! Minimal constructor to which leaves the connection uninitialized
         * \param source The presynaptic SpikingGroup
         * \param destination the postsynaptic NeuronGroup
         */
        STPeTMConnection(NeuronID rows, NeuronID cols);
        
        /*! Default constructor to which leaves the connection uninitialized
         * \param source The presynaptic SpikingGroup
         * \param destination the postsynaptic NeuronGroup
         * \param transmitter The transmitter type
         */
        STPeTMConnection(SpikingGroup * source, NeuronGroup * destination, TransmitterType transmitter=GLUT);
        
        /*! Default constructor to initialize connection from file
         * \param source The presynaptic SpikingGroup
         * \param destination the postsynaptic NeuronGroup
         * \param filename The file to load the connectivity from upon initialization
         * \param transmitter The transmitter type
         * \param name The connection name as it appears in debugging output
         */
        STPeTMConnection(SpikingGroup * source, NeuronGroup * destination, const char * filename , TransmitterType transmitter=GLUT);
        
        /*! Default constructor to initialize connection with random sparse connectivity
         * \param source The presynaptic SpikingGroup
         * \param destination the postsynaptic NeuronGroup
         * \param weight The default weight for connections
         * \param sparseness The probability of a connection for the sparse connectivity
         * \param transmitter The transmitter type
         * \param name The connection name as it appears in debugging output
         */
        STPeTMConnection(SpikingGroup * source, NeuronGroup * destination, AurynWeight weight, AurynFloat sparseness=0.05, TransmitterType transmitter=GLUT, string name="STPeTMConnection");
        
        /*! Setter for f, the multiplicative factor of synaptic facilitation. */
        void set_f(AurynFloat f_);
        
        /*! \brief Clone parameters  */
        void clone_parameters(STPeTMConnection * con);
                
        /*! Internal function to push spike attributes. */
        void push_attributes();
        
    };
    
}

#endif /* STPeTMConnection_h */
