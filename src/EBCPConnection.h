//
//  EBCPConnection.h
//  
//
//  Created by Alexandre Payeur on 2/22/19.
//

#ifndef EBCPCONNECTION_H_
#define EBCPCONNECTION_H_

#include "auryn_definitions.h"
#include "AurynVector.h"
#include "DuplexConnection.h"
#include "Trace.h"
#include "LinearTrace.h"
#include "SpikeDelay.h"
#include "BCPConnection.h"

namespace auryn {
    
    /*! \brief Plastic connection that depends on (filtered) presynaptic events,
     *  and postsynaptic events and bursts.
     *
     *  Implements a learning rule derived from the BCP connection class learning
     *  rule. Here, only presynaptic events (as opposed to whole spike trains) are
     *  taken into account.
     */
    
    class EBCPConnection : public BCPConnection
    {
        protected:
            Trace * tr_event_pre;
            AurynVector<unsigned short> * burst_state_pre;

            void propagate_forward(const NeuronID pre);
            void propagate_backward(const NeuronID translated_post, const AurynState valence);
            void compute_presyn_event_rate();

        public:
        
            EBCPConnection(SpikingGroup * source,
                           NeuronGroup * destination,
                           TransmitterType transmitter=GLUT);
        
            EBCPConnection(SpikingGroup * source,
                           NeuronGroup * destination,
                           TransmitterType transmitter,
                           AurynFloat eta,
                           AurynFloat maxweight,
                           AurynFloat tau_pre);
        
            EBCPConnection(SpikingGroup * source,
                           NeuronGroup * destination,
                           AurynWeight weight,
                           AurynFloat sparseness=0.05,
                           AurynFloat eta=1e-3,
                           AurynFloat maxweight=1.,
                           AurynFloat tau_pre=20.e-3,
                           TransmitterType transmitter=GLUT,
                           string name = "EBCPConnection" );
        
            virtual ~EBCPConnection();
            virtual void finalize();
            void free();
        
            virtual void propagate();
            virtual void evolve();

    };
    
}

#endif /*EBCPCONNECTION_H_*/
