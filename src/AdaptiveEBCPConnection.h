//
//  AdaptiveEBCPConnection.h
//  
//
//  Created by Alexandre Payeur on 5/13/19.
//

#ifndef AdaptiveEBCPCONNECTION_H_
#define AdaptiveEBCPCONNECTION_H_

#include "auryn_definitions.h"
#include "AurynVector.h"
#include "DuplexConnection.h"
#include "Trace.h"
#include "LinearTrace.h"
#include "SpikeDelay.h"
#include "EBCPConnection.h"

namespace auryn {
    
    /*! \brief Plastic connection that depends on presynaptic events,
     *  and postsynaptic events and bursts. Differ from EBCPConnection by
     *  the fact that A_plus and A_minus are activity-dependent.
     *
     *  Implements a learning rule derived from the BCP connection class learning
     *  rule. Here, only presynaptic events (as opposed to whole spike trains) are
     *  taken into account.
     */
    
    class AdaptiveEBCPConnection : public EBCPConnection
    {
        protected:
            void compute_burst_rate();

        public:
            AdaptiveEBCPConnection(SpikingGroup * source, NeuronGroup * destination,
                          AurynWeight weight, AurynFloat sparseness=0.05,
                          AurynFloat eta=1e-3,
                          AurynFloat maxweight=1.,
                          AurynFloat tau_pre=20.e-3,
                          TransmitterType transmitter=GLUT,
                          string name = "AdaptiveEBCPConnection" );
        
            //virtual ~AdaptiveEBCPConnection();
            virtual void finalize();
            //void free();
            virtual void propagate();
            //virtual void evolve();

    };
    
}

#endif /*AdaptiveEBCPCONNECTION_H_*/
