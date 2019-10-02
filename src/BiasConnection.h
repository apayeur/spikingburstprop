//
//  BiasConnection.h
//  
//
//  Created by Alexandre Payeur on 9/4/19.
//

#ifndef BIASCONNECTION_H_
#define BIASCONNECTION_H_

#include "auryn_definitions.h"
#include "AurynVector.h"
#include "DuplexConnection.h"
#include "Trace.h"
#include "LinearTrace.h"
#include "SpikeDelay.h"
#include "AdaptiveEBCPConnection.h"

namespace auryn {
    
    /*! \brief Plastic connection that depends on postsynaptic events and bursts.
     *  Differ from AdaptiveEBCPConnection by the fact that the learning rule is
     *  independent of presynaptic events.
     *
     */
    
    class BiasConnection : public AdaptiveEBCPConnection
    {
        
    protected:
        void propagate_backward(const NeuronID translated_post, const AurynState valence);
        void propagate_forward();
        void compute_burst_rate();
        
    public:
        BiasConnection(SpikingGroup * source,
                       NeuronGroup * destination,
                       TransmitterType transmitter=GLUT);
        
        BiasConnection(SpikingGroup * source,
                       NeuronGroup * destination,
                       TransmitterType transmitter,
                       AurynFloat eta,
                       AurynFloat maxweight);
        
        BiasConnection(SpikingGroup * source,
                               NeuronGroup * destination,
                               AurynWeight weight,
                               AurynFloat sparseness=0.05,
                               AurynFloat eta=1e-3,
                               AurynFloat maxweight=1.,
                               TransmitterType transmitter=GLUT,
                               string name = "BiasConnection" );
        
        virtual ~BiasConnection();
        virtual void finalize();
        void free();
        virtual void propagate();
        //virtual void evolve();
        
        void construct_identity_connection(AurynWeight weight_value);
        void verify_identity_connection();

    };

}

#endif /* BIASCONNECTION_H_ */
