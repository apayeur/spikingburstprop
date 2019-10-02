//
//  BiasIdentityConnection.hpp
//  
//
//  Created by Alexandre Payeur on 9/9/19.
//

#ifndef BIASIDENTITYCONNECTION_H_
#define BIASIDENTITYCONNECTION_H_

#include "auryn_definitions.h"
#include "AurynVector.h"
#include "IdentityConnection.h"
#include "Trace.h"
#include "LinearTrace.h"
#include "SpikeDelay.h"


namespace auryn {
    /*! \brief Plastic connection that depends on postsynaptic events and bursts.
     *
     *  Differ from AdaptiveEBCPConnection by the fact that the learning rule is
     *  independent of presynaptic events. Contrary to BiasConnection, here the
     *  connection matrix is the identity matrix.
     *
     */
    
    class BiasIdentityConnection : public IdentityConnection
    {
    private:
        void init(AurynFloat eta, AurynWeight w_init);
        
    protected:
        Trace * tr_post;    /*!< Trace that accumulates over postsynaptic spikes */
        Trace * tr_event;   /*!< Trace that accumulates over postsynaptic events */
        Trace * tr_burst;   /*!< Trace that accumulates over postsynaptic bursts */
        AurynVector<unsigned short> * burst_state;
        AurynVectorFloat * weights;
        AurynState burst_thr;
        AurynWeight wmin;
        AurynWeight wmax;
        NeuronID every;  // TODO: make every and offset protected in IdentityConnection
        NeuronID offset;
        
        void propagate_forward();
        void propagate_backward(const NeuronID translated_post, const AurynState valence);
        void compute_burst_rate();
        void free();
        AurynWeight on_pre(NeuronID post);
        
    public:
        AurynFloat eta_;    /*!< learning rate */
        double min_rate;    // minimum event rate (homeostatic constant)
        double max_rate;    // minimum event rate (homeostatic constant)
        bool stdp_active;   // whether plasticity is active or not
        
        BiasIdentityConnection(SpikingGroup * source,
                               NeuronGroup * destination,
                               AurynWeight weight = 1.0,
                               AurynFloat eta=1e-3,
                               TransmitterType transmitter=GLUT,
                               string name = "BiasIdentityConnection");
        
        virtual ~BiasIdentityConnection();
        virtual void finalize();
        
        virtual void propagate();
        virtual void evolve();
        
        virtual void set_post_trace_event_tau(AurynFloat x);
        virtual void set_post_trace_burst_tau(AurynFloat x);
        
        virtual void set_min_weight(AurynWeight value);
        virtual void set_max_weight(AurynWeight value);
        virtual AurynWeight get_min_weight();
        virtual AurynWeight get_max_weight();
        
        virtual void stats (AurynDouble &mean, AurynDouble &std, StateID zid=0);
        
        Trace * get_tr_event();
        Trace * get_tr_burst();

    };
    
}


#endif /* BIASIDENTITYCONNECTION_H_ */
