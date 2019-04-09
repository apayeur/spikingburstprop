//
//  TransmitBurstConnection.h
//  
//
//  Created by Alexandre Payeur on 3/4/19.
//

#ifndef TRANSMITBURSTCONNECTION_H_
#define TRANSMITBURSTCONNECTION_H_

#include "auryn_definitions.h"
#include "AurynVector.h"
#include "SparseConnection.h"
#include "Trace.h"
#include "LinearTrace.h"


namespace auryn {

    /*! \class TransmitBurstConnection
     *  \brief Connection that transmit presynaptic bursts.
     *
     *  This synaptic connection extracts bursts from the presynaptic
     *  spike trains and propagates them to the postsynaptic targets.
     *  Isolated spikes (i.e., not part of a burst)
     *  are ignored and not propagated. 
     *
     *  \see sim_transmitbursts.cpp
     */
    
    class TransmitBurstConnection : public SparseConnection {
    
    private:
        void init();
        
    protected:
        Trace * tr_pre;
        AurynVector<unsigned short> * burst_state;
        AurynState burst_thr;
        
        /*! \brief Called by the propagate function whenever a
         *  presynaptic burst is detected.
         */
        void propagate_forward(const NeuronID pre);
        
    public:
        
        /*! \brief Constructor that should be used with this class.*/
        TransmitBurstConnection(SpikingGroup * source,
                                NeuronGroup * destination,
                                AurynWeight weight,
                                AurynDouble sparseness=0.05,
                                TransmitterType transmitter=GLUT,
                                string name="TransmitBurstConnection");
        
        /*! \brief The default destructor */
        virtual ~TransmitBurstConnection();
        void free();
        
        /*! \brief Internally used propagate method
         *
         *  This method propagates spikes in the main simulation loop. Should usually not be
         *  called directly by the user.
         */
        virtual void propagate();
        
    };
    
}

#endif /* TRANSMITBURSTCONNECTION_H_ */
