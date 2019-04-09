//
//  TransmitEventConnection.hpp
//  
//
//  Created by Alexandre Payeur on 3/5/19.
//

#ifndef TRANSMITEVENTCONNECTION_H_
#define TRANSMITEVENTCONNECTION_H_

#include "auryn_definitions.h"
#include "AurynVector.h"
#include "SparseConnection.h"
#include "Trace.h"
#include "LinearTrace.h"
#include "TransmitBurstConnection.h"

namespace auryn {
    
    /*! \class TransmitEventConnection
     *  \brief Connection that transmit presynaptic events.
     *
     *  It can be viewed as an idealized short-term depressing (STD) synapse.
     */
    class TransmitEventConnection : public TransmitBurstConnection {
    public:
        /*! \brief Constructor that must be used.
         *
         *  Merely calls the parent class' constructor.
         */
         TransmitEventConnection(SpikingGroup * source, NeuronGroup * destination, AurynWeight weight, AurynDouble sparseness=0.05, TransmitterType transmitter=GLUT, string name="TransmitEventConnection");
        
        /*! \brief Internally used propagate method
         *
         * This method propagates spikes in the main simulation loop. Should usually not be called directly by the user.*/
        virtual void propagate();
    };
}
        
#endif /* TRANSMITEVENTCONNECTION_H_ */
