//
//  EventBurstPoissonGroup.hpp
//  
//
//  Created by Alexandre Payeur on 6/3/19.
//

#ifndef EVENTBURSTPOISSONGROUP_H_
#define EVENTBURSTPOISSONGROUP_H_

//#include "auryn.h"
#include "auryn_definitions.h"
#include "AurynVector.h"
#include "BurstPoissonGroup.h"
#include "System.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>


namespace auryn {
    
    /*! \brief This file implements a simplified two-compartment bursting neuron model.
     *
     *  The model mostly follows the dynamics of the
     *  NaudGroup, except that the bursts and events
     *  are stochastically generated.
     *
     * */
    
    class EventBurstPoissonGroup : public BurstPoissonGroup
    {
    protected:
        // auxiliary state vector for burst and somatic spike probability
        AurynStateVector * t_prob_soma;
        AurynStateVector * t_prob_dend;
        
    public:
        /*! \brief Default constructor.
         *
         * @param size the size of the group.
         * @param distmode Node distribution mode
         */
        EventBurstPoissonGroup( NeuronID size, NodeDistributionMode distmode=AUTO );
        virtual ~EventBurstPoissonGroup();
        
        void integrate_membrane();
        
        void check_thresholds();
        
        /*! Internally used evolve function. Called by System. */
        virtual void evolve();
        
        
    };
}

#endif /* EventBurstPoissonGroup_H_ */
