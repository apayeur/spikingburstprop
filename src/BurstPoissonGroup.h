//
//  BurstPoissonGroup.hpp
//  
//
//  Created by Alexandre Payeur on 4/7/19.
//

#ifndef BURSTPOISSONGROUP_H_
#define BURSTPOISSONGROUP_H_

//#include "auryn.h"
#include "auryn_definitions.h"
#include "AurynVector.h"
#include "NaudGroup.h"
#include "System.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>


namespace auryn {
    
    /*! \brief This file implements a simplified two-compartment bursting neuron model.
     *
     *  The model mostly follows the dynamics of the
     *  NaudGroup, except for the bursts that
     *  are stochastically generated.
     *
     * */
    
    class BurstPoissonGroup : public NaudGroup
    {
    protected:
        static boost::mt19937 gen;
        boost::random::uniform_real_distribution<> * dist;
        boost::variate_generator<boost::mt19937&, boost::random::uniform_real_distribution<> > * die;
        
        unsigned int salt;
        AurynVector<unsigned int> * refractory_state;
        AurynVector<unsigned int> * burst_state;

        unsigned int abs_ref_period; //!< absolute refractory period following somatic spike
        float burst_duration;
        
        
    public:
        /*! \brief Default constructor.
         *
         * @param size the size of the group.
         * @param distmode Node distribution mode
         */
        BurstPoissonGroup( NeuronID size, NodeDistributionMode distmode=AUTO );
        virtual ~BurstPoissonGroup();
        
        void integrate_membrane();
        
        void check_thresholds();
        
        /*! Internally used evolve function. Called by System. */
        virtual void evolve();
        
        /*! Use this to seed the random number generator. */
        void seed(unsigned int s);
        
    };
}

#endif /* BURSTPOISSONGROUP_H_ */
