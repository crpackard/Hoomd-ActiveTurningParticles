
#include "IntegrationMethodTwoStep.h"
#include "NeighborList.h"

#ifndef __TWO_STEP_VICSEK_BASE__
#define __TWO_STEP_VICSEK_BASE__


class PYBIND11_EXPORT TwoStepVicsekBase : public IntegrationMethodTwoStep
    {
    public:
        //! Constructs the integration method and associates it with the system
        TwoStepVicsekBase(std::shared_ptr<SystemDefinition> sysdef,
                           std::shared_ptr<ParticleGroup> group,
                           std::shared_ptr<NeighborList> nlist,
                           Scalar alpha,
                           Scalar beta,
                           Scalar v0,
                           Scalar delta,
                           Scalar bias,
                           Scalar tau,
                           Scalar sym,
                           unsigned int seed);
        virtual ~TwoStepVicsekBase();

    protected:
        std::shared_ptr<ParticleGroup> m_group;
        std::shared_ptr<NeighborList> m_nlist;
        Scalar m_alpha;
        Scalar m_beta;
        Scalar m_v0;
        Scalar m_delta;
        Scalar m_bias;
        Scalar m_tau;
        Scalar m_sym;
        unsigned int m_seed;
        Scalar m_deltaT;
    };

//! Exports the TwoStepLangevinBase class to python
void export_TwoStepVicsekBase(pybind11::module& m);

#endif // #ifndef __TWO_STEP_LANGEVIN_BASE__
