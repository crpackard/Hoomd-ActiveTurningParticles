
#include <hoomd/extern/pybind/include/pybind11/pybind11.h>
#include "hoomd/extern/pybind/include/pybind11/numpy.h"

#include <math.h>
#include <random>

#include "hoomd/RandomNumbers.h"
#include "hoomd/RNGIdentifiers.h"

#include "hoomd/VectorMath.h"
#include "hoomd/HOOMDMath.h"

#include "TwoStepVicsek.h"

#include <iostream>
namespace py = pybind11;
using namespace std;
using namespace hoomd;


TwoStepVicsek::TwoStepVicsek(
        std::shared_ptr<SystemDefinition> sysdef,
        std::shared_ptr<ParticleGroup> group,
        std::shared_ptr<NeighborList> nlist,
        Scalar alpha,
        Scalar beta,
        Scalar v0,
        Scalar delta,
        Scalar bias,
        Scalar tau,
        Scalar sym,
        unsigned int seed)
        : TwoStepVicsekBase(sysdef, group, nlist, alpha, beta, v0, delta, bias, tau, sym, seed),
          m_group(group), 
          m_nlist(nlist), 
          m_alpha(alpha),
          m_beta(beta),
          m_v0(v0), 
          m_delta(delta),
          m_bias(bias),
          m_tau(tau),
          m_sym(sym),
          m_seed(seed)
        {
        m_exec_conf->msg->notice(5) << "Constructing IntegrateVicsek" << endl;
        }

TwoStepVicsek::~TwoStepVicsek()
    {
    m_exec_conf->msg->notice(5) << "Destroying IntegrateVicsek" << endl;
    }


std::vector< std::string > TwoStepVicsek::getProvidedLogQuantities()
    {
    vector<string> result;
    return result;
    }

/*! \param quantity Name of the log quantity to get
    \param timestep Current time step of the simulation
    \param my_quantity_flag passed as false, changed to true if quantity logged here
*/
Scalar TwoStepVicsek::getLogValue(const std::string& quantity, unsigned int timestep, bool &my_quantity_flag)
    {
    return Scalar(0);
    }



void TwoStepVicsek::integrateStepOne(unsigned int timestep)
    {
    // start by updating the neighborlist (see md/NeighborList.h for documentation)
    assert(m_nlist);
    m_nlist->compute(timestep);
    ArrayHandle<unsigned int> h_n_neigh(m_nlist->getNNeighArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_nlist(m_nlist->getNListArray(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_head_list(m_nlist->getHeadList(), access_location::host, access_mode::read);
   
    // Retrieve particles' d.o.f. and index tags for colored Vicsek model (see hoomd/ParticleData.cc for documentation)
    ArrayHandle<unsigned int> h_tag(m_pdata->getTags(), access_location::host, access_mode::read);
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_orientation(m_pdata->getOrientationArray(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar3> h_inertia(m_pdata->getMomentsOfInertiaArray(), access_location::host, access_mode::readwrite);

    const BoxDim& box = m_pdata->getGlobalBox();
    unsigned int group_size = m_group->getNumMembers();
    
    // Create a list to hold updated particle orientations and auxiliary variables.
    double updated_orientations[group_size];
    double updated_omegas[group_size];

    // Loop through all particles and update each's angular variables
    for (unsigned int i = 0; i < group_size; i++)
        {

        // Generate a random noise acting on the angular momentum.
        unsigned int ptag = h_tag.data[i];
        RandomGenerator rng(RNGIdentifier::TwoStepLangevin, m_seed, ptag, timestep);
        hoomd::NormalDistribution<Scalar> normal(1, 0);
        unsigned int dt = 1;
        Scalar zeta = normal(rng);

        // Fetch the current angular momentum, orientation, and position of particle `i`.
        vec3<Scalar> I(h_inertia.data[i]);
        double omega_i = atan2(I.y, I.x);
        double theta_i = h_orientation.data[i].x;
        Scalar3 r_i = make_scalar3(h_pos.data[i].x, h_pos.data[i].y, 0);

        // Initialize variables to record total angular momentum and orientational interactions with neighbors.
        double W = 0;
        double J = 0;

        // Initialize variable to record total number of interacting neighbors.
        unsigned int N = 0;

        // Loop over the set of instantaneous neighbors.
        const unsigned int myHead = h_head_list.data[i];
        const unsigned int size = (unsigned int)h_n_neigh.data[i];
        for (unsigned int k = 0; k < size; k++)
            {
            // Fetch the current position of particle `j`.
            unsigned int j = h_nlist.data[myHead + k];
            Scalar3 r_j = make_scalar3(h_pos.data[j].x, h_pos.data[j].y, 0);
            // Evaluate if `j` is a neighbor of `i`.
            Scalar3 dr = r_i - r_j;
            dr = box.minImage(dr);
            double drsq = dot(dr,dr);
            if(drsq < 1 && j != i)
                {
                // Fetch the current angular momentum and orientation of particle `j`.
                vec3<Scalar> I(h_inertia.data[j]);
                double omega_j = atan2(I.y, I.x);
                double theta_j = h_orientation.data[j].x;
                W += slow::sin(omega_j-omega_i);
                J += slow::sin(m_sym*(theta_j-theta_i));
                N += 1;
                }
            }

        if(N==0)
            {
            updated_omegas[i] = omega_i - (dt/m_tau)*(omega_i-m_bias) + m_delta*zeta;
            updated_orientations[i] = theta_i + updated_omegas[i];
            }
        else
            {
            updated_omegas[i] = omega_i + (m_beta/N)*W - (dt/m_tau)*(omega_i-m_bias) + m_delta*zeta;
            updated_orientations[i] = theta_i + (m_alpha/N)*J + updated_omegas[i];
            }
        }

    // Update noise & orientation variable of every particle (instantaneously)
    for (unsigned int i=0; i<group_size; i++)
        {
        h_orientation.data[i].x = updated_orientations[i];
        h_inertia.data[i].x = slow::cos(updated_omegas[i]);
        h_inertia.data[i].y = slow::sin(updated_omegas[i]);
        h_inertia.data[i].z = 0;
        }
    }


void TwoStepVicsek::integrateStepTwo(unsigned int timestep)
    {
    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::readwrite);
    ArrayHandle<Scalar4> h_orientation(m_pdata->getOrientationArray(), access_location::host, access_mode::readwrite);
    ArrayHandle<int3> h_image(m_pdata->getImages(), access_location::host, access_mode::readwrite);
    const BoxDim& box = m_pdata->getBox();
    unsigned int group_size = m_group->getNumMembers();

    // update the position of each particle, and enforce period boundary condition
    for (unsigned int i = 0; i < group_size; i++)
        {
        double theta_i = h_orientation.data[i].x;
        h_pos.data[i].x += m_v0*slow::cos(theta_i);
        h_pos.data[i].y += m_v0*slow::sin(theta_i);
        box.wrap(h_pos.data[i], h_image.data[i]);
        }
    }


void export_TwoStepVicsek(py::module& m)
    {
    py::class_<TwoStepVicsek, std::shared_ptr<TwoStepVicsek> >(m, "TwoStepVicsek", py::base<TwoStepVicsekBase>())
        .def(py::init< std::shared_ptr<SystemDefinition>,
                       std::shared_ptr<ParticleGroup>,
                       std::shared_ptr<NeighborList>,
                       Scalar,
                       Scalar,
                       Scalar,
                       Scalar,
                       Scalar,
                       Scalar,
                       Scalar,
                       unsigned int>())
        ;
    }
