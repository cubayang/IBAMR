// Filename: IBHierarchyIntegrator.h
// Created on 12 Jul 2004 by Boyce Griffith
//
// Copyright (c) 2002-2014, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#ifndef included_IBHierarchyIntegrator
#define included_IBHierarchyIntegrator

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <stddef.h>
#include <string>
#include <vector>

#include "SAMRAI/pdat/CellVariable.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/math/HierarchyDataOpsReal.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/mesh/ChopAndPackLoadBalancer.h"

#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "ibamr/IBStrategy.h"
#include "ibamr/INSHierarchyIntegrator.h"
#include "ibamr/ibamr_enums.h"
#include "ibtk/CartGridFunction.h"
#include "ibtk/HierarchyIntegrator.h"
#include "ibtk/LMarkerSetVariable.h"
#include "ibtk/ibtk_utilities.h"

namespace IBTK
{
class RobinPhysBdryPatchStrategy;
} // namespace IBTK

namespace SAMRAI
{
namespace hier
{

class PatchLevel;

class Patch;

class PatchHierarchy;

class PatchHierarchy;
} // namespace hier
namespace mesh
{

class GriddingAlgorithm;
} // namespace mesh
namespace tbox
{
class Database;
} // namespace tbox
} // namespace SAMRAI

/////////////////////////////// CLASS DEFINITION /////////////////////////////

namespace IBAMR
{
/*!
 * \brief Class IBHierarchyIntegrator provides an abstract interface for a time
 * integrator for various versions of the immersed boundary method on an AMR
 * grid hierarchy, along with basic data management for variables defined on
 * that hierarchy.
 */
class IBHierarchyIntegrator : public IBTK::HierarchyIntegrator
{
public:
    friend class IBStrategy;

    /*!
     * The destructor for class IBHierarchyIntegrator unregisters the integrator
     * object with the restart manager when the object is so registered.
     */
    ~IBHierarchyIntegrator();

    /*!
     * Return a pointer to the IBStrategy object registered with this
     * integrator.
     */
    boost::shared_ptr<IBStrategy> getIBStrategy() const;

    /*!
     * Supply a body force (optional).
     */
    void registerBodyForceFunction(const boost::shared_ptr<IBTK::CartGridFunction>& F_fcn);

    /*!
     * Register a load balancer for non-uniform load balancing.
     */
    void registerLoadBalancer(const boost::shared_ptr<SAMRAI::mesh::ChopAndPackLoadBalancer>& load_balancer);

    /*!
     * Return a pointer to the fluid velocity variable.
     */
    boost::shared_ptr<SAMRAI::hier::Variable> getVelocityVariable() const;

    /*!
     * Return a pointer to the fluid pressure state variable.
     */
    boost::shared_ptr<SAMRAI::hier::Variable> getPressureVariable() const;

    /*!
     * Return a pointer to the body force variable.
     */
    boost::shared_ptr<SAMRAI::hier::Variable> getBodyForceVariable() const;

    /*!
     * Return a pointer to the source strength variable.
     */
    boost::shared_ptr<SAMRAI::hier::Variable> getFluidSourceVariable() const;

    /*!
     * Initialize the variables, basic communications algorithms, solvers, and
     * other data structures used by this time integrator object.
     *
     * This method is called automatically by initializePatchHierarchy() prior
     * to the construction of the patch hierarchy.  It is also possible for
     * users to make an explicit call to initializeHierarchyIntegrator() prior
     * to calling initializePatchHierarchy().
     */
    void initializeHierarchyIntegrator(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                                       const boost::shared_ptr<SAMRAI::mesh::GriddingAlgorithm>& gridding_alg);

    /*!
     * Initialize the AMR patch hierarchy and data defined on the hierarchy at
     * the start of a computation.  If the computation is begun from a restart
     * file, the patch hierarchy and patch data are read from the hierarchy
     * database.  Otherwise, the patch hierarchy and patch data are initialized
     * by the gridding algorithm associated with the integrator object.
     *
     * The implementation of this function assumes that the hierarchy exists
     * upon entry to the function, but that it contains no patch levels.  On
     * return from this function, the state of the integrator object will be
     * such that it is possible to step through time via the advanceHierarchy()
     * function.
     */
    void initializePatchHierarchy(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                                  const boost::shared_ptr<SAMRAI::mesh::GriddingAlgorithm>& gridding_alg);

    /*!
     * Regrid the hierarchy.
     */
    void regridHierarchy();

protected:
    /*!
     * The constructor for class IBHierarchyIntegrator sets some default values,
     * reads in configuration information from input and restart databases, and
     * registers the integrator object with the restart manager when requested.
     */
    IBHierarchyIntegrator(const std::string& object_name,
                          const boost::shared_ptr<SAMRAI::tbox::Database>& input_db,
                          const boost::shared_ptr<IBStrategy>& ib_method_ops,
                          const boost::shared_ptr<INSHierarchyIntegrator>& ins_hier_integrator,
                          bool register_for_restart = true);

    /*!
     * Function to determine whether regridding should occur at the current time
     * step.
     */
    bool atRegridPointSpecialized() const;

    /*!
     * Initialize data on a new level after it is inserted into an AMR patch
     * hierarchy by the gridding algorithm.
     */
    void initializeLevelDataSpecialized(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                                        int level_number,
                                        double init_data_time,
                                        bool can_be_refined,
                                        bool initial_time,
                                        const boost::shared_ptr<SAMRAI::hier::PatchLevel>& old_level,
                                        bool allocate_data);

    /*!
     * Reset cached hierarchy dependent data.
     */
    void resetHierarchyConfigurationSpecialized(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                                                int coarsest_level,
                                                int finest_level);

    /*!
     * Set integer tags to "one" in cells where refinement of the given level
     * should occur according to the magnitude of the fluid vorticity.
     */
    void applyGradientDetectorSpecialized(const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                                          int level_number,
                                          double error_data_time,
                                          int tag_index,
                                          bool initial_time,
                                          bool uses_richardson_extrapolation_too);

    /*!
     * Write out specialized object state to the given database.
     */
    void putToRestartSpecialized(const boost::shared_ptr<SAMRAI::tbox::Database>& db) const;

    /*
     * Boolean value that indicates whether the integrator has been initialized.
     */
    bool d_integrator_is_initialized;

    /*!
     * Enum indicating the time integration employed for the IB equations.
     */
    TimeSteppingType d_time_stepping_type;

    /*!
     * Flags to determine whether warnings or error messages should be emitted
     * when time step size changes are encountered.
     */
    bool d_error_on_dt_change, d_warn_on_dt_change;

    /*
     * The (optional) INSHierarchyIntegrator is used to provide time integration
     * capability for the incompressible Navier-Stokes equations.
     */
    boost::shared_ptr<INSHierarchyIntegrator> d_ins_hier_integrator;

    /*
     * The regrid CFL interval indicates the number of meshwidths a particle may
     * move in any coordinate direction between invocations of the regridding
     * process.
     *
     * NOTE: Currently, when the CFL-based regrid interval is specified, it is
     * always used instead of the fixed-step regrid interval.
     */
    double d_regrid_cfl_interval, d_regrid_cfl_estimate;

    /*
     * IB method implementation object.
     */
    boost::shared_ptr<IBStrategy> d_ib_method_ops;

    /*
     * Hierarchy operations objects.
     */
    boost::shared_ptr<SAMRAI::math::HierarchyDataOpsReal<double> > d_hier_velocity_data_ops;
    boost::shared_ptr<SAMRAI::math::HierarchyDataOpsReal<double> > d_hier_pressure_data_ops;
    boost::shared_ptr<SAMRAI::math::HierarchyCellDataOpsReal<double> > d_hier_cc_data_ops;

    /*
     * Eulerian variables.
     */
    boost::shared_ptr<SAMRAI::hier::Variable> d_u_var, d_p_var, d_f_var, d_q_var;
    int d_u_idx, d_p_idx, d_f_idx, d_f_current_idx, d_q_idx;
    boost::shared_ptr<SAMRAI::hier::VariableContext> d_ib_context;

    /*
     * Refine and coarsen algorithm data.
     */
    boost::shared_ptr<IBTK::RobinPhysBdryPatchStrategy> d_u_phys_bdry_op, d_p_phys_bdry_op;
    boost::shared_ptr<SAMRAI::xfer::RefineAlgorithm> d_u_ghostfill_alg, d_f_prolong_alg, d_p_ghostfill_alg,
        d_q_prolong_alg;
    boost::shared_ptr<SAMRAI::hier::RefineOperator> d_u_ghostfill_op, d_f_prolong_op, d_p_ghostfill_op, d_q_prolong_op;

    boost::shared_ptr<SAMRAI::xfer::CoarsenAlgorithm> d_u_coarsen_alg, d_p_coarsen_alg;
    boost::shared_ptr<SAMRAI::hier::CoarsenOperator> d_u_coarsen_op, d_p_coarsen_op;

    /*
     * Body force functions.
     */
    boost::shared_ptr<IBTK::CartGridFunction> d_body_force_fcn;

    /*
     * Nonuniform load balancing data structures.
     */
    boost::shared_ptr<SAMRAI::mesh::ChopAndPackLoadBalancer> d_load_balancer;
    boost::shared_ptr<SAMRAI::pdat::CellVariable<double> > d_workload_var;
    int d_workload_idx;

    /*
     * Lagrangian marker data structures.
     */
    boost::shared_ptr<IBTK::LMarkerSetVariable> d_mark_var;
    int d_mark_current_idx, d_mark_new_idx, d_mark_scratch_idx;
    std::vector<IBTK::Point> d_mark_init_posns;
    std::string d_mark_file_name;

    /*!
     * \brief A class to communicate the Eulerian body force computed by class
     * IBHierarchyIntegrator to the incompressible Navier-Stokes solver.
     */
    class IBEulerianForceFunction : public IBTK::CartGridFunction
    {
    public:
        /*!
         * \brief Constructor.
         */
        IBEulerianForceFunction(const IBHierarchyIntegrator* ib_solver);

        /*!
         * \brief Destructor.
         */
        ~IBEulerianForceFunction();

        /*!
         * \name Methods to set the data.
         */
        //\{

        /*!
         * \note This concrete IBTK::CartGridFunction is time-dependent.
         */
        bool isTimeDependent() const;

        /*!
         * \brief Set the data on the patch interiors on the specified levels of
         * the patch hierarchy.
         */
        void setDataOnPatchHierarchy(const int data_idx,
                                     const boost::shared_ptr<SAMRAI::hier::Variable>& var,
                                     const boost::shared_ptr<SAMRAI::hier::PatchHierarchy>& hierarchy,
                                     const double data_time,
                                     const bool initial_time = false,
                                     const int coarsest_ln = -1,
                                     const int finest_ln = -1);

        /*!
         * Set the data on the patch interior.
         */
        void setDataOnPatch(int data_idx,
                            const boost::shared_ptr<SAMRAI::hier::Variable>& var,
                            const boost::shared_ptr<SAMRAI::hier::Patch>& patch,
                            double data_time,
                            bool initial_time = false,
                            const boost::shared_ptr<SAMRAI::hier::PatchLevel>& level = NULL);

        //\}

    private:
        /*!
         * \brief Default constructor.
         *
         * \note This constructor is not implemented and should not be used.
         */
        IBEulerianForceFunction();

        /*!
         * \brief Copy constructor.
         *
         * \note This constructor is not implemented and should not be used.
         *
         * \param from The value to copy to this object.
         */
        IBEulerianForceFunction(const IBEulerianForceFunction& from);

        /*!
         * \brief Assignment operator.
         *
         * \note This operator is not implemented and should not be used.
         *
         * \param that The value to assign to this object.
         *
         * \return A reference to this object.
         */
        IBEulerianForceFunction& operator=(const IBEulerianForceFunction& that);

        const IBHierarchyIntegrator* const d_ib_solver;
    };

    friend class IBEulerianForceFunction;

    /*!
     * \brief A class to communicate the Eulerian fluid source-sink distribution
     * computed by class IBHierarchyIntegrator to the incompressible
     * Navier-Stokes solver.
     */
    class IBEulerianSourceFunction : public IBTK::CartGridFunction
    {
    public:
        /*!
         * \brief Constructor.
         */
        IBEulerianSourceFunction(const IBHierarchyIntegrator* ib_solver);

        /*!
         * \brief Destructor.
         */
        ~IBEulerianSourceFunction();

        /*!
         * \name Methods to set the data.
         */
        //\{

        /*!
         * \note This concrete IBTK::CartGridFunction is time-dependent.
         */
        bool isTimeDependent() const;

        /*!
         * Set the data on the patch interior.
         */
        void setDataOnPatch(int data_idx,
                            const boost::shared_ptr<SAMRAI::hier::Variable>& var,
                            const boost::shared_ptr<SAMRAI::hier::Patch>& patch,
                            double data_time,
                            bool initial_time = false,
                            const boost::shared_ptr<SAMRAI::hier::PatchLevel>& level = NULL);

        //\}

    private:
        /*!
         * \brief Default constructor.
         *
         * \note This constructor is not implemented and should not be used.
         */
        IBEulerianSourceFunction();

        /*!
         * \brief Copy constructor.
         *
         * \note This constructor is not implemented and should not be used.
         *
         * \param from The value to copy to this object.
         */
        IBEulerianSourceFunction(const IBEulerianSourceFunction& from);

        /*!
         * \brief Assignment operator.
         *
         * \note This operator is not implemented and should not be used.
         *
         * \param that The value to assign to this object.
         *
         * \return A reference to this object.
         */
        IBEulerianSourceFunction& operator=(const IBEulerianSourceFunction& that);

        const IBHierarchyIntegrator* const d_ib_solver;
    };

    friend class IBEulerianSourceFunction;

private:
    /*!
     * \brief Default constructor.
     *
     * \note This constructor is not implemented and should not be used.
     */
    IBHierarchyIntegrator();

    /*!
     * \brief Copy constructor.
     *
     * \note This constructor is not implemented and should not be used.
     *
     * \param from The value to copy to this object.
     */
    IBHierarchyIntegrator(const IBHierarchyIntegrator& from);

    /*!
     * \brief Assignment operator.
     *
     * \note This operator is not implemented and should not be used.
     *
     * \param that The value to assign to this object.
     *
     * \return A reference to this object.
     */
    IBHierarchyIntegrator& operator=(const IBHierarchyIntegrator& that);

    /*!
     * Read input values from a given database.
     */
    void getFromInput(const boost::shared_ptr<SAMRAI::tbox::Database>& db, bool is_from_restart);

    /*!
     * Read object state from the restart file and initialize class data
     * members.
     */
    void getFromRestart();
};
} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////

#endif //#ifndef included_IBHierarchyIntegrator
