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

// Config files
#include <IBAMR_config.h>
#include <IBTK_config.h>
#include <SAMRAI/SAMRAI_config.h>

// Headers for basic PETSc functions
#include <petscsys.h>

// Headers for basic SAMRAI objects
#include <SAMRAI/algs/HyperbolicLevelIntegrator.h>
#include <SAMRAI/geom/CartesianGridGeometry.h>
#include <SAMRAI/mesh/BergerRigoutsos.h>
#include <SAMRAI/mesh/ChopAndPackLoadBalancer.h>
#include <SAMRAI/mesh/StandardTagAndInitialize.h>
#include <SAMRAI/tbox/RestartManager.h>

// Headers for application-specific algorithm/data structure objects
#include <SAMRAI/algs/TimeRefinementIntegrator.h>
#include <SAMRAI/solv/LocationIndexRobinBcCoefs.h>
#include <ibamr/AdvectorPredictorCorrectorHyperbolicPatchOps.h>
#include <ibamr/AdvectorExplicitPredictorPatchOps.h>
#include <ibamr/app_namespaces.h>
#include <ibtk/AppInitializer.h>
#include <ibtk/HierarchyMathOps.h>
#include "QInit.h"
#include "UFunction.h"

/*******************************************************************************
 * For each run, the input filename and restart information (if needed) must   *
 * be given on the command line.  For non-restarted case, command line is:     *
 *                                                                             *
 *    executable <input file name>                                             *
 *                                                                             *
 * For restarted run, command line is:                                         *
 *                                                                             *
 *    executable <input file name> <restart directory> <restart number>        *
 *                                                                             *
 *******************************************************************************/
int main(int argc, char* argv[])
{
    // Initialize MPI and SAMRAI.
    SAMRAI_MPI::init(&argc, &argv);
    SAMRAI_MPI::setCallAbortInSerialInsteadOfExit();
    SAMRAIManager::startup();

    { // cleanup dynamically allocated objects prior to shutdown

        // Parse command line options, set some standard options from the input
        // file, initialize the restart database (if this is a restarted run),
        // and enable file logging.
        auto app_initializer = boost::make_shared<AppInitializer>(argc, argv, "advect.log");
        auto input_db = app_initializer->getInputDatabase();
        auto main_db = app_initializer->getComponentDatabase("Main");

        // Get various standard options set in the input file.
        const bool dump_viz_data = app_initializer->dumpVizData();
        const int viz_dump_interval = app_initializer->getVizDumpInterval();
        const bool uses_visit = dump_viz_data && app_initializer->getVisItDataWriter();

        const bool dump_restart_data = app_initializer->dumpRestartData();
        const int restart_dump_interval = app_initializer->getRestartDumpInterval();
        const string restart_dump_dirname = app_initializer->getRestartDumpDirectory();

        const bool dump_timer_data = app_initializer->dumpTimerData();
        const int timer_dump_interval = app_initializer->getTimerDumpInterval();

        // Get solver configuration options.
        bool using_refined_timestepping = false;
        if (main_db->keyExists("timestepping"))
        {
            string timestepping_method = main_db->getString("timestepping");
            if (timestepping_method == "SYNCHRONIZED")
            {
                using_refined_timestepping = false;
            }
            else
            {
                using_refined_timestepping = true;
            }
        }
        if (using_refined_timestepping)
        {
            pout << "using subcycled timestepping.\n";
        }
        else
        {
            pout << "NOT using subcycled timestepping.\n";
        }

        // Create major algorithm and data objects that comprise the
        // application.  These objects are configured from the input database
        // and, if this is a restarted run, from the restart database.
        auto explicit_predictor = boost::make_shared<AdvectorExplicitPredictorPatchOps>(
            "AdvectorExplicitPredictorPatchOps",
            app_initializer->getComponentDatabase("AdvectorExplicitPredictorPatchOps"));
        auto grid_geometry = boost::make_shared<CartesianGridGeometry>(
            DIM, "CartesianGeometry", app_initializer->getComponentDatabase("CartesianGeometry"));
        auto hyp_patch_ops = boost::make_shared<AdvectorPredictorCorrectorHyperbolicPatchOps>(
            "AdvectorPredictorCorrectorHyperbolicPatchOps",
            app_initializer->getComponentDatabase("AdvectorPredictorCorrectorHyperbolicPatchOps"), explicit_predictor,
            grid_geometry);
        auto hyp_level_integrator = boost::make_shared<HyperbolicLevelIntegrator>(
            "HyperbolicLevelIntegrator", app_initializer->getComponentDatabase("HyperbolicLevelIntegrator"),
            hyp_patch_ops.get(), using_refined_timestepping);
        auto patch_hierarchy = boost::make_shared<PatchHierarchy>("PatchHierarchy", grid_geometry);
        auto error_detector = boost::make_shared<StandardTagAndInitialize>(
            "StandardTagAndInitialize", hyp_level_integrator.get(),
            app_initializer->getComponentDatabase("StandardTagAndInitialize"));
        auto box_generator = boost::make_shared<BergerRigoutsos>(DIM);
        auto load_balancer = boost::make_shared<ChopAndPackLoadBalancer>(
            DIM, "ChopAndPackLoadBalancer", app_initializer->getComponentDatabase("ChopAndPackLoadBalancer"));
        auto gridding_algorithm = boost::make_shared<GriddingAlgorithm>(
            patch_hierarchy, "GriddingAlgorithm", app_initializer->getComponentDatabase("GriddingAlgorithm"),
            error_detector, box_generator, load_balancer);
        auto time_integrator = boost::make_shared<TimeRefinementIntegrator>(
            "TimeRefinementIntegrator", app_initializer->getComponentDatabase("TimeRefinementIntegrator"),
            patch_hierarchy, hyp_level_integrator, gridding_algorithm);

        // Setup the advection velocity.
        const bool u_is_div_free = main_db->getBoolWithDefault("u_is_div_free", false);
        if (u_is_div_free)
        {
            pout << "advection velocity u is discretely divergence free.\n";
        }
        else
        {
            pout << "advection velocity u is NOT discretely divergence free.\n";
        }
        auto u_var = boost::make_shared<FaceVariable<double> >(DIM, "u");
        UFunction u_fcn("UFunction", grid_geometry, app_initializer->getComponentDatabase("UFunction"));
        hyp_patch_ops->registerAdvectionVelocity(u_var);
        hyp_patch_ops->setAdvectionVelocityIsDivergenceFree(u_var, u_is_div_free);
        hyp_patch_ops->setAdvectionVelocityFunction(u_var, boost::shared_ptr<CartGridFunction>(&u_fcn, NullDeleter()));

        // Setup the advected quantity.
        const ConvectiveDifferencingType difference_form =
            IBAMR::string_to_enum<ConvectiveDifferencingType>(main_db->getStringWithDefault(
                "difference_form", IBAMR::enum_to_string<ConvectiveDifferencingType>(ADVECTIVE)));
        pout << "solving the advection equation in " << enum_to_string<ConvectiveDifferencingType>(difference_form)
             << " form.\n";
        auto Q_var = boost::make_shared<CellVariable<double> >(DIM, "Q");
        QInit Q_init("QInit", grid_geometry, app_initializer->getComponentDatabase("QInit"));
        auto physical_bc_coef = boost::make_shared<LocationIndexRobinBcCoefs>(
            DIM, "physical_bc_coef", app_initializer->getComponentDatabase("LocationIndexRobinBcCoefs"));
        hyp_patch_ops->registerTransportedQuantity(Q_var);
        hyp_patch_ops->setAdvectionVelocity(Q_var, u_var);
        hyp_patch_ops->setConvectiveDifferencingType(Q_var, difference_form);
        hyp_patch_ops->setInitialConditions(Q_var, boost::shared_ptr<CartGridFunction>(&Q_init, NullDeleter()));
        hyp_patch_ops->setPhysicalBcCoefs(Q_var, physical_bc_coef);

        // Set up visualization plot file writer.
        auto visit_data_writer = app_initializer->getVisItDataWriter();
        if (uses_visit) hyp_patch_ops->registerVisItDataWriter(visit_data_writer);

        // Initialize hierarchy configuration and data on all patches.
        double dt_now = time_integrator->initializeHierarchy();

        // Deallocate initialization objects.
        app_initializer.reset();

        // Print the input database contents to the log file.
        plog << "Input database:\n";
        input_db->printClassData(plog);

        // Write out initial visualization data.
        int iteration_num = time_integrator->getIntegratorStep();
        double loop_time = time_integrator->getIntegratorTime();
        if (dump_viz_data && uses_visit)
        {
            pout << "\n\nWriting visualization files...\n\n";
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
        }

        // Main time step loop.
        double loop_time_end = time_integrator->getEndTime();
        while (!MathUtilities<double>::equalEps(loop_time, loop_time_end) && time_integrator->stepsRemaining())
        {
            iteration_num = time_integrator->getIntegratorStep();
            loop_time = time_integrator->getIntegratorTime();

            pout << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "At beginning of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";

            double dt_new = time_integrator->advanceHierarchy(dt_now);
            loop_time += dt_now;
            dt_now = dt_new;

            pout << "\n";
            pout << "At end       of timestep # " << iteration_num << "\n";
            pout << "Simulation time is " << loop_time << "\n";
            pout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
            pout << "\n";

            // At specified intervals, write visualization and restart files,
            // and print out timer data.
            iteration_num += 1;
            const bool last_step = !time_integrator->stepsRemaining();
            if (dump_viz_data && uses_visit && (iteration_num % viz_dump_interval == 0 || last_step))
            {
                pout << "\nWriting visualization files...\n\n";
                visit_data_writer->writePlotData(patch_hierarchy, iteration_num, loop_time);
            }
            if (dump_restart_data && (iteration_num % restart_dump_interval == 0 || last_step))
            {
                pout << "\nWriting restart files...\n\nn";
                RestartManager::getManager()->writeRestartFile(restart_dump_dirname, iteration_num);
            }
            if (dump_timer_data && (iteration_num % timer_dump_interval == 0 || last_step))
            {
                pout << "\nWriting timer data...\n\n";
                TimerManager::getManager()->print(plog);
            }
        }

        // Determine the accuracy of the computed solution.
        pout << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n"
             << "Computing error norms.\n\n";

        auto var_db = VariableDatabase::getDatabase();
        const auto Q_ctx = hyp_level_integrator->getCurrentContext();
        const int Q_idx = var_db->mapVariableAndContextToIndex(Q_var, Q_ctx);
        const int Q_cloned_idx = var_db->registerClonedPatchDataIndex(Q_var, Q_idx);

        const int coarsest_ln = 0;
        const int finest_ln = patch_hierarchy->getFinestLevelNumber();
        for (int ln = coarsest_ln; ln <= finest_ln; ++ln)
        {
            patch_hierarchy->getPatchLevel(ln)->allocatePatchData(Q_cloned_idx, loop_time);
        }
        Q_init.setDataOnPatchHierarchy(Q_cloned_idx, Q_var, patch_hierarchy, loop_time);

        HierarchyMathOps hier_math_ops("HierarchyMathOps", patch_hierarchy);
        hier_math_ops.setPatchHierarchy(patch_hierarchy);
        hier_math_ops.resetLevels(coarsest_ln, finest_ln);
        const int wgt_idx = hier_math_ops.getCellWeightPatchDescriptorIndex();

        HierarchyCellDataOpsReal<double> hier_cc_data_ops(patch_hierarchy, coarsest_ln, finest_ln);
        hier_cc_data_ops.subtract(Q_idx, Q_idx, Q_cloned_idx);
        pout << "Error in " << Q_var->getName() << " at time " << loop_time << ":\n"
             << "  L1-norm:  " << hier_cc_data_ops.L1Norm(Q_idx, wgt_idx) << "\n"
             << "  L2-norm:  " << hier_cc_data_ops.L2Norm(Q_idx, wgt_idx) << "\n"
             << "  max-norm: " << hier_cc_data_ops.maxNorm(Q_idx, wgt_idx) << "\n"
             << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

        if (dump_viz_data && uses_visit)
        {
            visit_data_writer->writePlotData(patch_hierarchy, iteration_num + 1, loop_time);
        }
    }

    SAMRAIManager::shutdown();
    SAMRAI_MPI::finalize();
    return 0;
}
