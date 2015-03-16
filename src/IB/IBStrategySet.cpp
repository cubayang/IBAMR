// Filename: IBStrategySet.cpp
// Created on 08 Mar 2012 by Boyce Griffith
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

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <string>
#include <vector>

#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/mesh/ChopAndPackLoadBalancer.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "ibamr/IBStrategy.h"
#include "ibamr/IBStrategySet.h"
#include "ibamr/ibamr_utilities.h"
#include "ibamr/namespaces.h" // IWYU pragma: keep
#include "SAMRAI/tbox/Database.h"

namespace IBTK
{
class RobinPhysBdryPatchStrategy;
} // namespace IBTK

namespace IBAMR
{
class IBHierarchyIntegrator;
} // namespace IBAMR
namespace SAMRAI
{
namespace tbox
{
template <class TYPE>
class Array;
} // namespace tbox
namespace xfer
{

class CoarsenSchedule;

class RefineSchedule;
} // namespace xfer
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBAMR
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

IBStrategySet::~IBStrategySet()
{
    // intentionally blank
    return;
}

void IBStrategySet::registerIBHierarchyIntegrator(IBHierarchyIntegrator* ib_solver)
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->registerIBHierarchyIntegrator(ib_solver);
    }
    return;
}

void IBStrategySet::registerEulerianVariables()
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->registerEulerianVariables();
    }
    return;
}

void IBStrategySet::registerEulerianCommunicationAlgorithms()
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->registerEulerianCommunicationAlgorithms();
    }
    return;
}

const IntVector& IBStrategySet::getMinimumGhostCellWidth() const
{
    static IntVector ghost_cell_width = IntVector::getZero(DIM);
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        ghost_cell_width = IntVector::max(ghost_cell_width, (*cit)->getMinimumGhostCellWidth());
    }
    return ghost_cell_width;
}

void IBStrategySet::setupTagBuffer(std::vector<int>& tag_buffer, boost::shared_ptr<PatchHierarchy> hierarchy) const
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->setupTagBuffer(tag_buffer, hierarchy);
    }
    return;
}

void IBStrategySet::preprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->preprocessIntegrateData(current_time, new_time, num_cycles);
    }
    return;
}

void IBStrategySet::postprocessIntegrateData(double current_time, double new_time, int num_cycles)
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->postprocessIntegrateData(current_time, new_time, num_cycles);
    }
    return;
}

void IBStrategySet::updateFixedLEOperators()
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->updateFixedLEOperators();
    }
    return;
}

void IBStrategySet::interpolateVelocity(int u_data_idx,
                                        const std::vector<boost::shared_ptr<CoarsenSchedule> >& u_synch_scheds,
                                        const std::vector<boost::shared_ptr<RefineSchedule> >& u_ghost_fill_scheds,
                                        double data_time)
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->interpolateVelocity(u_data_idx, u_synch_scheds, u_ghost_fill_scheds, data_time);
    }
    return;
}

void IBStrategySet::IBStrategySet::eulerStep(double current_time, double new_time)
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->eulerStep(current_time, new_time);
    }
    return;
}

void IBStrategySet::midpointStep(double current_time, double new_time)
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->midpointStep(current_time, new_time);
    }
    return;
}

void IBStrategySet::trapezoidalStep(double current_time, double new_time)
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->trapezoidalStep(current_time, new_time);
    }
    return;
}

void IBStrategySet::computeLagrangianForce(double data_time)
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->computeLagrangianForce(data_time);
    }
    return;
}

void IBStrategySet::spreadForce(int f_data_idx,
                                RobinPhysBdryPatchStrategy* f_phys_bdry_op,
                                const std::vector<boost::shared_ptr<RefineSchedule> >& f_prolongation_scheds,
                                double data_time)
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->spreadForce(f_data_idx, f_phys_bdry_op, f_prolongation_scheds, data_time);
    }
    return;
}

bool IBStrategySet::hasFluidSources() const
{
    bool has_fluid_sources = false;
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        has_fluid_sources = has_fluid_sources || (*cit)->hasFluidSources();
    }
    return has_fluid_sources;
}

void IBStrategySet::computeLagrangianFluidSource(double data_time)
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->computeLagrangianFluidSource(data_time);
    }
    return;
}

void IBStrategySet::spreadFluidSource(int q_data_idx,
                                      const std::vector<boost::shared_ptr<RefineSchedule> >& q_prolongation_scheds,
                                      double data_time)
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->spreadFluidSource(q_data_idx, q_prolongation_scheds, data_time);
    }
    return;
}

void IBStrategySet::interpolatePressure(int p_data_idx,
                                        const std::vector<boost::shared_ptr<CoarsenSchedule> >& p_synch_scheds,
                                        const std::vector<boost::shared_ptr<RefineSchedule> >& p_ghost_fill_scheds,
                                        double data_time)
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->interpolatePressure(p_data_idx, p_synch_scheds, p_ghost_fill_scheds, data_time);
    }
    return;
}

void IBStrategySet::preprocessSolveFluidEquations(double current_time, double new_time, int cycle_num)
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->preprocessSolveFluidEquations(current_time, new_time, cycle_num);
    }
    return;
}

void IBStrategySet::postprocessSolveFluidEquations(double current_time, double new_time, int cycle_num)
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->postprocessSolveFluidEquations(current_time, new_time, cycle_num);
    }
    return;
}

void IBStrategySet::postprocessData()
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->postprocessData();
    }
    return;
}

void IBStrategySet::initializePatchHierarchy(boost::shared_ptr<PatchHierarchy> hierarchy,
                                             boost::shared_ptr<GriddingAlgorithm> gridding_alg,
                                             int u_data_idx,
                                             const std::vector<boost::shared_ptr<CoarsenSchedule> >& u_synch_scheds,
                                             const std::vector<boost::shared_ptr<RefineSchedule> >& u_ghost_fill_scheds,
                                             int integrator_step,
                                             double init_data_time,
                                             bool initial_time)
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->initializePatchHierarchy(hierarchy, gridding_alg, u_data_idx, u_synch_scheds, u_ghost_fill_scheds,
                                         integrator_step, init_data_time, initial_time);
    }
    return;
}

void IBStrategySet::registerLoadBalancer(boost::shared_ptr<ChopAndPackLoadBalancer> load_balancer,
                                         int workload_data_idx)
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->registerLoadBalancer(load_balancer, workload_data_idx);
    }
    return;
}

void IBStrategySet::updateWorkloadEstimates(boost::shared_ptr<PatchHierarchy> hierarchy, int workload_data_idx)
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->updateWorkloadEstimates(hierarchy, workload_data_idx);
    }
    return;
}

void IBStrategySet::beginDataRedistribution(boost::shared_ptr<PatchHierarchy> hierarchy,
                                            boost::shared_ptr<GriddingAlgorithm> gridding_alg)
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->beginDataRedistribution(hierarchy, gridding_alg);
    }
    return;
}

void IBStrategySet::endDataRedistribution(boost::shared_ptr<PatchHierarchy> hierarchy,
                                          boost::shared_ptr<GriddingAlgorithm> gridding_alg)
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->endDataRedistribution(hierarchy, gridding_alg);
    }
    return;
}

void IBStrategySet::initializeLevelData(const boost::shared_ptr<PatchHierarchy>& hierarchy,
                                        int level_number,
                                        double init_data_time,
                                        bool can_be_refined,
                                        bool initial_time,
                                        const boost::shared_ptr<PatchLevel>& old_level,
                                        bool allocate_data)
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->initializeLevelData(hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level,
                                    allocate_data);
    }
    return;
}

void IBStrategySet::resetHierarchyConfiguration(const boost::shared_ptr<PatchHierarchy>& hierarchy,
                                                int coarsest_level,
                                                int finest_level)
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);
    }
    return;
}

void IBStrategySet::applyGradientDetector(const boost::shared_ptr<PatchHierarchy>& hierarchy,
                                          int level_number,
                                          double error_data_time,
                                          int tag_index,
                                          bool initial_time,
                                          bool uses_richardson_extrapolation_too)
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->applyGradientDetector(hierarchy, level_number, error_data_time, tag_index, initial_time,
                                      uses_richardson_extrapolation_too);
    }
    return;
}

void IBStrategySet::putToRestart(const boost::shared_ptr<Database>& db) const
{
    for (auto cit = d_strategy_set.begin(); cit != d_strategy_set.end(); ++cit)
    {
        (*cit)->putToRestart(db);
    }
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBAMR

//////////////////////////////////////////////////////////////////////////////
