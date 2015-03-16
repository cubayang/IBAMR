// Filename: StandardTagAndInitStrategySet.cpp
// Created on 26 Jun 2007 by Boyce Griffith
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

#include <algorithm>
#include <limits>
#include <vector>

#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "ibtk/StandardTagAndInitStrategySet.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

StandardTagAndInitStrategySet::~StandardTagAndInitStrategySet()
{
    // intentionally blank
    return;
}

double StandardTagAndInitStrategySet::getLevelDt(const boost::shared_ptr<PatchLevel>& level,
                                                 const double dt_time,
                                                 const bool initial_time)
{
    double dt = std::numeric_limits<double>::max();
    for (auto it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        dt = std::min(dt, (*it)->getLevelDt(level, dt_time, initial_time));
    }
    return dt;
}

double StandardTagAndInitStrategySet::advanceLevel(const boost::shared_ptr<PatchLevel>& level,
                                                   const boost::shared_ptr<PatchHierarchy>& hierarchy,
                                                   const double current_time,
                                                   const double new_time,
                                                   const bool first_step,
                                                   const bool last_step,
                                                   const bool regrid_advance)
{
    double dt = std::numeric_limits<double>::max();
    for (auto it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        dt = std::min(
            dt, (*it)->advanceLevel(level, hierarchy, current_time, new_time, first_step, last_step, regrid_advance));
    }
    return dt;
}

void StandardTagAndInitStrategySet::resetTimeDependentData(const boost::shared_ptr<PatchLevel>& level,
                                                           const double new_time,
                                                           const bool can_be_refined)
{
    for (auto it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->resetTimeDependentData(level, new_time, can_be_refined);
    }
    return;
}

void StandardTagAndInitStrategySet::resetDataToPreadvanceState(const boost::shared_ptr<PatchLevel>& level)
{
    for (auto it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->resetDataToPreadvanceState(level);
    }
    return;
}

void StandardTagAndInitStrategySet::initializeLevelData(const boost::shared_ptr<PatchHierarchy>& hierarchy,
                                                        const int level_number,
                                                        const double init_data_time,
                                                        const bool can_be_refined,
                                                        const bool initial_time,
                                                        const boost::shared_ptr<PatchLevel>& old_level,
                                                        const bool allocate_data)
{
    for (auto it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->initializeLevelData(hierarchy, level_number, init_data_time, can_be_refined, initial_time, old_level,
                                   allocate_data);
    }
    return;
}

void StandardTagAndInitStrategySet::resetHierarchyConfiguration(const boost::shared_ptr<PatchHierarchy>& hierarchy,
                                                                const int coarsest_level,
                                                                const int finest_level)
{
    for (auto it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->resetHierarchyConfiguration(hierarchy, coarsest_level, finest_level);
    }
    return;
}

void StandardTagAndInitStrategySet::applyGradientDetector(const boost::shared_ptr<PatchHierarchy>& hierarchy,
                                                          const int level_number,
                                                          const double error_data_time,
                                                          const int tag_index,
                                                          const bool initial_time,
                                                          const bool uses_richardson_extrapolation_too)
{
    for (auto it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->applyGradientDetector(hierarchy, level_number, error_data_time, tag_index, initial_time,
                                     uses_richardson_extrapolation_too);
    }
    return;
}

void StandardTagAndInitStrategySet::applyRichardsonExtrapolation(const boost::shared_ptr<PatchLevel>& level,
                                                                 const double error_data_time,
                                                                 const int tag_index,
                                                                 const double deltat,
                                                                 const int error_coarsen_ratio,
                                                                 const bool initial_time,
                                                                 const bool uses_gradient_detector_too)
{
    for (auto it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->applyRichardsonExtrapolation(level, error_data_time, tag_index, deltat, error_coarsen_ratio,
                                            initial_time, uses_gradient_detector_too);
    }
    return;
}

void
StandardTagAndInitStrategySet::coarsenDataForRichardsonExtrapolation(const boost::shared_ptr<PatchHierarchy>& hierarchy,
                                                                     const int level_number,
                                                                     const boost::shared_ptr<PatchLevel>& coarser_level,
                                                                     const double coarsen_data_time,
                                                                     const bool before_advance)
{
    for (auto it = d_strategy_set.begin(); it != d_strategy_set.end(); ++it)
    {
        (*it)->coarsenDataForRichardsonExtrapolation(hierarchy, level_number, coarser_level, coarsen_data_time,
                                                     before_advance);
    }
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
