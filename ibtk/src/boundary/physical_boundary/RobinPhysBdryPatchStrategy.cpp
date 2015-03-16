// Filename: RobinPhysBdryPatchStrategy.cpp
// Created on 28 Apr 2012 by Boyce Griffith
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

#include <ostream>
#include <set>
#include <vector>

#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/IntVector.h"
#include "ibtk/RobinPhysBdryPatchStrategy.h"
#include "ibtk/ibtk_utilities.h"
#include "ibtk/namespaces.h" // IWYU pragma: keep
#include "SAMRAI/tbox/Utilities.h"

namespace SAMRAI
{
namespace hier
{

class Patch;
} // namespace hier
namespace solv
{

class RobinBcCoefStrategy;
} // namespace solv
} // namespace SAMRAI

/////////////////////////////// NAMESPACE ////////////////////////////////////

namespace IBTK
{
/////////////////////////////// STATIC ///////////////////////////////////////

/////////////////////////////// PUBLIC ///////////////////////////////////////

RobinPhysBdryPatchStrategy::RobinPhysBdryPatchStrategy()
    : RefinePatchStrategy(), d_patch_data_indices(), d_bc_coefs(), d_homogeneous_bc(false)
{
    // intentionally blank
    return;
}

RobinPhysBdryPatchStrategy::~RobinPhysBdryPatchStrategy()
{
    // intentionally blank
    return;
}

void RobinPhysBdryPatchStrategy::setPatchDataIndex(const int patch_data_index)
{
    std::set<int> patch_data_indices;
    patch_data_indices.insert(patch_data_index);
    setPatchDataIndices(patch_data_indices);
    return;
}

void RobinPhysBdryPatchStrategy::setPatchDataIndices(const std::set<int>& patch_data_indices)
{
    d_patch_data_indices.clear();
    d_patch_data_indices = patch_data_indices;
    return;
}

void RobinPhysBdryPatchStrategy::setPatchDataIndices(const ComponentSelector& patch_data_indices)
{
    std::set<int> patch_data_index_set;
    for (int l = 0; l < patch_data_indices.getSize(); ++l)
    {
        if (patch_data_indices.isSet(l))
        {
            const int patch_data_index = l;
            patch_data_index_set.insert(patch_data_index);
        }
    }
    setPatchDataIndices(patch_data_index_set);
    return;
}

void RobinPhysBdryPatchStrategy::setPhysicalBcCoef(boost::shared_ptr<RobinBcCoefStrategy> const bc_coef)
{
    setPhysicalBcCoefs(std::vector<boost::shared_ptr<RobinBcCoefStrategy>>(1, bc_coef));
    return;
}

void RobinPhysBdryPatchStrategy::setPhysicalBcCoefs(const std::vector<boost::shared_ptr<RobinBcCoefStrategy>>& bc_coefs)
{
    for (unsigned int l = 0; l < bc_coefs.size(); ++l)
    {
        TBOX_ASSERT(bc_coefs[l]);
    }
    d_bc_coefs = bc_coefs;
    return;
}

void RobinPhysBdryPatchStrategy::setHomogeneousBc(bool homogeneous_bc)
{
    d_homogeneous_bc = homogeneous_bc;
    return;
}

bool RobinPhysBdryPatchStrategy::getHomogeneousBc() const
{
    return d_homogeneous_bc;
}

void RobinPhysBdryPatchStrategy::preprocessRefine(Patch& /*fine*/,
                                                  const Patch& /*coarse*/,
                                                  const Box& /*fine_box*/,
                                                  const IntVector& /*ratio*/)
{
    // intentionally blank
    return;
}

void RobinPhysBdryPatchStrategy::postprocessRefine(Patch& /*fine*/,
                                                   const Patch& /*coarse*/,
                                                   const Box& /*fine_box*/,
                                                   const IntVector& /*ratio*/)
{
    // intentionally blank
    return;
}

void RobinPhysBdryPatchStrategy::accumulateFromPhysicalBoundaryData(Patch& /*patch*/,
                                                                    double /*fill_time*/,
                                                                    const IntVector& /*ghost_width_to_fill*/)
{
    TBOX_ERROR("RobinPhysBdryPatchStrategy::accumulateFromPhysicalBoundaryData(): unimplemented\n");
    return;
}

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////

/////////////////////////////// NAMESPACE ////////////////////////////////////

} // namespace IBTK

//////////////////////////////////////////////////////////////////////////////
