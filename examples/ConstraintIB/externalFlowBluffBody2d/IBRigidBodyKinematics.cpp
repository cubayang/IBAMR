// Filename: IBRigidBodyKinematics.h
// Created by Amneet Bhalla on 12/20/2011.
// Modified by Olivier Mesnard on 6/30/2015.

// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
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


#include <string>

#include "IBRigidBodyKinematics.h"
#include <tbox/Utilities.h>
#include <tbox/PIO.h>
#include <tbox/SAMRAI_MPI.h>
#include "ibamr/namespaces.h"
#include "muParser.h"


namespace IBAMR
{
  
namespace
{
  
static const double PII =   3.14159265358979323846264338327950288419716939937510;
  
} // namespace unknown


IBRigidBodyKinematics::IBRigidBodyKinematics(const std::string& object_name,
                                             Pointer<Database> input_db,
                                             LDataManager* l_data_manager,
                                             bool register_for_restart)
    : ConstraintIBKinematics(object_name, input_db, l_data_manager, register_for_restart),
      d_current_time(0.0),
      d_kinematics_vel(NDIM),
      d_shape(NDIM),
      d_center_of_mass(3),
      d_incremented_angle_from_reference_axis(3),
      d_tagged_pt_position(3),
      d_parser_time(new double),
      d_parser_posn(new double[NDIM])
{
    // Read-in kinematics velocity functions
    std::vector<std::string> kinematicsvel_function_strings;
    for (int d = 0; d < NDIM; ++d)
    {
        std::ostringstream stream;
        stream << "_function_" << d;
        const std::string postfix = stream.str();
        std::string key_name = "kinematics_velocity" + postfix;
      
        if (input_db->isString(key_name))
        {       
            kinematicsvel_function_strings.push_back(input_db->getString(key_name));
        }
        else
        {
            kinematicsvel_function_strings.push_back("0.0");
            TBOX_WARNING("IBRigidBodyKinematics::IBRigidBodyKinematics() :\n" 
                         << " no function corresponding to key ``" << key_name 
                         <<  "'' found for dimension = " << d 
                         << "; using kinematics_vel = 0.0. " << std::endl); 
        }
    
        d_kinematicsvel_parsers.push_back(new mu::Parser());
        d_kinematicsvel_parsers.back()->SetExpr(kinematicsvel_function_strings.back());
        d_all_parsers.push_back(d_kinematicsvel_parsers.back());
    }
    
    // Define constants and variables for the parsers.
    for (std::vector<mu::Parser*>::const_iterator cit = d_all_parsers.begin(); cit != d_all_parsers.end(); ++cit)
    {   
        // Various names for pi.
        (*cit)->DefineConst("pi", PII);
        (*cit)->DefineConst("Pi", PII);
        (*cit)->DefineConst("PI", PII);
      
        // Variables.
        (*cit)->DefineVar("T", d_parser_time);
        (*cit)->DefineVar("t", d_parser_time);
        for (int d = 0; d < NDIM; ++d)
        {
            std::ostringstream stream;
            stream << d;
            const std::string postfix = stream.str();
            (*cit)->DefineVar("X"  + postfix, &(d_parser_posn[d]));
            (*cit)->DefineVar("x"  + postfix, &(d_parser_posn[d]));
            (*cit)->DefineVar("X_" + postfix, &(d_parser_posn[d]));
            (*cit)->DefineVar("x_" + postfix, &(d_parser_posn[d]));
        }  
    }

    // Set the size of vectors.
    const StructureParameters& struct_param = getStructureParameters();
    const int coarsest_ln                   = struct_param.getCoarsestLevelNumber();
    const int finest_ln                     = struct_param.getFinestLevelNumber();
    const int total_levels                  = finest_ln - coarsest_ln + 1; 
    d_kinematics_vel.resize(total_levels);
    
    const std::vector<std::pair<int,int> >& idx_range = struct_param.getLagIdxRange();
    for(int ln = 0 ; ln < total_levels ; ++ln)
    {
        const int nodes_this_ln = idx_range[ln].second - idx_range[ln].first;
        d_kinematics_vel[ln].resize(NDIM);
        for(int d = 0; d < NDIM; ++d)
        {
           d_kinematics_vel[ln][d].resize(nodes_this_ln);
        }
    }
        
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) 
    {
        getFromRestart();
        // At restart current and new states point to the same state.
        setRigidBodySpecificVelocity(d_current_time, d_incremented_angle_from_reference_axis, 
                                     d_center_of_mass, d_tagged_pt_position);
        setShape(d_current_time, d_incremented_angle_from_reference_axis);
    }

    return;  
} // IBRigidBodyKinematics

IBRigidBodyKinematics::~IBRigidBodyKinematics()
{
    for (std::vector<mu::Parser*>::const_iterator cit = d_all_parsers.begin(); cit != d_all_parsers.end(); ++cit)
    {
        delete (*cit);
    }
    delete d_parser_time;
    delete[] d_parser_posn;
    
    return;
} // ~IBRigidBodyKinematics

void IBRigidBodyKinematics::putToDatabase(Pointer<Database> db)
{
    db->putDouble("d_current_time", d_current_time);
    db->putDoubleArray("d_center_of_mass", &d_center_of_mass[0], 3);
    db->putDoubleArray("d_incremented_angle_from_reference_axis", &d_incremented_angle_from_reference_axis[0], 3);
    db->putDoubleArray("d_tagged_pt_position", &d_tagged_pt_position[0], 3);
   
    return;
} // putToDatabase

void IBRigidBodyKinematics::getFromRestart()
{  
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to " << d_object_name 
                                 << " not found in restart file." << std::endl);
    }
    
    d_current_time = db->getDouble("d_current_time");
    db->getDoubleArray("d_center_of_mass", &d_center_of_mass[0], 3);
    db->getDoubleArray("d_incremented_angle_from_reference_axis", &d_incremented_angle_from_reference_axis[0], 3);
    db->getDoubleArray("d_tagged_pt_position", &d_tagged_pt_position[0], 3);
  
    return;
} // getFromRestart

void IBRigidBodyKinematics::setRigidBodySpecificVelocity(const double time,
                                                         const std::vector<double>& incremented_angle_from_reference_axis,
                                                         const std::vector<double>& center_of_mass,
                                                         const std::vector<double>& tagged_pt_position)
{
    *d_parser_time = time;
    std::vector<double> vel_parser(NDIM);
    for(int d = 0; d < NDIM; ++d) 
    {
        d_parser_posn[d] = center_of_mass[d];
    }
    for(int d = 0; d < NDIM; ++d) 
    {
        vel_parser[d] = d_kinematicsvel_parsers[d]->Eval();
    }

    static const StructureParameters& struct_param = getStructureParameters();
    static const int coarsest_ln = struct_param.getCoarsestLevelNumber();
    static const int finest_ln = struct_param.getFinestLevelNumber();
    static const int total_levels = finest_ln - coarsest_ln + 1; 
    static const std::vector<std::pair<int,int> >& idx_range = struct_param.getLagIdxRange();

    for(int ln = 0 ; ln < total_levels ; ++ln)
    {
        const int nodes_this_ln = idx_range[ln].second - idx_range[ln].first;
        for(int d = 0; d < NDIM; ++d)
        {
            for(int idx = 0; idx < nodes_this_ln; ++idx)
            {    
                d_kinematics_vel[ln][d][idx] = vel_parser[d];
            }
        }
    }
  
    return;
} // setRigidBodySpecificVelocity


void IBRigidBodyKinematics::setKinematicsVelocity(const double time,
                                                  const std::vector<double>& incremented_angle_from_reference_axis,
                                                  const std::vector<double>& center_of_mass,
                                                  const std::vector<double>& tagged_pt_position)
{  
    d_new_time = time;
    d_incremented_angle_from_reference_axis = incremented_angle_from_reference_axis;
    d_center_of_mass = center_of_mass;
    d_tagged_pt_position = tagged_pt_position;
    
    setRigidBodySpecificVelocity(d_new_time,d_incremented_angle_from_reference_axis,d_center_of_mass,d_tagged_pt_position);
    
    return;
} // setKinematicsVelocity


const std::vector<std::vector<double> >& IBRigidBodyKinematics::getKinematicsVelocity(const int level) const
{
    static const StructureParameters& struct_param = getStructureParameters();
    static const int coarsest_ln  = struct_param.getCoarsestLevelNumber();
    
    static const int finest_ln    = struct_param.getFinestLevelNumber();
    TBOX_ASSERT(coarsest_ln <= level && level <= finest_ln);
         
    return d_kinematics_vel[level - coarsest_ln];
} // getKinematicsVelocity

void IBRigidBodyKinematics::setShape(const double time, const std::vector<double>& incremented_angle_from_reference_axis)
{
    //intentionally left blank
    return;
} // setShape

const std::vector<std::vector<double> >& IBRigidBodyKinematics::getShape(const int /*level*/) const
{
    return d_shape; 
} // getShape

} // namespace IBAMR
