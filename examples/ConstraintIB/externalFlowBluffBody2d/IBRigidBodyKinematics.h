// Filename: IBRigidBodyKinematics.h
// Created by Amneet Bhalla on 1/10/2012
// Modified by Olivier Mesnard on 6/30/2015

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

// This is a concrete class which provides kinematics of a rigid solid to ConstraintIBMethod class.
     
 
#ifndef included_IBRigidBodyKinematics
#define included_IBRigidBodyKinematics

#include <iostream>
#include <vector>
#include <map>

#include "ibamr/ConstraintIBKinematics.h"


namespace mu
{
class Parser;
}// namespace mu

namespace IBAMR
{
/*!
 * \brief IBRigigBodyKinematics is a concrete class that provides definition for
 * the base ConstraintIBKinematics class.
 */
class IBRigidBodyKinematics : public ConstraintIBKinematics
{ 
public:
    /*!
     * \brief Constructor.
     */
    IBRigidBodyKinematics(const std::string& object_name,
                          SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                          IBTK::LDataManager* l_data_manager,
                          bool register_for_restart = true);
    
    /*!
     * Destructor.
     */
    virtual ~IBRigidBodyKinematics();
    
    /*!
     * \brief Set kinematics velocity for rigid body.
     * \see IBAMR::ConstraintIBKinematics::setKinematicsVelocity
     */
    virtual void setKinematicsVelocity(const double time,
                                       const std::vector<double>& incremented_angle_from_reference_axis,
                                       const std::vector<double>& center_of_mass,
                                       const std::vector<double>& tagged_pt_position);
    
    /*!
     * \brief Get the kinematics velocity for rigid body on the specified level.
     * \see IBAMR::ConstraintIBKinematics::getKinematicsVelocity
     */
    virtual const std::vector<std::vector<double> >& getKinematicsVelocity(const int level) const;    
  
    /*!
     * \brief Set the shape of the rigid body at the required time.
     * \see IBAMR::ConstraintIBKinematics::setShape
     */
    virtual void setShape(const double time, const std::vector<double>& incremented_angle_from_reference_axis);
    
    /*!
     * \brief Get the shape of the rigid body on the required level.
     */
    virtual const std::vector<std::vector<double> >& getShape(const int level) const;

    /*!
     * \brief Override the ConstraintIBKinematics base class method.
     */
    virtual void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);
    
private:
    /*!
     * \brief The default constructor is not implemented and should not be used.
     */
    IBRigidBodyKinematics();
    
    /*!
     * \brief The copy constructor is not implemented and should not be used.
     */
    IBRigidBodyKinematics(const IBRigidBodyKinematics& from);
  
    /*!
     * \brief The assignment operator is not implemented and should not be used.
     */
    IBRigidBodyKinematics& operator = (const IBRigidBodyKinematics& that);

    /*!
     * \brief Get data from restart.
     */
    void getFromRestart();

    /*!
     * \brief Set deformation kinematics velocity of the rogid body.
     */
    void setRigidBodySpecificVelocity(const double time,
                                      const std::vector<double>& incremented_angle_from_reference_axis,
                                      const std::vector<double>& center_of_mass,
                                      const std::vector<double>& tagged_pt_position);

    /*!
     * Current time (t) and new time (t+dt).
     */
    double d_current_time, d_new_time;

    /*!
     * Deformational velocity and shape vectors.
     */
    std::vector<std::vector<std::vector<double> > > d_kinematics_vel; 
    std::vector<std::vector<double> > d_shape;

    /*!
     * Save COM, tagged point position and incremented angle from reference axis for restarted runs.
     */
    std::vector<double> d_center_of_mass, d_incremented_angle_from_reference_axis, d_tagged_pt_position;

    /*!
     * The mu::Parser objects which evaluate the data-setting functions.
     */
    std::vector<mu::Parser*> d_kinematicsvel_parsers;
    std::vector<mu::Parser*> d_all_parsers;

    /*!
     * Time and position variables.
     */
    double* d_parser_time;    
    double* d_parser_posn;
   
}; // IBRigidBodyKinematics 
  
} // IBAMR
#endif //#ifndef included_IBRigidBodyKinematics
