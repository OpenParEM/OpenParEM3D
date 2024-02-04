////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//    OpenParEM3D - A fullwave 3D electromagnetic simulator.                  //
//    Copyright (C) 2022 Brian Young                                          //
//                                                                            //
//    This program is free software: you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation, either version 3 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.   //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

// Structure for maintaining and using project setup information in OpenParEM tools.
// This is a structure rather than a class so that it can be passed to C programs.

#ifndef PROJECT_H
#define PROJECT_H

#include "petscsys.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inputFrequency.h"

// for the refinement variable
#define NONE 1
#define ALPHA 2
#define BETA 3
#define ABSGAMMA 4
#define ABSZO 5
#define REZO 6
#define IMZO 7

char* removeNewLineChar (char *);
char *removeComment (char *);
int removeQuote (char *);
void removeQuotes (char *);

struct projectData {
   char* version_name;       // set in init_project 
   char* version_value;      // set in init_project

   // user settable

   int project_calculate_poynting;
   int project_save_fields;

   char *mesh_file;
   int mesh_order;
   int mesh_save_refined;
   double mesh_refinement_port_limit;
   double mesh_2D_refinement_fraction;
   double mesh_3D_refinement_fraction;
   double mesh_quality_limit;

   char *port_definition_file;

   char *materials_global_path;
   char *materials_global_name;
   char *materials_local_path;
   char *materials_local_name;
   char *materials_default_boundary;
   int materials_check_limits;

   char *refinement_frequency;
   char *refinement_variable;
   int refinement_iteration_min;
   int refinement_iteration_max;
   int refinement_required_passes;
   double refinement_tolerance;

   struct inputFrequencyPlan *inputFrequencyPlans;
   unsigned long int inputFrequencyPlansAllocated;
   unsigned long int inputFrequencyPlansCount;

   double reference_impedance;
   char *touchstone_version;
   char *touchstone_frequency_unit;
   char *touchstone_format;

   int solution_modes;
   double solution_temperature;
   double solution_2D_tolerance;
   double solution_3D_tolerance;
   int solution_iteration_limit;
   int solution_modes_buffer;
   char *solution_impedance_definition;
   char *solution_impedance_calculation;
   int solution_check_closed_loop;
   int solution_accurate_residual;
   int solution_shift_invert;
   int solution_use_initial_guess;
   double solution_shift_factor;

   int output_show_refining_mesh;
   int output_show_postprocessing;
   int output_show_iterations;
   int output_show_license;

   int test_create_cases;
   int test_show_detailed_cases;

   int debug_show_memory;
   int debug_show_project;
   int debug_show_frequency_plan;
   int debug_show_materials;
   int debug_show_port_definitions;
   int debug_show_impedance_details;
   int debug_save_port_fields;
   int debug_tempfiles_keep;

   int field_points_count;
   int field_points_allocated;
   double *field_points_x;
   double *field_points_y;
   double *field_points_z;

   // not user settable

   char *project_name;
};

#endif
