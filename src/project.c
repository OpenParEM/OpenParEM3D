////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//    OpenParEM3D - A fullwave 3D electromagnetic simulator.                  //
//    Copyright (C) 2024 Brian Young                                          //
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

// C version of the setup information data.
// Do this in C so that it can be passed to C programs.

#include "project.h"

// allocate memory and copy a string
char* allocCopyString (char *a) {
   char *b=NULL;

   if (a == NULL) return NULL;

   b=(char *)malloc((strlen(a)+1)*sizeof(char));
   if (b == NULL) return NULL;
   strcpy(b,a);

   return b;
}

char* removeNewLineChar (char *a) {
   int i;

   if (a == NULL) return NULL;

   i=0;
   while (i < strlen(a)) {
      if (a[i] == '\n') {a[i]='\0'; break;}
      i++;
   }

   return a;
}

char *removeComment (char *a) {
   long int i;

   if (a == NULL) return NULL;

   // remove comment 
   i=0;
   while (i < strlen(a)-1) {
      if (a[i] == '/' && a[i+1] == '/') {a[i]='\0'; break;}
      i++;
   }

   // remove trailing spaces
   i=strlen(a)-1;
   while (i >= 0) {
      if (a[i] != ' ' && a[i] != '\t') {a[i+1]='\0'; break;}
      i--;
   }

   return a;
}

int removeSpace (char *a) {
   int i,j,len;

   if (a == NULL) return 0;

   len=strlen(a);
   i=0;
   while (i < len) {
      if (a[i] == ' ' || a[i] == '\t') {
         j=i;
         while (j < len) {
            a[j]=a[j+1];
            j++;
         }
         return 1;
      }
      i++;
   }

   return 0;
}

void removeSpaces (char *a) {
   while (removeSpace(a));
}

int is_blank (char *a) {
   int i=0;
   while (i < strlen(a)) {
      if (a[i] != ' ' && a[i] != '\n' && a[i] != '\t') return 0;
      i++;
   }
   return 1;
}

int is_text (char *a) {
   if (a == NULL) return 0;
   if (strlen(a) == 0) return 0;
   return 1;
}

int is_true (char *a) {
   if (a == NULL) return 0;
   if (strcmp(a,"true") == 0) return 1;
   if (strcmp(a,"TRUE") == 0) return 1;
   if (strcmp(a,"1") == 0) return 1;
   if (strcmp(a,"True") == 0) return 1;
   return 0;
}

int is_false (char *a) {
   if (a == NULL) return 0;
   if (strcmp(a,"false") == 0) return 1;
   if (strcmp(a,"FALSE") == 0) return 1;
   if (strcmp(a,"0") == 0) return 1;
   if (strcmp(a,"False") == 0) return 1;
   return 0;
}

int is_bool (char *a) {
   if (is_true(a)) return 1;
   if (is_false(a)) return 1;
   return 0;
}

int is_int (char *a) {
   int found;

   if (a == NULL) return 0;
   if (strlen(a) == 0) return 0;

   // skip leading white space
   int i=0;
   while (i < strlen(a)) {
      if (a[i] != ' ' && a[i] != '\t') break;
      i++;
   }

   // skip a leading sign
   if (i < strlen(a) && (a[i] == '-' || a[i] == '+')) i++;

   while (i < strlen(a)) {
      found=0;

      if (a[i] == '0' ||
          a[i] == '1' ||
          a[i] == '2' ||
          a[i] == '3' ||
          a[i] == '4' ||
          a[i] == '5' ||
          a[i] == '6' ||
          a[i] == '7' ||
          a[i] == '8' ||
          a[i] == '9') found=1;

      if (!found) return 0;
      
      i++;
   }
   return 1;
}

int is_double (char *a) {
   int found;
   int found_e=-1;
   int found_period=-1;
   int found_sign=-1;

   if (a == NULL) return 0;
   if (strlen(a) == 0) return 0;

   // skip leading white space
   int i=0;
   while (i < strlen(a)) {
      if (a[i] != ' ' && a[i] != '\t') break;
      i++;
   }

   // skip a leading sign
   if (i < strlen(a) && (a[i] == '-' || a[i] == '+')) i++;

   while (i < strlen(a)) {
      found=0;
      if (a[i] == '0' ||
          a[i] == '1' ||
          a[i] == '2' ||
          a[i] == '3' ||
          a[i] == '4' ||
          a[i] == '5' ||
          a[i] == '6' ||
          a[i] == '7' ||
          a[i] == '8' ||
          a[i] == '9') found=1;

      if (a[i] == 'e' || a[i] == 'E') {
         if (found_e >= 0) return 0;
         found_e=i;
         found=1;
      }

      if (a[i] == '-' || a[i] == '+') {
         if (found_sign >= 0) return 0;
         found_sign=i;
         found=1;
      }

      if (a[i] == '.') {
         if (found_period >= 0) return 0;
         found_period=i;
         found=1;
      }

      if (!found) return 0;

      i++;
   }

   if (found_period >= 0 && found_e >= 0 && found_period > found_e) return 0;
   if (found_sign >= 0) {
      if (found_e >= 0) {
         if (found_sign != found_e+1) return 0;
      } else return 0;
   }
   if (found_sign == strlen(a)-1) return 0;
   if (found_e == strlen(a)-1) return 0;

   return 1;
}

int comma_count (char *a)
{
   int i,count;

   if (a == NULL) return 0;

   count=0;
   i=0;
   while (i < strlen(a)) {
      if (a[i] == ',') count++;
      i++;
   }
   return count;
}

int double_compare (double a, double b, double tol)
{
   if (a == b) return 1;
   if (a == 0 && fabs(b) < tol) return 1;
   if (b == 0 && fabs(a) < tol) return 1;
   if (fabs((a-b)/a) < tol) return 1;
   return 0;
}

int is_refinement_frequency (char *a)
{
   if (a == NULL) return 0;
   if (strcmp(a,"all") == 0) return 1;
   if (strcmp(a,"none") == 0) return 1;
   if (strcmp(a,"high") == 0) return 1;
   if (strcmp(a,"low") == 0) return 1;
   if (strcmp(a,"highlow") == 0) return 1;
   if (strcmp(a,"lowhigh") == 0) return 1;
   if (strcmp(a,"plan") == 0) return 1;
   return 0;
}

int is_refinement_variable (char *a)
{
   if (a == NULL) return 0;
   if (strcmp(a,"S") == 0) return 1;
   if (strcmp(a,"SorH") == 0) return 1;
   if (strcmp(a,"SandH") == 0) return 1;
   return 0;
}

int calculate_S (char *a)
{
   if (strcmp(a,"S") == 0) return 1;
   if (strcmp(a,"SorH") == 0) return 1;
   if (strcmp(a,"SandH") == 0) return 1;
   return 0;
}

int calculate_relative_error_on_S (char *a)
{
   if (strcmp(a,"S") == 0) return 1;
   if (strcmp(a,"SorH") == 0) return 1;
   if (strcmp(a,"SandH") == 0) return 1;
   return 0;
}

int calculate_absolute_error_on_S (char *a)
{
   return 0;
}

int is_converged (struct projectData *projData, double maxRelativeError, double maxAbsoluteError)
{
   if (strcmp(projData->refinement_variable,"S") == 0) {
      if (maxRelativeError >= 0 && maxRelativeError <= projData->refinement_relative_tolerance) return 1;
   }

   if (strcmp(projData->refinement_variable,"SorH") == 0) {
      if (maxRelativeError >= 0 && maxRelativeError <= projData->refinement_relative_tolerance) return 1;
      if (maxAbsoluteError >= 0 && maxAbsoluteError <= projData->refinement_absolute_tolerance) return 1;
   }

   if (strcmp(projData->refinement_variable,"SandH") == 0) {
      int passedRel=0;
      int passedAbs=0;
      if (maxRelativeError >= 0 && maxRelativeError <= projData->refinement_relative_tolerance) passedRel=1;
      if (maxAbsoluteError >= 0 && maxAbsoluteError <= projData->refinement_absolute_tolerance) passedAbs=1;
      if (passedRel && passedAbs) return 1;
   }

   return 0;
}

void add_inputFrequencyPlan (struct projectData *data, int type, double frequency, double start, double stop, double step, int pointsPerDecade, int lineNumber, int refine)
{
   // allocate more plans, if needed
   if (data->inputFrequencyPlansCount == data->inputFrequencyPlansAllocated) {
      data->inputFrequencyPlansAllocated+=5;
      data->inputFrequencyPlans=(struct inputFrequencyPlan *) realloc (data->inputFrequencyPlans,data->inputFrequencyPlansAllocated*sizeof(struct inputFrequencyPlan));
   }

   // assign the next plan and update the counter
   if (type == 0) {
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].type=type;
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].frequency=-1;
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].start=start;
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].stop=stop;
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].step=step;
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].pointsPerDecade=-1;
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].lineNumber=lineNumber;
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].refine=refine;
   } else if (type == 1) {
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].type=type;
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].frequency=-1;
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].start=start;
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].stop=stop;
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].step=-1;
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].pointsPerDecade=pointsPerDecade;
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].lineNumber=lineNumber;
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].refine=refine;
   } else if (type == 2) {
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].type=type;
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].frequency=frequency;
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].start=-1;
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].stop=-1;
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].step=-1;
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].pointsPerDecade=-1;
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].lineNumber=lineNumber;
      data->inputFrequencyPlans[data->inputFrequencyPlansCount].refine=refine;
   }
   data->inputFrequencyPlansCount++;
}

void init_project (struct projectData *data) {

   data->version_name=allocCopyString("#OpenParEM3Dproject");
   data->version_value=allocCopyString("1.0");

   data->project_name=allocCopyString("");
   data->project_calculate_poynting=0;
   data->project_save_fields=0;

   data->mesh_file=allocCopyString("");
   data->mesh_order=1;
   data->mesh_save_refined=0;
   data->mesh_2D_refinement_fraction=0.025;
   data->mesh_3D_refinement_fraction=0.005;
   data->mesh_quality_limit=20;

   data->port_definition_file=allocCopyString("");

   data->refinement_frequency=allocCopyString("highlow");
   data->refinement_iteration_min=1;
   data->refinement_iteration_max=10;
   data->refinement_required_passes=3;
   data->refinement_relative_tolerance=0.02;
   data->refinement_absolute_tolerance=1e-6;
   data->refinement_variable=allocCopyString("SandH");  

   data->materials_global_path=allocCopyString("../");
   data->materials_global_name=allocCopyString("global_materials");
   data->materials_local_path=allocCopyString("./");
   data->materials_local_name=allocCopyString("local_materials");
   data->materials_default_boundary=allocCopyString("PEC");
   data->materials_check_limits=1;

   data->inputFrequencyPlansAllocated=5;
   data->inputFrequencyPlansCount=0;
   data->inputFrequencyPlans=(struct inputFrequencyPlan *)malloc(data->inputFrequencyPlansAllocated*sizeof(struct inputFrequencyPlan));

   data->reference_impedance=50;
   data->touchstone_frequency_unit=allocCopyString("GHz");
   data->touchstone_format=allocCopyString("DB");

   data->solution_modes=0;                                  // This is an OpenParEM2D project variable - set to an invalid number
   data->solution_temperature=25; 
   data->solution_2D_tolerance=1e-13;
   data->solution_3D_tolerance=1e-13;
   data->solution_iteration_limit=5000;
   data->solution_modes_buffer=0;
   data->solution_impedance_definition=allocCopyString("invalid");   // OpenParEM2D parameter
   data->solution_impedance_calculation=allocCopyString("invalid");  // OpenParEM2D parameter
   data->solution_check_closed_loop=1;
   data->solution_check_homogeneous=1;
   data->solution_accurate_residual=0;
   data->solution_shift_invert=1;
   data->solution_use_initial_guess=1;
   data->solution_shift_factor=1;

   data->output_show_refining_mesh=0;
   data->output_show_postprocessing=0;
   data->output_show_iterations=0;
   data->output_show_license=0;

   data->test_create_cases=0;
   data->test_show_detailed_cases=0;

   data->debug_show_memory=0;
   data->debug_show_project=0;
   data->debug_show_frequency_plan=0;
   data->debug_show_materials=0;
   data->debug_show_port_definitions=0;
   data->debug_show_impedance_details=0;
   data->debug_save_port_fields=0;
   data->debug_skip_mixed_conversion=0;
   data->debug_skip_forced_reciprocity=0;
   data->debug_tempfiles_keep=0;

   data->field_points_count=0;
   data->field_points_allocated=0;
   data->field_points_x=NULL;
   data->field_points_y=NULL;
   data->field_points_z=NULL;
}

void free_project (struct projectData *data) {
   if (data == NULL) return;

   if (data->version_name) free (data->version_name);
   if (data->version_value) free (data->version_value);
   if (data->project_name) free (data->project_name);
   if (data->mesh_file) free (data->mesh_file);
   if (data->port_definition_file) free (data->port_definition_file);
   if (data->materials_global_path) free (data->materials_global_path);
   if (data->materials_global_name) free (data->materials_global_name);
   if (data->materials_local_path) free (data->materials_local_path);
   if (data->materials_local_name) free (data->materials_local_name);
   if (data->materials_default_boundary) free (data->materials_default_boundary);
   if (data->refinement_variable) free (data->refinement_variable);
   if (data->touchstone_frequency_unit) free (data->touchstone_frequency_unit);
   if (data->touchstone_format) free (data->touchstone_format);
   if (data->solution_impedance_definition) free (data->solution_impedance_definition);
   if (data->solution_impedance_calculation) free (data->solution_impedance_calculation);
   if (data->inputFrequencyPlans) free(data->inputFrequencyPlans);
   if (data->field_points_x) free(data->field_points_x);
   if (data->field_points_y) free(data->field_points_y);
   if (data->field_points_z) free(data->field_points_z);
}

void print_project (struct projectData *data, struct projectData *defaultData, const char *indent) {
   int i;
   int matched=0;
   char* logic[2];
   char* comment[2];

   if (data == NULL) return;


   logic[0]=allocCopyString("false");
   logic[1]=allocCopyString("true");

   comment[0]=allocCopyString("");
   comment[1]=allocCopyString("//");

   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%s %s\n",indent,data->version_name,data->version_value);
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s//Commented lines show either unspecified inputs that utilize the default values or specified inputs that match the default values.\n",indent);

   matched=0;  if (defaultData && data->project_calculate_poynting == defaultData->project_calculate_poynting) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sproject.calculate.poynting %s\n",indent,comment[matched],logic[data->project_calculate_poynting]);

   matched=0;  if (defaultData && data->project_save_fields == defaultData->project_save_fields) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sproject.save.fields %s\n",indent,comment[matched],logic[data->project_save_fields]);

   matched=0; if (defaultData && strcmp(data->mesh_file,defaultData->mesh_file) == 0) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%smesh.file %s\n",indent,comment[matched],data->mesh_file);

   matched=0; if (defaultData && data->mesh_order == defaultData->mesh_order) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%smesh.order %d\n",indent,comment[matched],data->mesh_order);

   matched=0;  if (defaultData && data->mesh_save_refined == defaultData->mesh_save_refined) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%smesh.save.refined %s\n",indent,comment[matched],logic[data->mesh_save_refined]);

   matched=0; if (defaultData && double_compare(data->mesh_3D_refinement_fraction,defaultData->mesh_3D_refinement_fraction,1e-14)) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%smesh.refinement.fraction %g\n",indent,comment[matched],data->mesh_3D_refinement_fraction);

   matched=0; if (defaultData && double_compare(data->mesh_quality_limit,defaultData->mesh_quality_limit,1e-14)) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%smesh.quality.limit %g\n",indent,comment[matched],data->mesh_quality_limit);

   matched=0; if (defaultData && strcmp(data->port_definition_file,defaultData->port_definition_file) == 0) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%smode.definition.file %s\n",indent,comment[matched],data->port_definition_file);

   matched=0; if (defaultData && strcmp(data->refinement_frequency,defaultData->refinement_frequency) == 0) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%srefinement.frequency %s\n",indent,comment[matched],data->refinement_frequency);

   matched=0; if (defaultData && data->refinement_iteration_min == defaultData->refinement_iteration_min) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%srefinement.iteration.min %d\n",indent,comment[matched],data->refinement_iteration_min);

   matched=0; if (defaultData && data->refinement_iteration_max == defaultData->refinement_iteration_max) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%srefinement.iteration.max %d\n",indent,comment[matched],data->refinement_iteration_max);

   matched=0; if (defaultData && data->refinement_required_passes == defaultData->refinement_required_passes) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%srefinement.required.passes %d\n",indent,comment[matched],data->refinement_required_passes);

   matched=0; if (defaultData && double_compare(data->refinement_relative_tolerance,defaultData->refinement_relative_tolerance,1e-14)) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%srefinement.relative.tolerance %g\n",indent,comment[matched],data->refinement_relative_tolerance);

   matched=0; if (defaultData && double_compare(data->refinement_absolute_tolerance,defaultData->refinement_absolute_tolerance,1e-14)) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%srefinement.absolute.tolerance %g\n",indent,comment[matched],data->refinement_absolute_tolerance);

   matched=0; if (defaultData && strcmp(data->refinement_variable,defaultData->refinement_variable) == 0) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%srefinement.variable %s\n",indent,comment[matched],data->refinement_variable);

   matched=0; if (defaultData && strcmp(data->materials_global_path,defaultData->materials_global_path) == 0) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%smaterials.global.path %s\n",indent,data->materials_global_path);

   matched=0; if (defaultData && strcmp(data->materials_global_name,defaultData->materials_global_name) == 0) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%smaterials.global.name %s\n",indent,data->materials_global_name);

   matched=0; if (defaultData && strcmp(data->materials_local_path,defaultData->materials_local_path) == 0) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%smaterials.local.path %s\n",indent,data->materials_local_path);

   matched=0; if (defaultData && strcmp(data->materials_local_name,defaultData->materials_local_name) == 0) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%smaterials.local.name %s\n",indent,data->materials_local_name);

   matched=0; if (defaultData && strcmp(data->materials_default_boundary,defaultData->materials_default_boundary) == 0) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%smaterials.default.boundary %s\n",indent,data->materials_default_boundary);

   matched=0;  if (defaultData && data->materials_check_limits == defaultData->materials_check_limits) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%smaterials.check.limits %s\n",indent,comment[matched],logic[data->materials_check_limits]);

   // no default frequency plans, so print all
   i=0;
   while (i < data->inputFrequencyPlansCount) {
      if (data->inputFrequencyPlans[i].type == 0) {
         if (data->inputFrequencyPlans[i].refine == 0) {
            prefix(); PetscPrintf(PETSC_COMM_WORLD,"%sfrequency.plan.linear %g,%g,%g\n",indent,
                        data->inputFrequencyPlans[i].start,data->inputFrequencyPlans[i].stop,data->inputFrequencyPlans[i].step);
         } else {
            prefix(); PetscPrintf(PETSC_COMM_WORLD,"%sfrequency.plan.linear.refine %g,%g,%g\n",indent,
                        data->inputFrequencyPlans[i].start,data->inputFrequencyPlans[i].stop,data->inputFrequencyPlans[i].step);
         }
      } else if (data->inputFrequencyPlans[i].type == 1) {
         if (data->inputFrequencyPlans[i].refine == 0) {
            prefix(); PetscPrintf(PETSC_COMM_WORLD,"%sfrequency.plan.log %g,%g,%d\n",indent,
                        data->inputFrequencyPlans[i].start,data->inputFrequencyPlans[i].stop,data->inputFrequencyPlans[i].pointsPerDecade);
         } else {
            prefix(); PetscPrintf(PETSC_COMM_WORLD,"%sfrequency.plan.log.refine %g,%g,%d\n",indent,
                        data->inputFrequencyPlans[i].start,data->inputFrequencyPlans[i].stop,data->inputFrequencyPlans[i].pointsPerDecade);
         }
      } else if (data->inputFrequencyPlans[i].type == 2) {
         if (data->inputFrequencyPlans[i].refine == 0) {
            prefix(); PetscPrintf(PETSC_COMM_WORLD,"%sfrequency.plan.point %g\n",indent,data->inputFrequencyPlans[i].frequency);
         } else {
            prefix(); PetscPrintf(PETSC_COMM_WORLD,"%sfrequency.plan.point.refine %g\n",indent,data->inputFrequencyPlans[i].frequency);
         }
      }
      i++;
   }

   matched=0; if (defaultData && double_compare(data->reference_impedance,defaultData->reference_impedance,1e-14)) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sreference.impedance %g\n",indent,comment[matched],data->reference_impedance);

   matched=0; if (defaultData && strcmp(data->touchstone_frequency_unit,defaultData->touchstone_frequency_unit) == 0) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%stouchstone.frequency.unit %s\n",indent,comment[matched],data->touchstone_frequency_unit);

   matched=0; if (defaultData && strcmp(data->touchstone_format,defaultData->touchstone_format) == 0) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%stouchstone.format %s\n",indent,comment[matched],data->touchstone_format);

   matched=0; if (defaultData && double_compare(data->solution_temperature,defaultData->solution_temperature,1e-14)) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%ssolution.temperature %g\n",indent,comment[matched],data->solution_temperature);

   matched=0; if (defaultData && double_compare(data->solution_2D_tolerance,defaultData->solution_2D_tolerance,1e-14)) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%ssolution.2D.tolerance %g\n",indent,comment[matched],data->solution_2D_tolerance);

   matched=0; if (defaultData && double_compare(data->solution_3D_tolerance,defaultData->solution_3D_tolerance,1e-14)) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%ssolution.3D.tolerance %g\n",indent,comment[matched],data->solution_3D_tolerance);

   matched=0; if (defaultData && data->solution_iteration_limit == defaultData->solution_iteration_limit) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%ssolution.iteration.limit %d\n",indent,comment[matched],data->solution_iteration_limit);

   matched=0; if (defaultData && data->solution_modes_buffer == defaultData->solution_modes_buffer) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%ssolution.modes.buffer %d\n",indent,comment[matched],data->solution_modes_buffer);

   matched=0;  if (defaultData && data->solution_check_closed_loop == defaultData->solution_check_closed_loop) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%ssolution.check.closed.loop %s\n",indent,comment[matched],logic[data->solution_check_closed_loop]);

   matched=0;  if (defaultData && data->solution_check_homogeneous == defaultData->solution_check_homogeneous) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%ssolution.check.homogeneous %s\n",indent,comment[matched],logic[data->solution_check_homogeneous]);

   matched=0;  if (defaultData && data->solution_accurate_residual == defaultData->solution_accurate_residual) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%ssolution.accurate.residual %s\n",indent,comment[matched],logic[data->solution_accurate_residual]);

   matched=0;  if (defaultData && data->solution_shift_invert == defaultData->solution_shift_invert) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%ssolution.shift.invert %s\n",indent,comment[matched],logic[data->solution_shift_invert]);

   matched=0;  if (defaultData && data->solution_use_initial_guess == defaultData->solution_use_initial_guess) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%ssolution.use.initial.guess %s\n",indent,comment[matched],logic[data->solution_use_initial_guess]);

   matched=0; if (defaultData && double_compare(data->solution_shift_factor,defaultData->solution_shift_factor,1e-14)) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%ssolution.shift.factor %g\n",indent,comment[matched],data->solution_shift_factor);

   matched=0;  if (defaultData && data->output_show_refining_mesh == defaultData->output_show_refining_mesh) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%soutput.show.refining.mesh %s\n",indent,comment[matched],logic[data->output_show_refining_mesh]);

   matched=0;  if (defaultData && data->output_show_postprocessing == defaultData->output_show_postprocessing) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%soutput.show.postprocessing %s\n",indent,comment[matched],logic[data->output_show_postprocessing]);

   matched=0;  if (defaultData && data->output_show_iterations == defaultData->output_show_iterations) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%soutput.show.iterations %s\n",indent,comment[matched],logic[data->output_show_iterations]);

   matched=0;  if (defaultData && data->output_show_license == defaultData->output_show_license) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%soutput.show.license %s\n",indent,comment[matched],logic[data->output_show_license]);

   matched=0;  if (defaultData && data->test_create_cases == defaultData->test_create_cases) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%stest.create.cases %s\n",indent,comment[matched],logic[data->test_create_cases]);

   matched=0;  if (defaultData && data->test_show_detailed_cases == defaultData->test_show_detailed_cases) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%stest.show.detailed.cases %s\n",indent,comment[matched],logic[data->test_show_detailed_cases]);

   matched=0;  if (defaultData && data->debug_show_memory == defaultData->debug_show_memory) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sdebug.show.memory %s\n",indent,comment[matched],logic[data->debug_show_memory]);

   matched=0;  if (defaultData && data->debug_show_project == defaultData->debug_show_project) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sdebug.show.project %s\n",indent,comment[matched],logic[data->debug_show_project]);

   matched=0;  if (defaultData && data->debug_show_frequency_plan == defaultData->debug_show_frequency_plan) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sdebug.show.frequency.plan %s\n",indent,comment[matched],logic[data->debug_show_frequency_plan]);

   matched=0;  if (defaultData && data->debug_show_materials == defaultData->debug_show_materials) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sdebug.show.materials %s\n",indent,comment[matched],logic[data->debug_show_materials]);

   matched=0;  if (defaultData && data->debug_show_port_definitions == defaultData->debug_show_port_definitions) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sdebug.show.port.definitions %s\n",indent,comment[matched],logic[data->debug_show_port_definitions]);

   matched=0;  if (defaultData && data->debug_show_impedance_details == defaultData->debug_show_impedance_details) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sdebug.show.impedance.details %s\n",indent,comment[matched],logic[data->debug_show_impedance_details]);

   matched=0;  if (defaultData && data->debug_save_port_fields == defaultData->debug_save_port_fields) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sdebug.save.port.fields %s\n",indent,comment[matched],logic[data->debug_save_port_fields]);

   matched=0;  if (defaultData && data->debug_skip_mixed_conversion == defaultData->debug_skip_mixed_conversion) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sdebug.skip.mixed.conversion %s\n",indent,comment[matched],logic[data->debug_skip_mixed_conversion]);

   matched=0;  if (defaultData && data->debug_skip_forced_reciprocity == defaultData->debug_skip_forced_reciprocity) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sdebug.skip.forced.reciprocity %s\n",indent,comment[matched],logic[data->debug_skip_forced_reciprocity]);

   matched=0;  if (defaultData && data->debug_tempfiles_keep == defaultData->debug_tempfiles_keep) matched=1;
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sdebug.tempfiles.keep %s\n",indent,comment[matched],logic[data->debug_tempfiles_keep]);

   // no default field points, so print all
   i=0;
   while (i < data->field_points_count) {
      prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sfield.point %g,%g,%g\n",indent,indent,data->field_points_x[i],data->field_points_y[i],data->field_points_z[i]);
      i++;
   }

   free(logic[0]);
   free(logic[1]);
   free(comment[0]);
   free(comment[1]);
}

int get_bool (char *a) {
   if (is_true(a)) return 1;
   return 0;
}

void print_invalid_entry (PetscErrorCode *error, int lineNumber, const char *indent) {
   prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3133: Invalid entry at line %d.\n",indent,indent,lineNumber);
   *error=1;
   return;
}

char* get_project_name (const char *filename) {
   char *b=NULL;
   int i;

   if (filename == NULL) return NULL;

   // copy filename into b (since filename is const)
   b=(char *)malloc((strlen(filename)+1)*sizeof(char));
   if (b == NULL) return NULL;
   strcpy(b,filename);

   // chop off the extension
   i=strlen(b)-1;
   while (i >= 0) {
      if (b[i] == '.') {b[i]='\0'; return b;}
      if (i == 0) break;
      i--;
   }

   return b;
}

int has_refinementFrequencyPlan (struct projectData *data) {
   int i;

   if (data == NULL) return 0;

   i=0;
   while (i < data->inputFrequencyPlansCount) {
      if (data->inputFrequencyPlans[i].refine) return 1;
      i++;
   }
   return 0;
}

PetscErrorCode load_project_file (const char *filename, struct projectData *data, const char* indent) {
   PetscMPIInt size,rank;
   FILE *fp=NULL;
   int openedFile,failedLoad;
   char *line=NULL;
   size_t len=0;
   ssize_t line_size;
   char *keyword=NULL;
   char *value=NULL;
   PetscErrorCode ierr=0;
   int lineCount=0;
   int commaCount;
   int lineIterationMax;
   int lineRefinementFrequency=0;
   int i;
   long unsigned int j;
   double x,y,z;
   int pointsPerDecade;
   double start,stop,step,frequency;
   int length;
   int planCount;
   int planType;
   double planFrequency;
   double planStart;
   double planStop;
   double planStep;
   int planPointsPerDecade;
   int planRefine;
   int planLineNumber;

   if (filename == NULL) return 1;

   MPI_Comm_size(PETSC_COMM_WORLD, &size);
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

   lineIterationMax=-1;

   // open the file
   openedFile=0;
   if (rank == 0) {
      fp=fopen(filename,"r");
      if (fp) openedFile=1;
   }
   ierr=MPI_Bcast (&openedFile,1,MPI_INT,0,PETSC_COMM_WORLD);
   if (! openedFile) {
      prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3173: Failed to open file \"%s\" for reading.\n",indent,indent,filename);
      return 1;
   }

   if (rank == 0) {
      // fp is open

      // check the version first

      line_size=getline(&line,&len,fp);
      while (line_size >= 0) {
         lineCount++;

         if (! is_blank(line) && len > 0) {
            line=removeNewLineChar(line);
            line=removeComment(line);

            keyword=strtok(line," ");

            if (keyword != NULL) {
               if (strcmp(keyword,data->version_name) == 0) {
                  value=strtok(NULL," ");
                  if (is_text(value)) {
                     if (strcmp(value,data->version_value) != 0) {
                        prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3134: Version mismatch. Expecting on first line: %s %s\n",
                                    indent,indent,data->version_name,data->version_value);
                        ierr=1;
                     }
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               } else {
                  prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3135: Missing version. Expecting on first line: %s %s\n",
                              indent,indent,data->version_name,data->version_value);
                  ierr=1;
               }
            }
            break;
         }
         if (line) {free(line); line=NULL;}
         line_size=getline(&line,&len,fp);
      }
   }

   // stop if an ierr is found so far
   failedLoad=0;
   if (rank == 0) failedLoad=ierr;
   ierr=MPI_Bcast (&failedLoad,1,MPI_INT,0,PETSC_COMM_WORLD);
   if (failedLoad) {
      if (rank == 0) fclose (fp);
      prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3165: Failed to load project.\n",indent,indent);
      return 1;
   }

   // continue
   if (rank == 0) {
      // fp is open

      // get the project name from the filename
      free(data->project_name);
      data->project_name=get_project_name (filename);

      // check everything else
      if (line) {free(line); line=NULL;}
      line_size=getline(&line,&len,fp);
      while (line_size >= 0) {
         lineCount++;

         if (! is_blank(line) && len > 0) {
            line=removeNewLineChar(line);
            line=removeComment(line);
            commaCount=comma_count(line);

            keyword=strtok(line," ");

            if (keyword != NULL) {

               if (strcmp(keyword,"project.calculate.poynting") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->project_calculate_poynting=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"project.save.fields") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->project_save_fields=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"mesh.file") == 0) {
                  value=strtok(NULL," ");
                  if (is_text(value)) {
                     free(data->mesh_file);
                     data->mesh_file=allocCopyString(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent); 
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"mesh.order") == 0) {
                  value=strtok(NULL," ");
                  if (is_int(value)) {
                     data->mesh_order=atoi(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                     if (data->mesh_order < 1) {
                        ierr=1;
                        prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3136: Value must be >= 1 at line %d.\n",indent,indent,lineCount);
                     }
                     if (data->mesh_order > 20) {
                        ierr=1;
                        prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3137: Value must be <= 20 at line %d.\n",indent,indent,lineCount);
                     }
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"mesh.save.refined") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->mesh_save_refined=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"mesh.refinement.fraction") == 0) {
                  value=strtok(NULL," ");
                  if (is_double(value)) {
                     data->mesh_3D_refinement_fraction=atof(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                     if (data->mesh_3D_refinement_fraction <= 0) {
                        ierr=1;
                        prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3143: Value must be > 0 at line %d.\n",indent,indent,lineCount);
                     }
                     if (data->mesh_3D_refinement_fraction > 1) {
                        ierr=1;
                        prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3144: Value must be <= 1 at line %d.\n",indent,indent,lineCount);
                     }
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"mesh.quality.limit") == 0) {
                  value=strtok(NULL," ");
                  if (is_double(value)) {
                     data->mesh_quality_limit=atof(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                     if (data->mesh_quality_limit <= 0) {
                        ierr=1;
                        prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3145: Value must be > 0 at line %d.\n",indent,indent,lineCount);
                     }
                     if (data->mesh_quality_limit > 1) {
                        prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sWarning: Value at line %d is large.\n",indent,indent,lineCount);
                     }
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"port.definition.file") == 0) {
                  value=strtok(NULL," ");
                  if (is_text(value)) {
                     free(data->port_definition_file);
                     data->port_definition_file=allocCopyString(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  }
               }

               else if (strcmp(keyword,"refinement.frequency") == 0) {
                  value=strtok(NULL," ");
                  if (is_text(value)) {
                     free(data->refinement_frequency);
                     data->refinement_frequency=allocCopyString(value);
                     lineRefinementFrequency=lineCount;
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else {
                     free(data->refinement_frequency);
                     data->refinement_frequency=allocCopyString("highlow");
                     ierr=1;
                     prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3146: Value must be \"all\", \"none\", \"high\", \"low\", \"highlow (default)\", \"lowhigh\", or \"plan\" at line %d.\n",
                                                            indent,indent,lineCount);
                  }
                  if (! is_refinement_frequency(data->refinement_frequency)) {
                     ierr=1;
                     prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3147: Value must be \"all\", \"none\", \"high\", \"low\", \"highlow (default)\", \"lowhigh\", or \"plan\" at line %d.\n",
                                                            indent,indent,lineCount);
                  }
               }

               else if (strcmp(keyword,"refinement.required.passes") == 0) {
                  value=strtok(NULL," ");
                  if (is_int(value)) {
                     data->refinement_required_passes=atoi(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                     if (data->refinement_required_passes < 1) {
                        ierr=1;
                        prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3148: Value must be >= 1 at line %d.\n",indent,indent,lineCount);
                     }
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"refinement.relative.tolerance") == 0) {
                  value=strtok(NULL," ");
                  if (is_double(value)) {
                     data->refinement_relative_tolerance=atof(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                     if (data->refinement_relative_tolerance <= 0) {
                        ierr=1;
                        prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3149: Value must be > 0 at line %d.\n",indent,indent,lineCount);
                     }
                     if (data->refinement_relative_tolerance > 1) {
                        ierr=1;
                        prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3150: Value must be <= 1 at line %d.\n",indent,indent,lineCount);
                     }
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"refinement.absolute.tolerance") == 0) {
                  value=strtok(NULL," ");
                  if (is_double(value)) {
                     data->refinement_absolute_tolerance=atof(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                     if (data->refinement_absolute_tolerance <= 0) {
                        ierr=1;
                        prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3056: Value must be > 0 at line %d.\n",indent,indent,lineCount);
                     }
                     if (data->refinement_absolute_tolerance > 10) {
                        ierr=1;
                        prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3139: Value must be <= 10 at line %d.\n",indent,indent,lineCount);
                     }
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"refinement.variable") == 0) {
                  value=strtok(NULL," ");
                  if (is_text(value)) {
                     free(data->refinement_variable);
                     data->refinement_variable=allocCopyString(value);
                     lineRefinementFrequency=lineCount;
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else {
                     free(data->refinement_variable);
                     data->refinement_variable=allocCopyString("E");
                     ierr=1;
                     prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3176: Value must be \"S\", \"SorH\", or \"SandH\" at line %d.\n",
                                                            indent,indent,lineCount);
                  }
                  if (! is_refinement_variable(data->refinement_variable)) {
                     ierr=1;
                     prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3177: Value must be \"S\", \"SorH\", or \"SandH\" at line %d.\n",
                                                            indent,indent,lineCount);
                  }
               }

               else if (strcmp(keyword,"refinement.iteration.min") == 0) {
                  value=strtok(NULL," ");
                  if (is_int(value)) {
                     data->refinement_iteration_min=atoi(value);
                     if (data->refinement_iteration_min < 1) {
                        ierr=1;
                        prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3151: Value must be >= 1 at line %d.\n",indent,indent,lineCount);
                     }
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"refinement.iteration.max") == 0) {
                  value=strtok(NULL," ");
                  if (is_int(value)) {
                     data->refinement_iteration_max=atoi(value);
                     lineIterationMax=lineCount;
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"materials.global.path") == 0) {
                  value=strtok(NULL," ");
                  if (is_text(value)) {
                     free(data->materials_global_path);
                     data->materials_global_path=allocCopyString(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else {
                     free(data->materials_global_path);
                     data->materials_global_path=allocCopyString("./");
                  }
               }

               else if (strcmp(keyword,"materials.global.name") == 0) {
                  value=strtok(NULL," ");
                  if (is_text(value)) {
                     free(data->materials_global_name);
                     data->materials_global_name=allocCopyString(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else {
                     free(data->materials_global_name);
                     data->materials_global_name=allocCopyString("");
                  }
               }

               else if (strcmp(keyword,"materials.local.path") == 0) {
                  value=strtok(NULL," ");
                  if (is_text(value)) {
                     free(data->materials_local_path);
                     data->materials_local_path=allocCopyString(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else {
                     free(data->materials_local_path);
                     data->materials_local_path=allocCopyString("./");
                  }
               }

               else if (strcmp(keyword,"materials.local.name") == 0) {
                  value=strtok(NULL," ");
                  if (is_text(value)) {
                     free(data->materials_local_name);
                     data->materials_local_name=allocCopyString(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else {
                     free(data->materials_local_name);
                     data->materials_local_name=allocCopyString("");
                  }
               }

               else if (strcmp(keyword,"materials.default.boundary") == 0) {
                  value=strtok(NULL," ");
                  if (is_text(value)) {
                     free(data->materials_default_boundary);
                     data->materials_default_boundary=allocCopyString(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else {
                     free(data->materials_default_boundary);
                     data->materials_default_boundary=allocCopyString("PEC");
                  }
               }

               else if (strcmp(keyword,"materials.check.limits") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->materials_check_limits=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"frequency.plan.linear") == 0 || strcmp(keyword,"frequency.plan.linear.refine") == 0) {
                  value=strtok(NULL,",");
                  if (is_double(value)) {
                     start=atof(value);
                     value=strtok(NULL,",");
                     if (is_double(value)) {
                        stop=atof(value);
                        value=strtok(NULL,",");
                        if (is_double(value)) {
                           step=atof(value);
                           if (strcmp(keyword,"frequency.plan.linear.refine") == 0) {
                              add_inputFrequencyPlan (data,0,-1,start,stop,step,-1,lineCount,1);
                           } else {
                              add_inputFrequencyPlan (data,0,-1,start,stop,step,-1,lineCount,0);
                           }
                           value=strtok(NULL,",");
                           if (value) print_invalid_entry (&ierr,lineCount,indent);
                        } else print_invalid_entry (&ierr,lineCount,indent);
                     } else print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"frequency.plan.log") == 0 || strcmp(keyword,"frequency.plan.log.refine") == 0) {
                  value=strtok(NULL,",");
                  if (is_double(value)) {
                     start=atof(value);
                     value=strtok(NULL,",");
                     if (is_double(value)) {
                        stop=atof(value);
                        value=strtok(NULL,",");
                        if (is_int(value)) {
                           pointsPerDecade=atoi(value);
                           if (strcmp(keyword,"frequency.plan.log.refine") == 0) {
                              add_inputFrequencyPlan (data,1,-1,start,stop,-1,pointsPerDecade,lineCount,1);
                           } else {
                              add_inputFrequencyPlan (data,1,-1,start,stop,-1,pointsPerDecade,lineCount,0);
                           }
                           value=strtok(NULL,",");
                           if (value) print_invalid_entry (&ierr,lineCount,indent);
                        } else print_invalid_entry (&ierr,lineCount,indent);
                     } else print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"frequency.plan.point") == 0 || strcmp(keyword,"frequency.plan.point.refine") == 0) {
                  value=strtok(NULL,",");
                  if (is_double(value)) {
                     frequency=atof(value);
                     if (strcmp(keyword,"frequency.plan.point.refine") == 0) {
                        add_inputFrequencyPlan (data,2,frequency,-1,-1,-1,-1,lineCount,1);
                     } else {
                        add_inputFrequencyPlan (data,2,frequency,-1,-1,-1,-1,lineCount,0);
                     }
                     value=strtok(NULL,",");
                     if (value) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"reference.impedance") == 0) {
                  value=strtok(NULL," ");
                  if (is_double(value)) {
                     data->reference_impedance=atof(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                     if (data->reference_impedance < 0) {
                        ierr=1;
                        prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3153: Value must be >= 0 at line %d.\n",indent,indent,lineCount);
                     }
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"touchstone.frequency.unit") == 0) {
                  value=strtok(NULL," ");
                  if (is_text(value)) {
                     free(data->touchstone_frequency_unit);
                     data->touchstone_frequency_unit=allocCopyString(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                     if (!(strcmp(data->touchstone_frequency_unit,"Hz") == 0) && 
                         !(strcmp(data->touchstone_frequency_unit,"kHz") == 0) &&
                         !(strcmp(data->touchstone_frequency_unit,"MHz") == 0) &&
                         !(strcmp(data->touchstone_frequency_unit,"GHz") == 0)) {
                        ierr=1;
                        prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3155: Value must be \"Hz\", \"kHz\", \"MHz\", or \"GHz\" at line %d.\n",indent,indent,lineCount);
                     }
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"touchstone.format") == 0) {
                  value=strtok(NULL," ");
                  if (is_text(value)) {
                     free(data->touchstone_format);
                     data->touchstone_format=allocCopyString(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                     if (!(strcmp(data->touchstone_format,"DB") == 0) && 
                         !(strcmp(data->touchstone_format,"MA") == 0) &&
                         !(strcmp(data->touchstone_format,"RI") == 0)) {
                        ierr=1;
                        prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3156: Value must be \"DB\", \"MA\", or \"RI\" at line %d.\n",indent,indent,lineCount);
                     }
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"solution.temperature") == 0) {
                  value=strtok(NULL," ");
                  if (is_double(value)) {
                     data->solution_temperature=atof(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                     if (data->solution_temperature < -273.15) {
                        ierr=1;
                        prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3157: Value must be >= -273.15 at line %d.\n",indent,indent,lineCount);
                     } 
                     if (data->solution_temperature > 1e5) {
                        ierr=1;
                        prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3158: Value must be < 1e5 at line %d.\n",indent,indent,lineCount);
                     } 
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"solution.2D.tolerance") == 0) {
                  value=strtok(NULL," ");
                  if (is_double(value)) {
                     data->solution_2D_tolerance=atof(value);
                     //lineSolutionTolerance=lineCount;
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                     if (data->solution_2D_tolerance <= 0) {
                        ierr=1;
                        prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3159: Value must be > 0 at line %d.\n",indent,indent,lineCount);
                     } 
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"solution.3D.tolerance") == 0) {
                  value=strtok(NULL," ");
                  if (is_double(value)) {
                     data->solution_3D_tolerance=atof(value);
                     //lineSolutionTolerance=lineCount;
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                     if (data->solution_3D_tolerance <= 0) {
                        ierr=1;
                        prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3160: Value must be > 0 at line %d.\n",indent,indent,lineCount);
                     }
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"solution.iteration.limit") == 0) {
                  value=strtok(NULL," ");
                  if (is_int(value)) {
                     data->solution_iteration_limit=atoi(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                     if (data->solution_iteration_limit < 1) {
                        ierr=1;
                        prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3161: Value must be >= 1 at line %d.\n",indent,indent,lineCount);
                     } 
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"solution.modes.buffer") == 0) {
                  value=strtok(NULL," ");
                  if (is_int(value)) {
                     data->solution_modes_buffer=atoi(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                     if (data->solution_modes_buffer < 0) {
                        ierr=1;
                        prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3162: Value must be >= 0 at line %d.\n",indent,indent,lineCount);
                     }
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"solution.check.closed.loop") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->solution_check_closed_loop=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"solution.check.homogeneous") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->solution_check_homogeneous=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"solution.accurate.residual") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->solution_accurate_residual=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"solution.shift.invert") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->solution_shift_invert=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"solution.use.initial.guess") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->solution_use_initial_guess=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"solution.shift.factor") == 0) {
                  value=strtok(NULL," ");
                  if (is_double(value)) {
                     data->solution_shift_factor=atof(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                     if (data->solution_shift_factor < 1) {
                        ierr=1;
                        prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3163: Value must be >= 1 at line %d.\n",indent,indent,lineCount);
                     }
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"output.show.refining.mesh") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->output_show_refining_mesh=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"output.show.postprocessing") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->output_show_postprocessing=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"output.show.iterations") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->output_show_iterations=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"output.show.license") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->output_show_license=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"test.create.cases") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->test_create_cases=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"test.show.detailed.cases") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->test_show_detailed_cases=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"debug.show.memory") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->debug_show_memory=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"debug.show.project") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->debug_show_project=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"debug.show.frequency.plan") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->debug_show_frequency_plan=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"debug.show.materials") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->debug_show_materials=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"debug.show.port.definitions") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->debug_show_port_definitions=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"debug.show.impedance.details") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->debug_show_impedance_details=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"debug.save.port.fields") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->debug_save_port_fields=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"debug.skip.mixed.conversion") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->debug_skip_mixed_conversion=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"debug.skip.forced.reciprocity") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->debug_skip_forced_reciprocity=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"debug.tempfiles.keep") == 0) {
                  value=strtok(NULL," ");
                  if (is_bool(value)) {
                     data->debug_tempfiles_keep=get_bool(value);
                     value=strtok(NULL," ");
                     if (is_text(value)) print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else if (strcmp(keyword,"field.point") == 0) {
                  if (commaCount == 2) {
                     value=strtok(NULL," ,");
                     if (is_double(value)) {
                        x=atof(value);
                        value=strtok(NULL," ,");
                        if (is_double(value)) {
                           y=atof(value);
                           value=strtok(NULL," ,");
                           if (is_double(value)) {
                              z=atof(value);
                              if (data->field_points_allocated == 0) {
                                 data->field_points_allocated=256;
                                 data->field_points_x=(double *)malloc(data->field_points_allocated*sizeof(double));
                                 data->field_points_y=(double *)malloc(data->field_points_allocated*sizeof(double));
                                 data->field_points_z=(double *)malloc(data->field_points_allocated*sizeof(double));
                              }
                              if (data->field_points_count == data->field_points_allocated) {
                                 data->field_points_allocated+=256;
                                 data->field_points_x=(double *)realloc(data->field_points_x,data->field_points_allocated*sizeof(double));
                                 data->field_points_y=(double *)realloc(data->field_points_y,data->field_points_allocated*sizeof(double));
                                 data->field_points_z=(double *)realloc(data->field_points_z,data->field_points_allocated*sizeof(double));
                              }
                              data->field_points_x[data->field_points_count]=x;
                              data->field_points_y[data->field_points_count]=y;
                              data->field_points_z[data->field_points_count]=z;
                              data->field_points_count++;
                           } else print_invalid_entry (&ierr,lineCount,indent);
                        } else print_invalid_entry (&ierr,lineCount,indent);
                     } else print_invalid_entry (&ierr,lineCount,indent);
                  } else print_invalid_entry (&ierr,lineCount,indent);
               }

               else {ierr=1; prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3164: Invalid entry at line %d.\n",indent,indent,lineCount);}
            }
         } 

         if (line) {free(line); line=NULL;}
         line_size=getline(&line,&len,fp);
      }

      if (line) {free(line); line=NULL;}
      fclose (fp);
   }

   // stop if an ierr is found so far
   failedLoad=0;
   if (rank == 0) failedLoad=ierr;
   ierr=MPI_Bcast (&failedLoad,1,MPI_INT,0,PETSC_COMM_WORLD);
   if (failedLoad) {
      prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3188: Failed to load project.\n",indent,indent);
      return 1;
   }

   // some consistency checks

   if (rank == 0) {

      if (data->inputFrequencyPlansCount == 0) {
         ierr=1;
         prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3166: At least one frequency plan is required.\n",indent,indent);
      } else {
         i=0;
         while (i < data->inputFrequencyPlansCount) {
            if (data->inputFrequencyPlans[i].type == 0 || data->inputFrequencyPlans[i].type == 1) {
               if (data->inputFrequencyPlans[i].stop < data->inputFrequencyPlans[i].start) {
                  ierr=1;
                  prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3167: The stopping frequency must be >= to the starting frequency at line %d.\n",indent,indent,data->inputFrequencyPlans[i].lineNumber);
               }
            }
            if (data->inputFrequencyPlans[i].type == 0 && data->inputFrequencyPlans[i].step <= 0) {
               ierr=1;
               prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3168: The frequency step must be positive at line %d.\n",indent,indent,data->inputFrequencyPlans[i].lineNumber);
            }
            if (data->inputFrequencyPlans[i].type == 1 && data->inputFrequencyPlans[i].pointsPerDecade < 1) {
               ierr=1;
               prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3169: The number of points per decade must be >= 1 at line %d.\n",indent,indent,data->inputFrequencyPlans[i].lineNumber);
            }
            i++;
         }

         if (strcmp(data->refinement_frequency,"plan") == 0 && ! has_refinementFrequencyPlan(data)) {
            ierr=1;
            prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3170: A refinement plan is called for at line %d but no plan is provided.\n",indent,indent,lineRefinementFrequency);
         }

         if (strcmp(data->refinement_frequency,"plan") != 0 && has_refinementFrequencyPlan(data)) {
            ierr=1;
            prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3171: A refinement plan is not called for at line %d but a refinement plan or plans is provided.\n",indent,indent,lineRefinementFrequency);
         }

      }

      if (data->refinement_iteration_max < data->refinement_iteration_min) {
         ierr=1;//
         prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3172: The maximum iteration limit of %d must be >= to the minimum iteration limit of %d",indent,indent,data->refinement_iteration_max,data->refinement_iteration_min);
         if (lineIterationMax >= 0) {prefix(); PetscPrintf(PETSC_COMM_WORLD," at line %d.\n",lineIterationMax);}
         else {prefix(); PetscPrintf(PETSC_COMM_WORLD,".\n");}
      }

   }

   // stop if an ierr is found so far
   failedLoad=0;
   if (rank == 0) failedLoad=ierr;
   ierr=MPI_Bcast (&failedLoad,1,MPI_INT,0,PETSC_COMM_WORLD);
   if (failedLoad) {
      prefix(); PetscPrintf(PETSC_COMM_WORLD,"%s%sERROR3187: Failed to load project.\n",indent,indent);
      return 1;
   }

   // send to other ranks

   if (rank == 0) {
      i=1;
      while (i < size) {

         ierr=MPI_Send(&(data->project_calculate_poynting),1,MPI_INT,i,1000001,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->project_save_fields),1,MPI_INT,i,1000002,PETSC_COMM_WORLD);

         length=strlen(data->mesh_file);
         ierr=MPI_Send (&length,1,MPI_INT,i,1000003,PETSC_COMM_WORLD);
         ierr=MPI_Send(data->mesh_file,length,MPI_CHAR,i,1000004,PETSC_COMM_WORLD);

         ierr=MPI_Send(&(data->mesh_order),1,MPI_INT,i,1000005,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->mesh_save_refined),1,MPI_INT,i,1000083,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->mesh_2D_refinement_fraction),1,MPI_DOUBLE,i,1000008,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->mesh_3D_refinement_fraction),1,MPI_DOUBLE,i,1000009,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->mesh_quality_limit),1,MPI_DOUBLE,i,1000010,PETSC_COMM_WORLD);

         length=strlen(data->port_definition_file);
         ierr=MPI_Send (&length,1,MPI_INT,i,1000011,PETSC_COMM_WORLD);
         ierr=MPI_Send(data->port_definition_file,length,MPI_CHAR,i,1000012,PETSC_COMM_WORLD);

         length=strlen(data->materials_global_path);
         ierr=MPI_Send (&length,1,MPI_INT,i,1000013,PETSC_COMM_WORLD);
         ierr=MPI_Send(data->materials_global_path,length,MPI_CHAR,i,1000014,PETSC_COMM_WORLD);

         length=strlen(data->materials_global_name);
         ierr=MPI_Send (&length,1,MPI_INT,i,1000015,PETSC_COMM_WORLD);
         ierr=MPI_Send(data->materials_global_name,length,MPI_CHAR,i,1000016,PETSC_COMM_WORLD);

         length=strlen(data->materials_local_path);
         ierr=MPI_Send (&length,1,MPI_INT,i,1000017,PETSC_COMM_WORLD);
         ierr=MPI_Send(data->materials_local_path,length,MPI_CHAR,i,1000018,PETSC_COMM_WORLD);

         length=strlen(data->materials_local_name);
         ierr=MPI_Send (&length,1,MPI_INT,i,1000019,PETSC_COMM_WORLD);
         ierr=MPI_Send(data->materials_local_name,length,MPI_CHAR,i,1000020,PETSC_COMM_WORLD);

         length=strlen(data->materials_default_boundary);
         ierr=MPI_Send (&length,1,MPI_INT,i,1000021,PETSC_COMM_WORLD);
         ierr=MPI_Send(data->materials_default_boundary,length,MPI_CHAR,i,1000022,PETSC_COMM_WORLD);

         ierr=MPI_Send(&(data->materials_check_limits),1,MPI_INT,i,1000023,PETSC_COMM_WORLD);

         length=strlen(data->refinement_frequency);
         ierr=MPI_Send (&length,1,MPI_INT,i,1000024,PETSC_COMM_WORLD);
         ierr=MPI_Send(data->refinement_frequency,length,MPI_CHAR,i,1000025,PETSC_COMM_WORLD);

         length=strlen(data->refinement_variable);
         ierr=MPI_Send (&length,1,MPI_INT,i,1000026,PETSC_COMM_WORLD);
         ierr=MPI_Send(data->refinement_variable,length,MPI_CHAR,i,1000027,PETSC_COMM_WORLD);

         ierr=MPI_Send(&(data->refinement_iteration_min),1,MPI_INT,i,1000028,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->refinement_iteration_max),1,MPI_INT,i,1000029,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->refinement_required_passes),1,MPI_INT,i,1000030,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->refinement_relative_tolerance),1,MPI_DOUBLE,i,1000031,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->refinement_absolute_tolerance),1,MPI_DOUBLE,i,1000084,PETSC_COMM_WORLD);

         ierr=MPI_Send(&(data->inputFrequencyPlansCount),1,MPI_UNSIGNED_LONG,i,1000032,PETSC_COMM_WORLD);

         j=0;
         while (j < data->inputFrequencyPlansCount) {
            ierr=MPI_Send(&(data->inputFrequencyPlans[j].type),1,MPI_INT,i,1000033,PETSC_COMM_WORLD);
            ierr=MPI_Send(&(data->inputFrequencyPlans[j].frequency),1,MPI_DOUBLE,i,1000034,PETSC_COMM_WORLD);
            ierr=MPI_Send(&(data->inputFrequencyPlans[j].start),1,MPI_DOUBLE,i,1000035,PETSC_COMM_WORLD);
            ierr=MPI_Send(&(data->inputFrequencyPlans[j].stop),1,MPI_DOUBLE,i,1000036,PETSC_COMM_WORLD);
            ierr=MPI_Send(&(data->inputFrequencyPlans[j].step),1,MPI_DOUBLE,i,1000037,PETSC_COMM_WORLD);
            ierr=MPI_Send(&(data->inputFrequencyPlans[j].pointsPerDecade),1,MPI_INT,i,1000038,PETSC_COMM_WORLD);
            ierr=MPI_Send(&(data->inputFrequencyPlans[j].refine),1,MPI_INT,i,1000039,PETSC_COMM_WORLD);
            ierr=MPI_Send(&(data->inputFrequencyPlans[j].lineNumber),1,MPI_INT,i,1000040,PETSC_COMM_WORLD);
            j++;
         }

         ierr=MPI_Send(&(data->reference_impedance),1,MPI_DOUBLE,i,1000041,PETSC_COMM_WORLD);

         length=strlen(data->touchstone_frequency_unit);
         ierr=MPI_Send (&length,1,MPI_INT,i,1000044,PETSC_COMM_WORLD);
         ierr=MPI_Send(data->touchstone_frequency_unit,length,MPI_CHAR,i,1000045,PETSC_COMM_WORLD);

         length=strlen(data->touchstone_format);
         ierr=MPI_Send (&length,1,MPI_INT,i,1000046,PETSC_COMM_WORLD);
         ierr=MPI_Send(data->touchstone_format,length,MPI_CHAR,i,1000047,PETSC_COMM_WORLD);

         ierr=MPI_Send(&(data->solution_modes),1,MPI_INT,i,1000048,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->solution_temperature),1,MPI_DOUBLE,i,1000049,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->solution_2D_tolerance),1,MPI_DOUBLE,i,1000050,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->solution_3D_tolerance),1,MPI_DOUBLE,i,1000051,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->solution_iteration_limit),1,MPI_INT,i,1000052,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->solution_modes_buffer),1,MPI_INT,i,1000053,PETSC_COMM_WORLD);

         length=strlen(data->solution_impedance_definition);
         ierr=MPI_Send (&length,1,MPI_INT,i,1000054,PETSC_COMM_WORLD);
         ierr=MPI_Send(data->solution_impedance_definition,length,MPI_CHAR,i,1000055,PETSC_COMM_WORLD);

         length=strlen(data->solution_impedance_calculation);
         ierr=MPI_Send (&length,1,MPI_INT,i,1000056,PETSC_COMM_WORLD);
         ierr=MPI_Send(data->solution_impedance_calculation,length,MPI_CHAR,i,1000057,PETSC_COMM_WORLD);

         ierr=MPI_Send(&(data->solution_check_closed_loop),1,MPI_INT,i,1000058,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->solution_check_homogeneous),1,MPI_INT,i,1000087,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->solution_accurate_residual),1,MPI_INT,i,1000059,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->solution_shift_invert),1,MPI_INT,i,1000060,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->solution_use_initial_guess),1,MPI_INT,i,1000061,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->solution_shift_factor),1,MPI_DOUBLE,i,1000062,PETSC_COMM_WORLD);

         ierr=MPI_Send(&(data->output_show_refining_mesh),1,MPI_INT,i,1000063,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->output_show_postprocessing),1,MPI_INT,i,1000064,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->output_show_iterations),1,MPI_INT,i,1000065,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->output_show_license),1,MPI_INT,i,1000066,PETSC_COMM_WORLD);

         ierr=MPI_Send(&(data->test_create_cases),1,MPI_INT,i,1000067,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->test_show_detailed_cases),1,MPI_INT,i,1000068,PETSC_COMM_WORLD);

         ierr=MPI_Send(&(data->debug_show_memory),1,MPI_INT,i,1000069,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->debug_show_project),1,MPI_INT,i,1000070,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->debug_show_frequency_plan),1,MPI_INT,i,1000071,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->debug_show_materials),1,MPI_INT,i,1000072,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->debug_show_port_definitions),1,MPI_INT,i,1000073,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->debug_show_impedance_details),1,MPI_INT,i,1000074,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->debug_save_port_fields),1,MPI_INT,i,1000075,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->debug_skip_mixed_conversion),1,MPI_INT,i,1000085,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->debug_skip_forced_reciprocity),1,MPI_INT,i,1000086,PETSC_COMM_WORLD);
         ierr=MPI_Send(&(data->debug_tempfiles_keep),1,MPI_INT,i,1000076,PETSC_COMM_WORLD);

         ierr=MPI_Send(&(data->field_points_count),1,MPI_INT,i,1000077,PETSC_COMM_WORLD);
         ierr=MPI_Send(data->field_points_x,data->field_points_count,MPI_DOUBLE,i,1000078,PETSC_COMM_WORLD);
         ierr=MPI_Send(data->field_points_y,data->field_points_count,MPI_DOUBLE,i,1000079,PETSC_COMM_WORLD);
         ierr=MPI_Send(data->field_points_z,data->field_points_count,MPI_DOUBLE,i,1000080,PETSC_COMM_WORLD);

         length=strlen(data->project_name);
         ierr=MPI_Send (&length,1,MPI_INT,i,1000081,PETSC_COMM_WORLD);
         ierr=MPI_Send(data->project_name,length,MPI_CHAR,i,1000082,PETSC_COMM_WORLD);

         i++;
      }

   } else {

      ierr=MPI_Recv(&(data->project_calculate_poynting),1,MPI_INT,0,1000001,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->project_save_fields),1,MPI_INT,0,1000002,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

      ierr=MPI_Recv(&length,1,MPI_INT,0,1000003,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      if (data->mesh_file) free(data->mesh_file);
      data->mesh_file=(char *) malloc((length+1)*sizeof(char));
      ierr=MPI_Recv(data->mesh_file,length,MPI_CHAR,0,1000004,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      data->mesh_file[length]='\0';

      ierr=MPI_Recv(&(data->mesh_order),1,MPI_INT,0,1000005,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->mesh_save_refined),1,MPI_INT,0,1000083,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->mesh_2D_refinement_fraction),1,MPI_DOUBLE,0,1000008,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->mesh_3D_refinement_fraction),1,MPI_DOUBLE,0,1000009,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->mesh_quality_limit),1,MPI_DOUBLE,0,1000010,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

      ierr=MPI_Recv(&length,1,MPI_INT,0,1000011,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      if (data->port_definition_file) free(data->port_definition_file);
      data->port_definition_file=(char *) malloc((length+1)*sizeof(char));
      ierr=MPI_Recv(data->port_definition_file,length,MPI_CHAR,0,1000012,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      data->port_definition_file[length]='\0';

      ierr=MPI_Recv(&length,1,MPI_INT,0,1000013,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      if (data->materials_global_path) free(data->materials_global_path);
      data->materials_global_path=(char *) malloc((length+1)*sizeof(char));
      ierr=MPI_Recv(data->materials_global_path,length,MPI_CHAR,0,1000014,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      data->materials_global_path[length]='\0';

      ierr=MPI_Recv(&length,1,MPI_INT,0,1000015,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      if (data->materials_global_name) free(data->materials_global_name);
      data->materials_global_name=(char *) malloc((length+1)*sizeof(char));
      ierr=MPI_Recv(data->materials_global_name,length,MPI_CHAR,0,1000016,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      data->materials_global_name[length]='\0';

      ierr=MPI_Recv(&length,1,MPI_INT,0,1000017,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      if (data->materials_local_path) free(data->materials_local_path);
      data->materials_local_path=(char *) malloc((length+1)*sizeof(char));
      ierr=MPI_Recv(data->materials_local_path,length,MPI_CHAR,0,1000018,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      data->materials_local_path[length]='\0';

      ierr=MPI_Recv(&length,1,MPI_INT,0,1000019,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      if (data->materials_local_name) free(data->materials_local_name);
      data->materials_local_name=(char *) malloc((length+1)*sizeof(char));
      ierr=MPI_Recv(data->materials_local_name,length,MPI_CHAR,0,1000020,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      data->materials_local_name[length]='\0';

      ierr=MPI_Recv(&length,1,MPI_INT,0,1000021,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      if (data->materials_default_boundary) free(data->materials_default_boundary);
      data->materials_default_boundary=(char *) malloc((length+1)*sizeof(char));
      ierr=MPI_Recv(data->materials_default_boundary,length,MPI_CHAR,0,1000022,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      data->materials_default_boundary[length]='\0';

      ierr=MPI_Recv(&(data->materials_check_limits),1,MPI_INT,0,1000023,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

      ierr=MPI_Recv(&length,1,MPI_INT,0,1000024,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      if (data->refinement_frequency) free(data->refinement_frequency);
      data->refinement_frequency=(char *) malloc((length+1)*sizeof(char));
      ierr=MPI_Recv(data->refinement_frequency,length,MPI_CHAR,0,1000025,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      data->refinement_frequency[length]='\0';

      ierr=MPI_Recv(&length,1,MPI_INT,0,1000026,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      if (data->refinement_variable) free(data->refinement_variable);
      data->refinement_variable=(char *) malloc((length+1)*sizeof(char));
      ierr=MPI_Recv(data->refinement_variable,length,MPI_CHAR,0,1000027,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      data->refinement_variable[length]='\0';

      ierr=MPI_Recv(&(data->refinement_iteration_min),1,MPI_INT,0,1000028,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->refinement_iteration_max),1,MPI_INT,0,1000029,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->refinement_required_passes),1,MPI_INT,0,1000030,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->refinement_relative_tolerance),1,MPI_DOUBLE,0,1000031,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->refinement_absolute_tolerance),1,MPI_DOUBLE,0,1000084,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

      ierr=MPI_Recv(&planCount,1,MPI_UNSIGNED_LONG,0,1000032,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

      j=0;
      while (j < planCount) {
         ierr=MPI_Recv(&planType,1,MPI_INT,0,1000033,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
         ierr=MPI_Recv(&planFrequency,1,MPI_DOUBLE,0,1000034,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
         ierr=MPI_Recv(&planStart,1,MPI_DOUBLE,0,1000035,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
         ierr=MPI_Recv(&planStop,1,MPI_DOUBLE,0,1000036,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
         ierr=MPI_Recv(&planStep,1,MPI_DOUBLE,0,1000037,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
         ierr=MPI_Recv(&planPointsPerDecade,1,MPI_INT,0,1000038,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
         ierr=MPI_Recv(&planRefine,1,MPI_INT,0,1000039,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
         ierr=MPI_Recv(&planLineNumber,1,MPI_INT,0,1000040,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

         add_inputFrequencyPlan (data,planType,planFrequency,planStart,planStop,planStep,planPointsPerDecade,planLineNumber,planRefine);
         
         j++;
      }

      ierr=MPI_Recv(&(data->reference_impedance),1,MPI_DOUBLE,0,1000041,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

      ierr=MPI_Recv(&length,1,MPI_INT,0,1000044,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      if (data->touchstone_frequency_unit) free(data->touchstone_frequency_unit);
      data->touchstone_frequency_unit=(char *) malloc((length+1)*sizeof(char));
      ierr=MPI_Recv(data->touchstone_frequency_unit,length,MPI_CHAR,0,1000045,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      data->touchstone_frequency_unit[length]='\0';

      ierr=MPI_Recv(&length,1,MPI_INT,0,1000046,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      if (data->touchstone_format) free(data->touchstone_format);
      data->touchstone_format=(char *) malloc((length+1)*sizeof(char));
      ierr=MPI_Recv(data->touchstone_format,length,MPI_CHAR,0,1000047,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      data->touchstone_format[length]='\0';

      ierr=MPI_Recv(&(data->solution_modes),1,MPI_INT,0,1000048,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->solution_temperature),1,MPI_DOUBLE,0,1000049,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->solution_2D_tolerance),1,MPI_DOUBLE,0,1000050,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->solution_3D_tolerance),1,MPI_DOUBLE,0,1000051,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->solution_iteration_limit),1,MPI_INT,0,1000052,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->solution_modes_buffer),1,MPI_INT,0,1000053,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

      ierr=MPI_Recv(&length,1,MPI_INT,0,1000054,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      if (data->solution_impedance_definition) free(data->solution_impedance_definition);
      data->solution_impedance_definition=(char *) malloc((length+1)*sizeof(char));
      ierr=MPI_Recv(data->solution_impedance_definition,length,MPI_CHAR,0,1000055,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      data->solution_impedance_definition[length]='\0';

      ierr=MPI_Recv(&length,1,MPI_INT,0,1000056,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      if (data->solution_impedance_calculation) free(data->solution_impedance_calculation);
      data->solution_impedance_calculation=(char *) malloc((length+1)*sizeof(char));
      ierr=MPI_Recv(data->solution_impedance_calculation,length,MPI_CHAR,0,1000057,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      data->solution_impedance_calculation[length]='\0';

      ierr=MPI_Recv(&(data->solution_check_closed_loop),1,MPI_INT,0,1000058,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->solution_check_homogeneous),1,MPI_INT,0,1000087,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->solution_accurate_residual),1,MPI_INT,0,1000059,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->solution_shift_invert),1,MPI_INT,0,1000060,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->solution_use_initial_guess),1,MPI_INT,0,1000061,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->solution_shift_factor),1,MPI_DOUBLE,0,1000062,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

      ierr=MPI_Recv(&(data->output_show_refining_mesh),1,MPI_INT,0,1000063,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->output_show_postprocessing),1,MPI_INT,0,1000064,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->output_show_iterations),1,MPI_INT,0,1000065,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->output_show_license),1,MPI_INT,0,1000066,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

      ierr=MPI_Recv(&(data->test_create_cases),1,MPI_INT,0,1000067,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->test_show_detailed_cases),1,MPI_INT,0,1000068,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

      ierr=MPI_Recv(&(data->debug_show_memory),1,MPI_INT,0,1000069,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->debug_show_project),1,MPI_INT,0,1000070,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->debug_show_frequency_plan),1,MPI_INT,0,1000071,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->debug_show_materials),1,MPI_INT,0,1000072,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->debug_show_port_definitions),1,MPI_INT,0,1000073,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->debug_show_impedance_details),1,MPI_INT,0,1000074,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->debug_save_port_fields),1,MPI_INT,0,1000075,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->debug_skip_mixed_conversion),1,MPI_INT,0,1000085,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->debug_skip_forced_reciprocity),1,MPI_INT,0,1000086,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      ierr=MPI_Recv(&(data->debug_tempfiles_keep),1,MPI_INT,0,1000076,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

      ierr=MPI_Recv(&(data->field_points_count),1,MPI_INT,0,1000077,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      data->field_points_allocated=data->field_points_count;

      if (data->field_points_x) free(data->field_points_x);
      data->field_points_x=(double *) malloc(data->field_points_count*sizeof(double));
      ierr=MPI_Recv(data->field_points_x,data->field_points_count,MPI_DOUBLE,0,1000078,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

      if (data->field_points_y) free(data->field_points_y);
      data->field_points_y=(double *) malloc(data->field_points_count*sizeof(double));
      ierr=MPI_Recv(data->field_points_y,data->field_points_count,MPI_DOUBLE,0,1000079,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

      if (data->field_points_z) free(data->field_points_z);
      data->field_points_z=(double *) malloc(data->field_points_count*sizeof(double));
      ierr=MPI_Recv(data->field_points_z,data->field_points_count,MPI_DOUBLE,0,1000080,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);

      ierr=MPI_Recv(&length,1,MPI_INT,0,1000081,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      if (data->project_name) free(data->project_name);
      data->project_name=(char *) malloc((length+1)*sizeof(char));
      ierr=MPI_Recv(data->project_name,length,MPI_CHAR,0,1000082,PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
      data->project_name[length]='\0';
 
   }

   return ierr;
}

