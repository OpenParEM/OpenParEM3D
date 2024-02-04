1. check on project files

   (a) find . -name ".*.lock" -print

       Remove any locks found.
       Verify that a valid job is not running before removing a lock file.

   (b) proj_search has output.show.license | grep true
       proj_search has debug.show.mode.definitions | grep true
       proj_search has debug.show.impedance.details | grep true
       proj_search has debug.tempfiles.keep | grep true
       proj_search has project.save.fields | grep true
       proj_search has test.create.cases | grep true
       proj_search has debug.save.port.fields | grep true
       proj_search has debug.tempfiles.keep | grep true

       Generally, these should be false to minimize disk usage.

2. To set up a project:

      mpirun -np N OpenParEM2D *.proj
      mv *_prototype_test_cases.csv *_test_cases.csv
      Edit the file *_test_cases.csv as desired.
      Set test.create.cases to false

   To check a setup project for errors:
   process3D.sh project.proj N

   Adjust *.proj settings and *_derived_test_cases.csv error limits as desired to complete the setup.
   Setup is complete when all tests pass.

3. To check a project for errors:

   process3D.sh project.proj N

4. To check errors in all projects:

   mv regression.log regression.log.archive
   regression3D.sh >& regression.log

   from the regression directory.  

   Check "regression.log" for "ERROR", "NOT CONVERGED", and "terminated".
   Check "regression_results.csv" for "FAIL".  Waive small errors.

5. To add projects to the regression suite, add a line to

   regression_case_list.txt



