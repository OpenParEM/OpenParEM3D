#ifndef LINEARSYSTEM_H
#define LINEARSYSTEM_H

#include <petsc.h>
#include <complex.h>
#include "project.h"
#include <_hypre_parcsr_mv.h>
#include "triplet.h"

//double* allocReaddof (char *, char *, size_t *);
//void printdof (double *, size_t);
FILE* openDataFile (const char *, const char *, char *, int);
int loadDataLine (FILE *, struct dataTriplet *, int);
int loadDataFileStats (const char *, const char *, char *, PetscInt *, PetscInt *, PetscInt *);
int loadDataFile (const char *, const char *, char *, Mat *, PetscInt, PetscInt, int, int, double, PetscMPIInt);

#endif

