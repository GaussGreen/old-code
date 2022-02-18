/* ========================================================================
   FILENAME:     pdestruct.h
   PURPOSE:     all the structures required to play with PDE's
   ======================================================================== */
#ifndef PDE_H_STRUCT
#define PDE_H_STRUCT

#include "pde_h_types.h"
#include "utallhdr.h"

/* A Generic type for the functions that compute the coefficients of the Differential Linear
 * Operator*/
typedef double (*PdeCoeffModelFunc)(void*, double[3], double);

/* A Generic type for the function that compute all sorts of thing for
   the change of variables and the localization of the Differential Linear Operator	 */
typedef Err (*PdeSetUpModelFunc)(void*);

/*
/* A Generic type for the function that compute all sorts of thing for
   the change of variables and the localization of the Differential Linear Operator
typedef  Err (*PdeSetUpModelFunc)(struct SrtPdeObject*);

/* A Generic type for the updating of info pointer of the Differential Linear Operator
typedef   Err (*PdeAttachInfoModelFunc)(struct SrtPdeObject*, void *);
*/

/* --------------------------------------------------------------------------------------------------------------------------------------------
        TYPE                :   SrtPdeGridObject
        DESCRIPTION   :   A structure that represents the PDE grid and the operator at a given time
step
/* A pointer to the discretised PDE  (containing the grid, the operator, ...)
   --------------------------------------------------------------------------------------------------------------------------------------------*/
typedef struct
{
    int dimension;
    int nb_payoffs;
    int nb_coeffs;

    /* Number of space points in each dimension */
    /* Indices go from 0 to num_X + 1 */
    long num_x;
    long num_y;
    long num_z;

    /* The Value functions : un[coeff][Z][Y][X][payoff] */
    double***** new_value;
    double***** old_value;
    double***** temp_value;

    /* The discrete values of  X, Y, Z  */
    double* vec_x;
    double* vec_y;
    double* vec_z;

    /* The increments of the values from X(i) to X(i+1) */
    double* inc_x;
    double* inc_y;
    double* inc_z;

    /* The operator for ADI (tridiagonal matrixes) */
    double*** AX;
    double*** BX;
    double*** CX;
    double*** AX1;
    double*** BX1;
    double*** CX1;

    double*** AY;
    double*** BY;
    double*** CY;
    double*** AY1;
    double*** BY1;
    double*** CY1;

    double*** AZ;
    double*** BZ;
    double*** CZ;
    double*** AZ1;
    double*** BZ1;
    double*** CZ1;

    double*** EXY;
    double*** EXZ;
    double*** EYZ;

    double*** L1;
    double*** D1;
    double*** L2;
    double*** D2;
    double*** L3;
    double*** D3;

} SrtPdeGridObject;

/* --------------------------------------------------------------------------------------------------------------------------------------------
        TYPE                :   SrtPdeObject
        DESCRIPTION   :   A structure to discretise and solve the PDE at a given time step
   --------------------------------------------------------------------------------------------------------------------------------------------*/
typedef struct
{
    /* The Dimension of the problem: one, two or three space variables */
    int dimension;

    /* The 	grid */
    SrtPdeGridObject* grid;

    /* For the Mesh style:	regular, exponential,...*/
    SrtPdeMeshType mesh_type_x;
    SrtPdeMeshType mesh_type_y;
    SrtPdeMeshType mesh_type_z;

    /* For the Boundary conditions:	Neuman, Dirichlet,... */
    SrtPdeBoundaryCond bound_cond_x;
    SrtPdeBoundaryCond bound_cond_y;
    SrtPdeBoundaryCond bound_cond_z;

    /* The Solving Scheme */
    SrtPdeSolvingScheme scheme;

    double theta;

    /* When computing differential of the solution wrt parameters,  a source member comes along:
     * same L, new PDE */
    int nb_second_members;

    /* The Number of payoffs to work on (also num of columns in the Grfn tableau,...)  */
    int nb_payoffs;

    /* For all the time information, the coefficients for the linear operator,...*/
    void* info;

    /* Functions to compute the coeffficients of the linear differential operator at a givent point
     * (x,y,z) */
    PdeCoeffModelFunc coeff_0;
    PdeCoeffModelFunc coeff_X;
    PdeCoeffModelFunc coeff_XX;
    PdeCoeffModelFunc coeff_Y;
    PdeCoeffModelFunc coeff_YY;
    PdeCoeffModelFunc coeff_Z;
    PdeCoeffModelFunc coeff_ZZ;
    PdeCoeffModelFunc coeff_XY;
    PdeCoeffModelFunc coeff_XZ;
    PdeCoeffModelFunc coeff_YZ;

    /* Localization function the extrema values
            PdeSetUpModelFunc             pde_set_up;  */

    /* The Extrema values */
    double min_x;
    double min_y;
    double min_z;
    double max_x;
    double max_y;
    double max_z;

    /*	PdeAttachInfoModelFunc		pde_attach_info_model; */

} SrtPdeObject, *SrtPdeObjectPtr;

typedef struct
{
    long num_mesh;
    long num_step;
    long num_pde;

    double h;
    double time_step;
    double log_spot_min;
    double spot;
    double spot_min;
    double spot_max;
    long   spot_index;

    double* log_spots;
    double* spots;

    double* time;

    double* basevol;
    double* beta;
    double* gamma;
    double* omega;

    double* voldrift;
    double* vovol;

    double* correlation_und_vol;
    double* correlation_und_rate;

    double* df;
    double* forward;

} SrtBasicPdeInf;

#endif
/*-------------------------------------------------------------------*/
/*--- End of File ---*/
