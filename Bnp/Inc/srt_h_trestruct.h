#ifndef SRT_H_TRESTRUCT_H         
#define SRT_H_TRESTRUCT_H

/* Common constants for tree discretisation */
#define TWOSQRTTWO					2.828427124746
#define DELTAMINBOUND				1.154700538379	/** 2/sqrt(3) **/
#define TWOFACDELTAMINBOUND			1.43			/** +- sqrt(2) **/
#define U_FACTOR					2.0				/*1.73*/ /*1.25*/ /*1.1547006*/
#define DELTAMAXBOUND				2.0
#define DELTAOPTRATIO				1.732050808
#define TREEMAXSTEP					500
#define	CUTTREEATSTDEV				5
#define NUMSTDEVINSPACING			1.414213562373 /* sqrt(2) */

/*	Maximum number of nodes in Nono tree for a given time step:
	return error if exceeded */
#define MAXNODE						400000			

/* Centering convention --- node closest in x or closest in log x */
#define CLOSESTINX 0
#define CLOSESTINLNX 1

/*-----------------------------------------------------------------------------
	TYPE            :SrtTwoFacMmtStruct 
	AUTHOR          :K L Chau (Modified AS)
	DESCRIPTION     :Store the first and second moments at each node for a two
					factor tree
-----------------------------------------------------------------------------*/

typedef struct 
{
	double E_X1;
	double E_X2;
	double E_X1X1;
	double E_X2X2;
	double E_X1X2;
} SrtTwoFacMmtStruct;

/*-----------------------------------------------------------------------------
  TYPE            :SrtTrinTreNdInf 
  AUTHOR          :O. Van Eyseren
  DESCRIPTION     :Information about a particular node in a trinomial tree
				  (whether it is LGM, CHeyette, Log...
-----------------------------------------------------------------------------*/

typedef struct {
	SrtSample 	cur_sam;
	SrtSample 	drift_sam;
   	double 		p[3];
	double 		df;				/* df to next time point */  
	long		mid_son_index;
    double 		var_at_sam;
 
	long		son_index[3];
	double		son_state_var_value[3];
   	double		state_var_variance;
	double		state_var_expectation;
	SrtMdlType	mdl_type;
     } SrtTrinTreNdInf;

/*-----------------------------------------------------------------------------

  TYPE            :SrtNonoTreeNodeInfo 
  AUTHOR          :K L Chau (Modified AS)
  DESCRIPTION     :Information about a particular node in a nonomial tree
-----------------------------------------------------------------------------*/

typedef struct 
{
	SrtSample 			cur_sam;			/* X1 and X2 (original basis) */
   	double 				p[3][3];			/* Probas */
	double 				df;					/* Df to next time point */
	double				forward[2];			/* EX1 and EX2 (dual basis) */
	long				  son_index[2][3];	/* Son Indexes (dual basis) */
	double				son_statevar[2][3];	/* Son State Vars (dual basis) */

	/* These are used in srt_f_twoundtre.c */
	/* They are useless for lgm2d tree */
	SrtTwoFacMmtStruct	moments;				
	SrtSample 			drift_sam;
	double				var_at_sam[2];
	double				correl;
	double				delta_t;
	long				mid_son_index[2];

}  SrtNonoTreeNodeInfo;

/*-----------------------------------------------------------------------------

  TYPE              : SrtPentoTreeNodeInfo 
  AUTHOR         : EF
  DESCRIPTION : Information about a particular tree node 
-----------------------------------------------------------------------------*/

typedef struct 
{
	SrtSample 		 cur_sam;			   /* X1 and X2 (original basis) */
   	double 				p[5];			          /* Probas */
	double 				df;					       /* Df to next time point */
	double				forward[2];			    /* EX1 and EX2 (dual basis) */
	long				  son_index[2][3];	  /* Son Indexes (dual basis) */
	double				son_statevar[2][3];	/* Son State Vars (dual basis) */
	
	double              quanto_fwd;
}  
SrtPentoTreeNodeInfo;

/*-----------------------------------------------------------------------------
  TYPE            :SrtTreStpInf 
  AUTHOR          :O. Van Eyseren
  DESCRIPTION     :Information about a particular time step in a tree
-----------------------------------------------------------------------------*/

typedef struct{
	long	max_x_index;
	long	min_x_index;
	double	logxmin;	/* Used in lognormal models */
	double	xmin;		/* Used in normal models */
	double	u;			/* Multiplicative or additive*/
} SrtTreStpInf;

/*-----------------------------------------------------------------------------
  TYPE            :SrtCheTreInf 
  AUTHOR          :E.Auld
  DESCRIPTION     :Information about a particular time in a cheyette tree
-----------------------------------------------------------------------------*/

typedef struct{
	long	max_r_index;
	long	min_r_index;
	long	max_phi_index;
	long	min_phi_index;
	double	logrmin;		/* Used in lognormal models */
	double	rmin;			/* Used in normal models */
	double	u;				/* Multiplicative or additive*/
} SrtCheTreInf;

/*-----------------------------------------------------------------
  TYPE            :SrtTwoFacTreInf 
  AUTHOR          :K L Chau
  DESCRIPTION     :information about a particular time in a nonomial tree
-------------------------------------------------------------------*/

typedef struct
{
	long	max_x_index[2];
	long	min_x_index[2];
	double	logxmin[2];		/* Used in lognormal models */
	double	xmin[2];		/* Used in normal models */
	double  xmax[2];
	double	u[2];			/* Multiplicative or additive*/
} SrtTwoFacTreInf, SrtNonoTreStpInf;


/*-----------------------------------------------------------------
  TYPE               :SrtDiagTwoFacTreeInfo 
  AUTHOR          :Antoine Savine
  DESCRIPTION  :Tree Geometry for "diagonalised" Nono used for LGM2F 
-------------------------------------------------------------------*/

typedef struct
{
	/* Dual basis = eigenvectors of the local cov matrix * sqrt (eigenvalues) */
	/* Coordinates of dual basis in terms of original basis */
	double	dual_basis[2][2];
	double	inverse_basis[2][2];
	/* Spacing and indexation are described in the dual basis */
	double	spacing[2];
	long	max_index[2];

} SrtDiagTwoFacTreeInfo;

#endif
