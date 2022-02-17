/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE         : randgen.cpp                                                 */
/*                                                                            */
/* DESCRIPTION  : Random generators for Monte carlo Kernels                   */
/*                                                                            */
/*                Sep. 10th 1997                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/
#ifndef RANDGEN_H
#define RANDGEN_H 
 

extern double myRipleyRandom(void);

extern void init_Rand(int deterministe, int meth);

extern void HammersleyAndPremiers(int dimInteg, long nbTir);

extern double DinfiStar(int dimInteg, long nbTir);

extern double Dinfi(int dimInteg, long nbTir);

extern double Jinfi(int dimInteg, long nbTir);

extern long IntValue(double v);

extern int* G_PREMIER;
extern int PREM_PREV_SIZE;
extern double** G_HAMMERSLEY;
extern long G_NB_TIRAGES;
extern int G_NB_DIM;


extern void FreeHammersley(void);



#endif
/*----------------------------------------------------------------------------*/
