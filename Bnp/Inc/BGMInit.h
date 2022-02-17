#ifndef BGMINIT_H
#define BGMINIT_H

/*	Specific parameters of model	*/

Err SrtFreeBGMUnd(void *undPtr);

/* Function to free an underlying object in a linked list  */
Err srt_f_BGMdescfree(SrtUndDesc *spec_desc);

Err SrtInitBGMUnd(char *undName, char *ycName, char *cRefRateCode,
                  double *shifts, int num_of_factors, int num_of_tenors,
                  double **VolMatrix, double **histcorrel);

Err SrtInitBGMMidatUnd(char *ycName, char *cRefRateCode, double *shifts,
                       int num_of_factors, int num_of_tenors,
                       double **VolMatrix, double **histcorrel, SrtUndPtr *und);

Err srt_f_BGMSABRdescfree(SrtUndDesc *spec_desc);

Err SrtInitBGMSABRUnd(char *undName, char *ycName, char *cRefRateCode,
                      int MaxNumPeriod, int MaturityInPeriod, double **TSLibor,
                      double **Correl, long *IStartLiquidATM,
                      long *IEndLiquidATM, int NumLiquidATM, int NumMatCube,
                      int NumUndCube, long *IMatCube, long *IUndCube,
                      double ***ATMSens, long *IStartLiquidSABR,
                      long *IEndLiquidSABR, int NumLiquidSABR,
                      double ***AlphaSens, double ***RhoSens, double bumpATM,
                      double bumpalpha, double bumprho);

#endif /* END BGMInit.h */