/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : newton.h                                                     */
/*                                                                            */
/* DESCRIPTION : Header for use of various numerical analysis routines        */
/*                                                                            */
/* DATE        : Wed Apr 17 1996                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/
 
 

#ifndef _NEWTON_H
#define _NEWTON_H



class ARM_Vector;
class ARM_Matrix;


typedef int (*PFUNC)(ARM_Vector *, ARM_Matrix *, ARM_Vector *, void **);

extern	double newtonRoot(
	double (*)(double &, double &, double, void **), 
	double, 
	double, 
	void **, 
	double, 
	int);

extern int	multiDimNewtonRoot(
    PFUNC		function,
	ARM_Vector 	*root,
	void		**fixedParams, 
	double		tol, 
	int			itermax);





#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
