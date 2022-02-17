/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file barriere_bs.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_BARRIERE_BS_H
#define _GP_CF_BARRIERE_BS_H

#include "firsttoinc.h"
#include "gpbase/port.h"


CC_BEGIN_NAMESPACE(ARM)



//////////////////////////////////////////////////////////////////////////////////////////////
///  
///			Begining Exportable  Pricing Functions 
///
//////////////////////////////////////////////////////////////////////////////////////////////



double Export_BS_EuroBarriere(double f,double k, double b, double r, double v, double t,double rate,
								   int callput, int inout, int updown);

double Export_BS_EuroBarriere_ImpliedVol(double f,double k, double b, double r, double opt, double t,double rate,
										 int callput, int inout, int updown);

double Export_BS_EuroDoubleBarriere(double f,double k, double bup, double bdown, double v, double t,double r,double b,
								   int callput);

double Export_BS_EuroDoubleBarriere_ImpliedVol(double f,double k, double bup, double bdown, double opt, double t,double r,double b,
								   int callput);

double Export_BS_PartialTime_Start_SingleBarrier(double f,double k, double barrier, double rebate, double v,
												 double bendtime, double t,int callput, int optype);

double Export_BS_PartialTime_End_SingleBarrier(double f,double k, double barrier, double rebate, double v,
												 double bstarttime, double t,int callput, int optype);

double Export_BS_SingleBarrier_2Asset(double f1,double k1, double f2,double k2, double v1, double v2,double corr,
												double t,int callput, int optype);

///////////////////////////////////////////////////////////////////////
///  
///			1st Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_BS_EuroBarriere(int i,double f,double k, double b, double r, double v, double t,double rate,
								   int callput, int inout, int updown);
double Export_BS_EuroDoubleBarriere(int i,double f,double k, double bup, double bdown, double v, double t,
								   int callput);
double Export_BS_PartialTime_Start_SingleBarrier(int i,double f,double k, double barrier, double rebate, double v,
												 double bendtime, double t,int callput, int optype);
double Export_BS_PartialTime_End_SingleBarrier(int i,double f,double k, double barrier, double rebate, double v,
												 double bstarttime, double t,int callput, int optype);

double Export_BS_SingleBarrier_2Asset(int i,double f1,double k1, double f2,double k2, double v1, double v2,double corr,
												double t,int callput, int optype);


///////////////////////////////////////////////////////////////////////
///  
///			2nd Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_BS_EuroBarriere(int i,int j,double f,double k, double b, double r, double v, double t,double rate,
								   int callput, int inout, int updown);
double Export_BS_EuroDoubleBarriere(int i, int j,double f,double k, double bup, double bdown, double v, double t,
								   int callput);
double Export_BS_PartialTime_Start_SingleBarrier(int i, int j,double f,double k, double barrier, double rebate, double v,
												 double bendtime, double t,int callput, int optype);
double Export_BS_PartialTime_End_SingleBarrier(int i, int j,double f,double k, double barrier, double rebate, double v,
												 double bstarttime, double t,int callput, int optype);

double Export_BS_SingleBarrier_2Asset(int i, int j,double f1,double k1, double f2,double k2, double v1, double v2,double corr,
												double t,int callput, int optype);



///////////////////////////////////////////////////////////////////////
///  
///			End Exportable Pricing Functions
///
///////////////////////////////////////////////////////////////////////




CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

