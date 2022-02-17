/*
 * $Log: interpol.h,v $
 * Revision 1.21  2003/08/26 11:30:53  ebenhamou
 * remove commented code
 *
 * Revision 1.20  2003/07/31 16:34:56  mab
 * in : inline int indexBeforeValue(double* xxElt, int size, double x)
 * instead of throwing an exception return -1
 *
 * Revision 1.19  2003/07/31 16:07:58  jmprie
 * modif definition de InterpolateZc()
 *
 * Revision 1.18  2003/07/25 15:41:57  emezzine
 * Correct SearchIdx()
 *
 * Revision 1.17  2003/07/15 14:33:28  jmprie
 * ajout d'un flag pour traiter correctectement de l'interpolation au
 * tout debut de la courbe ds InterpolateZc()
 *
 * Revision 1.16  2003/07/11 12:17:56  emezzine
 * Add SearchIndex()
 *
 * Revision 1.15  2003/06/18 14:27:57  mab
 * Formatting
 *
 * Revision 1.14  2003/06/18 13:02:13  jmprie
 * ajout de  PolynomialInterpolation
 *
 * Revision 1.13  2003/05/30 08:34:35  jmprie
 * creation des InterpolateForwardRate() et InterpolateZc qui implementent
 *  des interpolations particulieres a partir une liste de tx forwards pour
 * calculer soit d'autres tx forward soit des Zc
 *
 * Revision 1.12  2003/04/29 13:03:33  mab
 * Correction in : inline int locateIndex(double* xxElt, int size, double x)
 *
 * Revision 1.11  2003/04/23 13:09:37  mab
 * in : inline int locateIndex(double* xxElt, int size, double x)
 * Treatement case : ( size == 1 )
 *
 * Revision 1.10  2003/03/24 10:35:24  mab
 * improvements in locateIndex
 * ADDED : boolTest
 *
 * Revision 1.9  2003/03/21 17:54:41  mab
 * int locateIndex
 * static j replaced static j = 0
 *
 * Revision 1.8  2003/02/10 17:56:48  mab
 * Added : double triangularInterpol(..)
 *
 * Revision 1.7  2001/05/23 17:46:19  smysona
 * correction locateIndex
 *
 * Revision 1.6  2001/04/27 09:37:09  smysona
 *  Palier pour les interpol de vol
 *
 * Revision 1.5  2001/04/23 09:28:36  smysona
 *  transfer de indexBeforeValue vers le header (inline)
 *
 * Revision 1.4  2001/04/03 12:05:14  nicolasm
 * Proto avec double*
 *
 * Revision 1.3  1999/02/15 09:45:47  ypilchen
 * Rajout de : extern
 * :
 *
 * Revision 1.2  1998/11/27 12:04:44  nicolasm
 * ajout de la fonction linInterpol qui prend 5 double en argument
 *
 */

/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : interpol.h                                                   */
/*                                                                            */
/* DESCRIPTION : Interpol. utilities                                          */
/*                                                                            */
/* DATE        : Wed Apr 17 1996                                              */
/*                                                                            */
/*----------------------------------------------------------------------------*/
#ifndef _INTERPOL_H
#define _INTERPOL_H



#include "linalg.h"




inline  int     locateIndex(const double* xxElt, int size, double x);
inline  int     locateIndex(const double* xxElt, int size, double x,int prevlocation);
inline  int     indexBeforeValue(const double *, int, double);

inline int      SearchIdx(ARM_Vector* xxElt, double x, double );

inline int      ExactSearchValue(ARM_Vector* xList,double x,int prevIdx,
                                 double xTol,char* errMsge);

extern  int     indexAfterValue(ARM_Vector *, double);

extern  double  linInterpol(ARM_Vector *, double, ARM_Vector *);
extern double linInterpolCol(const ARM_Vector &, double, const ARM_Matrix&,int col);
extern  double  linInterpol(double *, int, double, double *);
extern  double SlopelinInterpol(double* xElt, int xSize, double xval, double* yElt);
extern  double  linInterpol2(ARM_Vector *, double, ARM_Vector *, 
                             bool step = true);
extern  double  linInterpol2Col(const ARM_Vector &, double, const ARM_Matrix&, int col,
                             bool step = true);
extern  double  D1LinInterpol(ARM_Vector *, double, ARM_Vector *);
extern  double  D2LinInterpol(ARM_Vector *, double, ARM_Vector *);
extern  double  twoDimLinInterpol(ARM_Vector *, ARM_Vector *, 
                                  double, double, ARM_Vector *);
extern double linInterpol(double x, double x1, double y1, double x2,
                          double y2);

extern double triangularInterpol(ARM_Vector* x, ARM_Vector* y,
                                 ARM_Matrix* z, double xval, 
                                 double yval);

extern double SplineInterpolateFunc(ARM_Vector* vect1, ARM_Vector* vect2,
                                    double inX,
                                    ARM_Vector* SecondDerivCalc = NULL,
                                    int keep2Der = 0);

/*  Fonctions d'interpolation ds une liste de taux forward visant a
    calculer soit des forwards de meme tenor soit des zero-coupons  */
     
extern int FindNextFwdIdx(double calcDate,double start,
                          ARM_Vector* ycFwdResetDate,
                          ARM_Vector* ycFwdStartDate,
                          double terminalDate);

extern ARM_Matrix* InterpolateForwardRate(double calcDate,
                                          ARM_Vector* ycFwdResetDate,
                                          ARM_Vector* ycFwdStartDate,
                                          double ycTerminalDate,
                                          double ycRateConvert,
                                          ARM_Matrix* forward,
                                          ARM_Vector* rateStartDate);

extern ARM_Matrix* InterpolateZc(   double calcDate,
                                    ARM_Vector* ycFwdResetDate,
                                    ARM_Vector* ycFwdStartDate,
                                    double ycTerminalDate,
                                    int ycFwdDayCount,
                                    ARM_Matrix* forward,
                                    ARM_Matrix* savAcc,
                                    ARM_Vector* zcDate,
                                    void* miscData[]);

/*---- Calcul des coefs du polynome passant par m+1 points ----*/
extern ARM_Vector* PolynomialInterpolation(ARM_Vector* x,
                                           ARM_Vector* fx);

// Interpolation based on Lagrange Polynom
extern double QuadraticLagrangeInterpol(double x,
                                        double x0, double y0,
                                        double x1, double y1,
                                        double x2, double y2);

// Interpolation based on Lagrange Polynom N=4
extern double LagrangeInterpol(double x,
                               double x0, double y0,
                               double x1, double y1,
                               double x2, double y2,
                               double x3, double y3);

inline int indexBeforeValue(ARM_Vector* xx, double x)
{
    return indexBeforeValue(xx->GetElt(), xx->GetSize(), x);
}


inline int locateIndex(ARM_Vector* xx, double x)
{
    return locateIndex(xx->GetElt(), xx->GetSize(), x);
}

inline int locateIndex(ARM_Vector* xx, double x,int prevlocation)
{
    return locateIndex(xx->GetElt(), xx->GetSize(), x,prevlocation);
}



/*----------------------------------------------------------------------------*
  Assuming that the coeffs of xx are sorted in increasing order, this routine
    returns the largest index j such that xx->Elt(j) <= x.
    
	handles the case of x <= than xx->Elt(0) and x>= xx->Elt(size-1)
*----------------------------------------------------------------------------*/

inline int indexBeforeValue(const double* xxElt, int size, double x)
{
    int k;


    if ( x < xxElt[0] - K_DOUBLE_TOL )
    {
       return(-1);
    }

    if ( x >= xxElt[size-1] ) 
       return(size-1);
    
    k = locateIndex(xxElt, size, x);
    
    return(k);
}



inline int locateIndex(const double* xxElt, int size, double x)
{
    int ascnd, ju, jm, jl;

    static int j = 0;


    if ( size == 1 )
    {
       if ( x < xxElt[0]-K_DOUBLE_TOL )
       {
          return(int(-1));
       } 
       else
       {
          return(int(0));
       }
    }

    jl = -1;
    ju = size;

    ascnd = xxElt[size-1] > xxElt[0];

    if ( (j < size) && (j>=0) && (( x >= xxElt[j] - K_DOUBLE_TOL ) == ascnd ))
       jl = j;


    while ( ju-jl > 1) 
    {
        jm = (ju+jl)/2;

        int boolTest;
       
        boolTest = (( x >= xxElt[jm]-K_DOUBLE_TOL ) == ascnd ) 
                   ||
                   (( x == xxElt[jm] ) == ascnd );

        if (boolTest)
        {
           jl = jm;
        }
        else
        {
           ju = jm;
        }
    }

    j = jl;

    return(j);
}


inline int locateIndex(const double* xxElt, int size, double x,int prevlocation)
{
    int ascnd, ju, jm, jl;

    // static int j = 0;
	int j=prevlocation; 

    if ( size == 1 )
    {
       if ( x < xxElt[0]-K_DOUBLE_TOL )
       {
          return(int(-1));
       } 
       else
       {
          return(int(0));
       }
    }

    jl = -1;
    ju = size;

    ascnd = xxElt[size-1] > xxElt[0];

    if ( (j < size) && (j>=0) && (( x >= xxElt[j] - K_DOUBLE_TOL ) == ascnd ))
       jl = j;


    while ( ju-jl > 1) 
    {
        jm = (ju+jl)/2;

        int boolTest;
       
        boolTest = (( x >= xxElt[jm]-K_DOUBLE_TOL ) == ascnd ) 
                   ||
                   (( x == xxElt[jm] ) == ascnd );

        if (boolTest)
        {
           jl = jm;
        }
        else
        {
           ju = jm;
        }
    }

    j = jl;

    return(j);
}


int SearchIdx(ARM_Vector* xxElt,double x, double tol_)
{
    int ascnd, ju, jm, jl;
    int size = xxElt->GetSize();

    static int jj = 0;
    int j = jj;

    if ( size == 1 )
    {
       if ( x < xxElt->Elt(0) -K_DOUBLE_TOL )
       {
          return(int(-1));
       } 
       else
       {
          return(int(0));
       }
    }
  
    jl = -1;
    ju = size;

    ascnd = xxElt->Elt(size-1) > xxElt->Elt(0);

    if ( (j < size) && (j>=0) && (( x >= xxElt->Elt(j) - tol_ ) == ascnd ))
       jl = j;


    while ( ju-jl > 1) 
    {
        jm = (ju+jl)/2;

        if ((( x >= xxElt->Elt(jm) - tol_) == ascnd) || (x == xxElt->Elt(jm)) )
        {
           jl = jm;
        }
        else
        {
           ju = jm;
        }
    }

    j = jl;
    jj=j;

    return(j);
}


inline int ExactSearchValue(ARM_Vector* xList,double x,int prevIdx,
                            double xTol,char* errMsge)
{
    int idx,nbX=(xList ? xList->GetSize() : 0);
    double curX;

	for(idx=prevIdx;idx<nbX;++idx)
    {
        curX=xList->Elt(idx);
		if(curX - xTol <= x && x <= curX + xTol)
			break;
    }
    if(idx==nbX)
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,errMsge);
    }
 
    return idx;
}




#endif
/*----------------------------------------------------------------------------*/
/*---- Endif ----*/
