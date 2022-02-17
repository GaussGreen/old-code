/*
 * $Log: interpol.cpp,v $
 * Revision 1.19  2003/09/24 11:18:44  jmprie
 * traitement en corrigeant le drift ds le cas ou une date de Zc est entre la
 * sliceDate et la startDate du forward "spot" rajoute
 *
 * Revision 1.18  2003/08/26 11:30:33  ebenhamou
 * remove commented code
 *
 * Revision 1.17  2003/07/31 16:08:55  jmprie
 * ajout d'un tableau de void * pour passer des parametres supplementaires
 * pour une interpolation plus juste avec un forward spot et une correction
 * en drift
 *
 * Revision 1.16  2003/07/17 14:41:17  jmprie
 * ds InterpolateZc() modif interpolation qd t est entre deux starts de fwd de
 * diffusion
 *
 * Revision 1.15  2003/07/15 14:36:32  jmprie
 *  ajout d'un flag pour traiter correctectement de l'interpolation au
 * tout debut de la courbe ds InterpolateZc()
 *
 * Revision 1.14  2003/06/30 07:42:13  jmprie
 *  correction fuite memoire ds PolynomialInterpolation
 *
 * Revision 1.13  2003/06/24 14:22:00  jmprie
 * correction bug ds calcul du polynome de degre m passant par m+1 pts
 *
 * Revision 1.12  2003/06/18 14:28:10  mab
 * Spline Interpol added
 *
 * Revision 1.11  2003/06/18 13:01:15  jmprie
 * correction bug ds FindFwdNextIdx() et ajout de PolynomialInterpolation
 *
 * Revision 1.10  2003/05/30 08:32:43  jmprie
 * creation des InterpolateForwardRate() et InterpolateZc qui implementent
 * des interpolations particulieres a partir une liste de tx forwards pour
 * calculer soit d'autres tx forward soit des Zc
 *
 * Revision 1.9  2003/02/10 17:57:37  mab
 * Added : double triangularInterpol(..)
 *
 */


/*----------------------------------------------------------------------------*
    interpol.cpp
 
    Interpolation routines
 
*----------------------------------------------------------------------------*/


#include <stdio.h>
#include <math.h>


#include "expt.h"
#include "armglob.h"
#include "dates.h"

#include "linalg.h"
#include "interpol.h"

//#include "zerocurv.h"
#include "currency.h"






/*----------------------------------------------------------------------------*
    Assuming the coeffs of xx are sorted in either increasing or 
    decreasing order, this routine returns the index j such that x is in
    [xx->Elt(j), xx->Elt(j+1)[.
*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*
    Assuming that the coeffs of xx are sorted in increasing order, this routine
    returns the smallest index j such that xx->Elt(j) >= x.
*----------------------------------------------------------------------------*/

int indexAfterValue(ARM_Vector* xx, double x)
{
    int    k;

    
    if ( x > xx->Elt(xx->GetSize()-1) + K_DOUBLE_TOL) 
    {
       throw Exception(__LINE__, __FILE__, ERR_CONDITION_NOT_MEET,
                       "x is greater than all vector elements");
    }

    if ( x <= xx->Elt(0) ) 
       return(0);
    
    k = locateIndex(xx, x);
    
    if ( x == xx->Elt(k) )
    {
       return(k);
    }
    else 
    {
       return(k+1);
    }
}



/*----------------------------------------------------------------------------*
    Returns value of linear interpolation of y = f(x) at point xval. Returns
    kHugeDouble if failed.
*----------------------------------------------------------------------------*/

double linInterpol(ARM_Vector* x, double xval, ARM_Vector* y)
{
    double* xElt = x->GetElt();
    double* yElt = y->GetElt();

    
    if (( x == NULL ) || ( y == NULL ))
    {
       throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
                "The two vectors should be not NULL");

       return(-99999999999);
    }

    int xSize = x->GetSize();
    int ySize = y->GetSize();


    if ( xSize != ySize )
    {
       throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
                          "The two vectors should have the same size");

       return(-99999999999);
    }

    return(linInterpol(xElt, xSize, xval, yElt));
}



double linInterpol(double* xElt, int xSize, double xval, double* yElt)
{
    double yval, lambda;
    int    index;

    
    // Get index of point right before xval
    // or the appropriate index if xval does not belong to (xmin, xmax)

	if(xSize==1) return (yElt)[0];

    if ( xval < xElt[0]-K_DOUBLE_TOL ) 
    {
       index = 0;
    }
    else if ( xval > xElt[xSize-1]+K_DOUBLE_TOL) 
    {
       index = xSize-2;
    }
    else 
    {
       index = indexBeforeValue(xElt, xSize, xval);
    }

    //  do linear interpolation

    if ( xval == (xElt)[index] ) 
    {
       return((yElt)[index]);
    }

    lambda = ((xElt)[index+1]-xval)/((xElt)[index+1]-(xElt)[index]);

    yval = lambda*(yElt)[index]+(1.0-lambda)*(yElt)[index+1];

    return(yval);
}


double linInterpolCol(const ARM_Vector& xElt, double xval, const ARM_Matrix &yy,int col)
{
    double yval, lambda;
    int    index;

    
    // Get index of point right before xval
    // or the appropriate index if xval does not belong to (xmin, xmax)

	if(xElt.size()==1) return yy.Elt(0,col);

    if ( xval < xElt[0]-K_DOUBLE_TOL ) 
    {
       index = 0;
    }
    else if ( xval > xElt[xElt.size()-1]+K_DOUBLE_TOL) 
    {
       index = xElt.size()-2;
    }
    else 
    {
       index = indexBeforeValue(xElt.GetElt(), xElt.size(), xval);
    }

    //  do linear interpolation

    if ( xval == xElt[index] ) 
    {
       // return((yElt)[index]);
		return yy.Elt(index,col); 
    }

    lambda = (xElt[index+1]-xval)/(xElt[index+1]-xElt[index]);

    // yval = lambda*(yElt)[index]+(1.0-lambda)*(yElt)[index+1];
	yval = lambda* yy.Elt(index,col) +(1.0-lambda)*yy.Elt(index+1,col);

    return(yval);
}


double SlopelinInterpol(double* xElt, int xSize, double xval, double* yElt)
{
    double lambda;
    int    index;

    if(xSize==1) return (yElt)[0];
    // Get index of point right before xval
    // or the appropriate index if xval does not belong to (xmin, xmax)


    if ( xval < xElt[0]-K_DOUBLE_TOL ) 
    {
       index = 0;
    }
    else if ( xval > xElt[xSize-1]+K_DOUBLE_TOL) 
    {
       index = xSize-2;
    }
    else 
    {
       index = indexBeforeValue(xElt, xSize, xval);
    }

    //  do linear interpolation

    lambda = ((yElt)[index+1]-(yElt)[index])/((xElt)[index+1]-(xElt)[index]);

    return(lambda);
}


/*----------------------------------------------------------------------------*
    Returns value of linear interpolation of y = f(x) at point xval. Returns
    kHugeDouble if failed.

*----------------------------------------------------------------------------*/
double  linInterpol2(ARM_Vector* x, double xval, ARM_Vector* y, bool step)
{  
    if ( x->GetSize() != y->GetSize() )
    {
       throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
                       "The two vectors should have the same size");
    }
 
    // case Size = 1
    long size = x->GetSize();

    if ( size == 1 )
    {
       return(y->Elt(0));
    }
 
    if (step)
    {
       if ( xval >= x->Elt(size-1) )
          return y->Elt(size-1);

       if ( xval <= x->Elt(0) )
          return y->Elt(0);
    }
 
	return linInterpol(x,xval,y);

    // Search the two point nearest xval
    // Recherche des points les plus proches (n encadrent pas forcement)
    // Search sup1 and sup2 as x(sup2) < x(sup1) < xval
    // equal to -1 if not find
    
/*    dist = K_HUGE_DOUBLE;          // distance x(sup1) from xval
    distprec = K_HUGE_DOUBLE;      // distance x(sup2) from xval
    sup1 = -1;
    sup2 = -1;

    for (index = 0; index < size; index++)
    {
	   valtmp = x->Elt(index);

       if (valtmp >= xval)
       { 
          if ( fabs(valtmp-xval) < dist )
          {
             sup2 = sup1;
             distprec = dist;

             sup1 = index;
             dist = fabs(valtmp-xval);
          }
          else if ( fabs(valtmp-xval) < distprec )
          {
             sup2 = index;
             distprec = fabs(valtmp-xval);
          }
       }
    } 

    // Search inf1 and inf2 as x(inf2) > x(inf1) > xval
    // equal to -1 if not find
    
    dist = K_HUGE_DOUBLE;          // distance x(inf1) from xval
    distprec = K_HUGE_DOUBLE;      // distance x(inf2) from xval
    inf1 = -1;
    inf2 = -1;
 
    for (index = 0; index < size; index++)
    {
	   valtmp = x->Elt(index);

       if ( valtmp <= xval )
       { 
          if ( fabs(valtmp-xval) < dist )
          {
             inf2 = inf1;
             distprec = dist;

             inf1 = index;
             dist = fabs(valtmp - xval);
          }
          else if ( fabs(valtmp-xval) < distprec )
          {
              inf2 = index;

              distprec = fabs(valtmp - xval);
          }
       }
    }

    // cas xval > x->Elt(size)

    if ( inf1 == -1 )
    {
       inf1 = sup2;
    }

    // cas xval < x->Elt(0)

    if ( sup1 == -1 )
    {
       sup1 = inf2;
    }

    //  do linear interpolation

    double x1 = x->Elt(inf1);
    double x2 = x->Elt(sup1);

    double y1 = y->Elt(inf1);
    double y2 = y->Elt(sup1);
 
    if ( x1 == x2 )
    {
       if ( y1 == y2 )
          return(y1);
       else 
          throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
                             "Interpolation failed  -  Bad data");
    }

    if ( xval == x1 )
       return(y1);
     
    if ( xval == x2 )
       return(y2);
     
    lambda = (x2-xval)/(x2-x1);

    yval = lambda*y1+(1.0-lambda)*y2;
 
    return(yval);
*/}

double  linInterpol2Col(const ARM_Vector& x, double xval, const ARM_Matrix& yy, int col,bool step)
{  
    if ( x.GetSize() != yy.GetNumLines() )
    {
       throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
                       "The two vectors should have the same size");
    }
 
    // case Size = 1
    long size = x.GetSize();

    if ( size == 1 )
    {
       return(yy.Elt(0,col));
    }
 
    if (step)
    {
       if ( xval >= x.Elt(size-1) )
          return yy.Elt(size-1,col);

       if ( xval <= x.Elt(0) )
          return yy.Elt(0,col);
    }
 
	return linInterpolCol(x,xval,yy,col);
}

/*----------------------------------------------------------------------------*
    Returns first derivative dy/dt of linear interpolation of y = f(x) 
     at point xval. 
*----------------------------------------------------------------------------*/

double D1LinInterpol(ARM_Vector* x, double xval, ARM_Vector* y)
{
    double yp;
    int    index, size;
    

    if ( x->GetSize() != y->GetSize() )
    {
       throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
                       "The two vectors should have the same size");
    }
    
    // Calculate the first derivative

    size = x->GetSize() - 1;

    if ( xval <= (*x)[0] + K_DOUBLE_TOL )
    {
       yp = (-3.0*(*y)[0]+4.0*(*y)[1]-(*y)[2]) 
               /((double)(-3.0*(*x)[0]+4.0*(*x)[1]-(*x)[2]));
    }
    else if ( xval >= (*x)[size] - K_DOUBLE_TOL)
    {
       yp = (3.0*(*y)[size]-4.0*(*y)[size-1]+(*y)[size-2]) 
             /((double)(3.0*(*x)[size]+4.0*(*x)[size-1]-(*x)[size-2]));
    }
    else 
    {
       index = indexBeforeValue(x, xval);

       if ( xval == (*x)[index] )
       {
          yp = ((*y)[index+1]-(*y)[index-1]) 
                /((double) ((*x)[index+1]-(*x)[index-1]));
       }
       else 
       {
          yp = ((*y)[index+1]
               -(*y)[index])/((double) ((*x)[index+1]-(*x)[index]));
       }
    }

    return(yp);
}



/*----------------------------------------------------------------------------*
    Returns 2nd derivative of linear interpolation of y = f(x) at point xval. 
    Returns kHugeDouble if failed.
*----------------------------------------------------------------------------*/

double D2LinInterpol(ARM_Vector *x, double xval, ARM_Vector *y)
{
    double y2, d1, d2;
    int    index, size;


    if (x->GetSize() != y->GetSize())
    {
       throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
                             "The two vectors should have the same size");
    }
    
    // Calculate the 2nd derivative

    size = x->GetSize() - 1;

    if ( xval <= (*x)[0] - K_DOUBLE_TOL )
    {
       y2 = (-3.0*D1LinInterpol(x, (*x)[0], y) 
            +4.0*D1LinInterpol(x, (*x)[1], y)-D1LinInterpol(x, (*x)[2], y))
             /( (double) (-3.0*(*x)[0]+4.0*(*x)[1]-(*x)[2]));
    }
    else if ( xval >= (*x)[x->GetSize()-1] + K_DOUBLE_TOL)
    {
       y2 = (3.0*D1LinInterpol(x, (*x)[size], y) 
             - 4.0*D1LinInterpol(x, (*x)[size-1], y) 
             + D1LinInterpol(x, (*x)[size-2], y))
             / ( (double) (3.0*(*x)[size]+4.0*(*x)[size-1]-(*x)[size-2]));
    }
    else 
    {
       index = indexBeforeValue(x, xval);
       d1 = (*x)[index] - (*x)[index-1];
       d2 = (*x)[index+1] - (*x)[index];

       y2 = (d1*(*y)[index+1]+d2*(*y)[index-1]-d1*d2*(*y)[index])
                     /(0.5*d1*d2*(d1+d2));
    }

    return(y2);
}



/*----------------------------------------------------------------------------*
    Returns value of linear interpolation of z = f(x,y) at point (xval,yval).
    z is stored as a vector and z(i,j) = f(x(i), y(j)) is defined as
    z->Elt(i*y->GetSize() + j).
    Returns kHugeDouble if failed.
*----------------------------------------------------------------------------*/

double twoDimLinInterpol(ARM_Vector* x, ARM_Vector* y, double xval, 
                         double yval, ARM_Vector *z)
{
    double lambda1, lambda2, val;
    int    i1, i2, j1, j2, sizex, sizey;
    
    sizex = x->GetSize();
    sizey = y->GetSize();
    

    if ( z->GetSize() != sizex * sizey )
    {
       throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
                             "Invalid z size");
    }

    if ( xval < (*x)[0] || xval > (*x)[sizex-1] ||
            yval < (*y)[0] || yval > (*y)[sizey-1] )
    {
       throw Exception(__LINE__, __FILE__, ERR_CONDITION_NOT_MEET,
             "x and/or y do not belong to the interpollation region");
    }

    i1 = indexBeforeValue(x, xval);
    j1 = indexBeforeValue(y, yval);
    i2 = i1 + 1;
    j2 = j1 + 1;
    
    
    if (xval == (*x)[i1] && yval == (*y)[j1]) 
    {
       val = (*z)[j1 + i1 * sizey];
    } 
    else if (xval == (*x)[i1] && yval != (*y)[j1]) 
    {
       lambda2 = ((*y)[j2] - yval) / ((*y)[j2] - (*y)[j1]);

       val = lambda2 * (*z)[j1 + i1 * sizey]
              + (1.0-lambda2) * (*z)[j2 + i1 * sizey];
    }
    else if (xval != (*x)[i1] && yval == (*y)[j1]) 
    {
       lambda1 = ((*x)[i2] - xval) / ((*x)[i2] - (*x)[i1]);

       val = lambda1 * (*z)[j1 + i1 * sizey]
              + (1.0-lambda1) * (*z)[j1 + i2 * sizey];
    }
    else 
    {
       lambda1 = ((*x)[i2] - xval) / ((*x)[i2] - (*x)[i1]);
       lambda2 = ((*y)[j2] - yval) / ((*y)[j2] - (*y)[j1]);

       val = lambda1 * lambda2 * (*z)[j1 + i1 * sizey]
          + lambda1 * (1.0-lambda2) * (*z)[j2 + i1 * sizey]
          + (1.0-lambda1) * lambda2 * (*z)[j1 + i2 * sizey]
          + (1.0-lambda1) * (1.0-lambda2) * (*z)[j2 + i2 * sizey];
    }

    return(val);
}

double linInterpol(double x, double x1, double y1, double x2, double y2)
{
        if (x1 == x2) return y1;
 
        return ((x2 - x)*y1 + (x - x1)*y2)/(x2 - x1);
}



// Triangular interpolation
double triangularInterpol(ARM_Vector* x, ARM_Vector* y,
                          ARM_Matrix* z, double xval, 
                          double yval)
{
    if (( x->GetSize() != z->GetNumLines() ) 
        || 
        ( y->GetSize() != z->GetNumCols() )
       )
    {
       throw Exception(__LINE__, __FILE__, ERR_SQUARE_OR_SIZE_PB,
                       "Interpolation failed  -  Bad data");
    }

    double zval = 0.;

    long xSize = x->GetSize();
    long ySize = y->GetSize();

    long xindex, yindex;

    ARM_Vector* vRowOrCol;

    if ( xval < x->Elt(0)-K_DOUBLE_TOL )
    {
       xindex = -1;
    }
    else if ( xval > x->Elt(xSize -1)+K_DOUBLE_TOL )
    {
       xindex = xSize-1;
    }
    else
    {
       xindex = indexBeforeValue(x, xval);
    }

    if ( yval < y->Elt(0)-K_DOUBLE_TOL )
    {
       yindex = -1;
    }
    else if ( yval > y->Elt(ySize-1)+K_DOUBLE_TOL )
    {
       yindex = ySize-1;
    }
    else
    {
       yindex = indexBeforeValue(y, yval);
    }

    if ( xindex == -1 )
    {
       if ( yindex == -1 )
       {
          return z->Elt(0,0);
       }
       else if ( yindex == ySize-1 )
       {
          return z->Elt(0,ySize-1);
       }
       else
       {
          vRowOrCol = z->GetRow(0);

          zval = linInterpol(y, yval, vRowOrCol);

          if (vRowOrCol)
          {
             delete vRowOrCol;

             vRowOrCol = NULL;
          }

          return zval;
       }
    }
    else if ( xindex == xSize-1 )
    {
       if ( yindex == -1 )
       {
          return z->Elt(xSize-1,0);
       }
       else if ( yindex == ySize-1 )
       {
          return z->Elt(xSize-1,ySize-1);
       }
       else
       {
          vRowOrCol = z->GetRow(xSize-1);

          zval = linInterpol(y, yval, vRowOrCol);

          if (vRowOrCol)
          {
             delete vRowOrCol;

             vRowOrCol = NULL;
          }

          return zval;
       }
    }
    else
    {
       if ( yindex == -1 )
       {
          vRowOrCol = z->GetColumn(0);

          zval = linInterpol(x, xval, vRowOrCol);

          if (vRowOrCol)
          {
             delete vRowOrCol;

             vRowOrCol = NULL;
          }

          return zval;
       }
       else if ( yindex == ySize-1 )
       {
          vRowOrCol = z->GetColumn(ySize-1);

          zval = linInterpol(x, xval, vRowOrCol);

          if (vRowOrCol)
          {
             delete vRowOrCol;

             vRowOrCol = NULL;
          }

          return zval;
       }
       else
       {
          double y3 = y->Elt(yindex)
                      +(y->Elt(yindex+1)-y->Elt(yindex))
                       /(x->Elt(xindex+1)-x->Elt(xindex))*(xval-x->Elt(xindex));

          double z3 = z->Elt(xindex,yindex)
                      +(z->Elt(xindex+1,yindex+1)
                      -z->Elt(xindex,yindex))/(x->Elt(xindex+1)
                      -x->Elt(xindex))*(xval-x->Elt(xindex));

          if ( y3 >= yval )
          {
             double lDeltaY = 0.;

             // lower triangle
             double z1 = z->Elt(xindex,yindex)
                         +(z->Elt(xindex+1,yindex)
                         -z->Elt(xindex,yindex))/(x->Elt(xindex+1)
                         -x->Elt(xindex))*(xval-x->Elt(xindex));

             if ( y3 != y->Elt(yindex) )
                lDeltaY = (z3-z1)/(y3-y->Elt(yindex));

             zval = z1+lDeltaY*(yval-y->Elt(yindex));
          }
          else
          {
             double lDeltaY = 0.;

             // upper triangle
             double z2 = z->Elt(xindex, yindex+1)+(z->Elt(xindex+1,yindex+1)
                         -z->Elt(xindex,yindex+1))/(x->Elt(xindex+1)
                         -x->Elt(xindex))*(xval-x->Elt(xindex));

             if ( y3 != y->Elt(yindex+1) )
                lDeltaY = (z3-z2)/(y3-y->Elt(yindex+1));

             zval = z2+lDeltaY*(yval-y->Elt(yindex+1));
          }

          return zval;
       }
    }
}



/*--------------------------------------------------------------------------*/
/*  Renvoie l'indice du 1er fwd intercepte par start ds la liste donnee     */
/*--------------------------------------------------------------------------*/
int FindNextFwdIdx(double calcDate,
                   double start,
                   ARM_Vector* ycFwdResetDate,
                   ARM_Vector* ycFwdStartDate,
                   double ycTerminalDate)
{
    int firstFwdIdx;
    int nbFwd = ycFwdStartDate->GetSize();

    if ( calcDate > ycFwdResetDate->Elt(nbFwd-1)+K_NEW_DOUBLE_TOL )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "No more forward alive at the event date");
    }

    for (firstFwdIdx = 0; firstFwdIdx < nbFwd; firstFwdIdx++)
    {
        if (calcDate - K_NEW_DOUBLE_TOL <= ycFwdResetDate->Elt(firstFwdIdx) &&
            start-K_NEW_DOUBLE_TOL <= ycFwdStartDate->Elt(firstFwdIdx))
        {
            if (start+K_NEW_DOUBLE_TOL < ycFwdStartDate->Elt(firstFwdIdx))
            {
                int prevIdx = MAX(0,firstFwdIdx-1);

                if ( calcDate-K_NEW_DOUBLE_TOL <= ycFwdResetDate->Elt(prevIdx) )
                {
                   // Le precedent fwd fixe ou fixera
                   firstFwdIdx = prevIdx;
                }
            }
            break;
        }
    }

    if (firstFwdIdx >= nbFwd)
    {
        if(start > ycTerminalDate + K_NEW_DOUBLE_TOL)
            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "Can't find any alive forward at the event date");
        else
            firstFwdIdx=nbFwd-1;

    }

    return firstFwdIdx;
}

/*--------------------------------------------------------------------------*/
/*  Interpolation des forwards directement a partir d'une liste de forwards */
/*  implicitement compatibles i.e. monetaire et de meme tenor               */
/*--------------------------------------------------------------------------*/
ARM_Matrix* InterpolateForwardRate( double calcDate,
                                    ARM_Vector* ycFwdResetDate,
                                    ARM_Vector* ycFwdStartDate,
                                    double ycTerminalDate,
                                    double ycRateConvert,
                                    ARM_Matrix* forward,
                                    ARM_Vector* rateStartDate)
{
    int nbRateDate=rateStartDate->GetSize();

    // Test de bonne utilisation
    if(rateStartDate->Elt(nbRateDate-1) > ycTerminalDate + K_NEW_DOUBLE_TOL)
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "Can't interpolate a forward rate starting after the terminal date");
    }

    int nbNode=forward->GetNumLines();
    int nbFwd=ycFwdResetDate->GetSize();

    ARM_Matrix* result = new ARM_Matrix(nbNode,nbRateDate,0.0);

    int nodeIdx,nextFwdIdx;
    double t,ratio;

    // Cherche du 1er forward intercepte
    int firstFwdIdx=FindNextFwdIdx(calcDate,rateStartDate->Elt(0),
        ycFwdResetDate,ycFwdStartDate,ycTerminalDate);

    // Interpolation directe des forwards
    int fwdIdx=nbFwd-1;
    for(int rateDateIdx=nbRateDate-1;rateDateIdx>=0;--rateDateIdx)
    {
        t=rateStartDate->Elt(rateDateIdx);
        while(t < ycFwdStartDate->Elt(fwdIdx) - K_NEW_DOUBLE_TOL)
        {
            // Positionnement par rapport aux forwards de diffusion
            if(fwdIdx > firstFwdIdx) --fwdIdx;
            else break;
        }

        if(ycFwdStartDate->Elt(nbFwd-1)-K_NEW_DOUBLE_TOL <= t)
        {
            // On est dans la zone du dernier forward
            for(nodeIdx=0;nodeIdx<nbNode;++nodeIdx)
                result->Elt(nodeIdx,rateDateIdx) = ycRateConvert *
                forward->Elt(nodeIdx,nbFwd-1);
        }
        else if(fwdIdx==firstFwdIdx &&
            t < ycFwdStartDate->Elt(firstFwdIdx) - K_NEW_DOUBLE_TOL)
        {
            // On est entre la sliceDate et la start du premier forward
            // i.e. sliceDate-eps <= t < Start(firstFwdIdx)-eps

            for(nodeIdx=0;nodeIdx<nbNode;++nodeIdx)
                result->Elt(nodeIdx,rateDateIdx) = ycRateConvert *
                forward->Elt(nodeIdx,fwdIdx);
        }
        else if(ycFwdStartDate->Elt(fwdIdx) - K_NEW_DOUBLE_TOL <= t &&
            t <= ycFwdStartDate->Elt(fwdIdx) + K_NEW_DOUBLE_TOL)
        {
            // On est juste sur une date de Start de fwd

            for(nodeIdx=0;nodeIdx<nbNode;++nodeIdx)
                result->Elt(nodeIdx,rateDateIdx) = ycRateConvert *
                forward->Elt(nodeIdx,fwdIdx);
        }
        else
        {
            // On est strictement entre deux forwards
            // i.e. Start(fwdIdx)+eps < t < Start(nextFwdIdx)-eps

            nextFwdIdx=fwdIdx+1;

            ratio=(t-ycFwdStartDate->Elt(fwdIdx))/
                (ycFwdStartDate->Elt(nextFwdIdx) -
                ycFwdStartDate->Elt(fwdIdx));

            for(nodeIdx=0;nodeIdx<nbNode;++nodeIdx)
            {
                result->Elt(nodeIdx,rateDateIdx) = ycRateConvert *
                    ( (1-ratio)*forward->Elt(nodeIdx,fwdIdx) +
                    ratio*forward->Elt(nodeIdx,nextFwdIdx) );
            }
        }
    } // for nb fwds a calculer

    return result;
}

/*--------------------------------------------------------------------------*/
/*  Interpolation des Zc a partir d'une liste de forwards et des savings    */
/*  account associes                                                        */
/*--------------------------------------------------------------------------*/
ARM_Matrix* InterpolateZc(double calcDate,
                          ARM_Vector* ycFwdResetDate,
                          ARM_Vector* ycFwdStartDate,
                          double ycTerminalDate,
                          int ycFwdDayCount,
                          ARM_Matrix* forward,
                          ARM_Matrix* savAcc,
                          ARM_Vector* zcDate,
						  void* miscData[]){return NULL;}
//{
//    int nbZcDate=zcDate->GetSize();
//
//    // Test de bonne utilisation
//    if(zcDate->Elt(nbZcDate-1) > ycTerminalDate + K_NEW_DOUBLE_TOL)
//    {
//       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
//            "Can't interpolate a Zc maturing after the terminal date");
//    }
//
//    // Recuperation des donnees additionnelles
//    ARM_ZeroCurve* zeroCurve = (ARM_ZeroCurve *) (miscData[0]);
//    ARM_Vector* spotForward = (ARM_Vector *) (miscData[1]);
//    ARM_Vector* forward0 = (ARM_Vector *) (miscData[2]);
//    ARM_Vector* shift = (ARM_Vector *) (miscData[3]);
//    ARM_Vector* proba = (ARM_Vector *) (miscData[4]);
//
//    bool isSpotForward = (spotForward!=NULL);
//    bool isRateAdjMethod = (zeroCurve!=NULL && forward0!=NULL && shift!=NULL);
//    bool isDriftAdjMethod = (isRateAdjMethod && proba!=NULL);
//    bool isSavingAccountMethod = (zeroCurve!=NULL && !isRateAdjMethod);
//
//    // Methode non implementee car la correction en taux sans recalage
//    // de la partie deterministe du drift est actuellement suffisante
//    isDriftAdjMethod=false;
//
//    double spotForward0=0.0;
//    if(forward0)
//        spotForward0=forward0->Elt(0);
//    double spotShift=0.0;
//    if(shift)
//        spotShift=shift->Elt(0);
//    double asOfDate=0.0;
//    if(zeroCurve)
//    {
//        asOfDate=zeroCurve->GetAsOfDate().GetJulian();
//        if(zeroCurve->GetCurrencyUnit()->GetLiborIndexDayCount() != ycFwdDayCount)
//        {
//           throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
//                "Inconsistency between forward rates and zero curve day count");
//        }
//    }
//
//    int nbNode=forward->GetNumLines();
//    int nbFwd=ycFwdResetDate->GetSize();
//
//    ARM_Matrix* result = new ARM_Matrix(nbNode,nbZcDate,0.0);
//
//    int nodeIdx,nextFwdIdx;
//    double t,delta,ratio,x,x0;
//    double interCoef,zcT,zcPrev,zcNext;
//    double startDate,nextStartDate,fwdShift,xShift;
//
//    // Cherche du 1er forward intercepte
//    int firstFwdIdx=FindNextFwdIdx(calcDate,zcDate->Elt(0),
//        ycFwdResetDate,ycFwdStartDate,ycTerminalDate);
//
//    double firstStartDate=ycFwdStartDate->Elt(firstFwdIdx);
//    ARM_Date firstStart(firstStartDate);
//    double spotStartDate = calcDate + firstStartDate-
//        ycFwdResetDate->Elt(firstFwdIdx);
//    if(firstStartDate < spotStartDate - K_NEW_DOUBLE_TOL)
//        spotStartDate=firstStartDate;
//
//    // Calcul des prix de Zc
//    int fwdIdx=nbFwd-1;
//    for(int zcDateIdx=nbZcDate-1;zcDateIdx>=0;--zcDateIdx)
//    {
//        t=zcDate->Elt(zcDateIdx);
//
//
//        // Positionnement par rapport aux forwards de diffusion
//        while(t < ycFwdStartDate->Elt(fwdIdx) - K_NEW_DOUBLE_TOL)
//        {
//            if(fwdIdx > firstFwdIdx) --fwdIdx;
//            else break;
//        }
//
//
//        if(ycTerminalDate-K_NEW_DOUBLE_TOL <= t)
//        {
//            // On est juste sur la date terminale
//            // i.e. T* - eps <= t <= T* + eps
//            for(nodeIdx=0;nodeIdx<nbNode;++nodeIdx)
//                result->Elt(nodeIdx,zcDateIdx)=1.0;
//        }
//
//        else if(ycFwdStartDate->Elt(nbFwd-1)-K_NEW_DOUBLE_TOL <= t)
//        {
//            if(t<=ycFwdStartDate->Elt(nbFwd-1)+K_NEW_DOUBLE_TOL)
//            {
//                // Juste sur la Start du dernier fwd
//                for(nodeIdx=0;nodeIdx<nbNode;++nodeIdx)
//                    result->Elt(nodeIdx,zcDateIdx)=savAcc->Elt(nodeIdx,nbFwd-1);
//            }
//            else
//            {
//                // On interpole car Start(nbFwd-1) + eps < t < T* - eps
//                delta=0.01*CountYears(ycFwdDayCount,t,ycTerminalDate);
//
//                if(isSavingAccountMethod)
//                {
//                    // La partie stochastique du saving account reste egale
//                    // a celle du dernier saving account connu
//                    zcT=zeroCurve->DiscountPrice((t-asOfDate)/K_YEAR_LEN);
//                    zcPrev=zeroCurve->DiscountPrice(
//                        (ycFwdStartDate->Elt(nbFwd-1)-asOfDate)/K_YEAR_LEN);
//                    interCoef=zcT/zcPrev;
//                    for(nodeIdx=0;nodeIdx<nbNode;++nodeIdx)
//                        result->Elt(nodeIdx,zcDateIdx)=interCoef*
//                            savAcc->Elt(nodeIdx,nbFwd-1);
//
//                }
//                else if(isRateAdjMethod)
//                {
//                    // Correction en taux ou en drift equivalente ds ce cas
//                    // car le taux periode t->T* est martingale sous Q*
//                    fwdShift=shift->Elt(nbFwd-1);
//                    x0 = zeroCurve->ForwardYield((ARM_Date)t,
//                        (ARM_Date)ycTerminalDate,K_COMP_PROP,K_ADJUSTED);
//                    xShift=fwdShift;
//                    interCoef=(x0+xShift) / (forward0->Elt(nbFwd-1)+fwdShift);
//
//                    for(nodeIdx=0;nodeIdx<nbNode;++nodeIdx)
//                    {
//                        x=(forward->Elt(nodeIdx,nbFwd-1)+fwdShift)*interCoef -
//                          xShift;
//                        result->Elt(nodeIdx,zcDateIdx) = 1+delta*x;
//                    }
//                }
//                else
//                {
//                    // Le tx fwd est prix egal au dernier connu
//                    for(nodeIdx=0;nodeIdx<nbNode;++nodeIdx)
//                    {
//                        result->Elt(nodeIdx,zcDateIdx)=
//                            (1.0 + delta*forward->Elt(nodeIdx,nbFwd-1));
//
//                    }
//                }
//            }
//        } // if sur ou apres la start de dernier fwd de diffusion
//
//        else if(fwdIdx==firstFwdIdx &&
//            t < firstStartDate - K_NEW_DOUBLE_TOL)
//        {
//            // On est entre la sliceDate et la start du premier forward
//            // i.e. sliceDate-eps <= t < Start(firstFwdIdx)-eps
//
//            delta=0.01*CountYears(ycFwdDayCount,t,firstStartDate);
//            if(isSavingAccountMethod)
//            {
//                // La partie stochastique du saving account reste egale
//                // a celle du prochain saving account connu
//                zcT=zeroCurve->DiscountPrice((t-asOfDate)/K_YEAR_LEN);
//                zcNext=zeroCurve->DiscountPrice((firstStartDate-asOfDate)/
//                    K_YEAR_LEN);
//                interCoef=zcT/zcNext;
//                for(nodeIdx=0;nodeIdx<nbNode;++nodeIdx)
//                {
//                    result->Elt(nodeIdx,zcDateIdx)=interCoef*
//                        savAcc->Elt(nodeIdx,firstFwdIdx);
//                }
//            }
//            else if(isSpotForward)
//            {
//                ratio=(t-spotStartDate)/(firstStartDate - spotStartDate);
//                if(t <= spotStartDate + K_NEW_DOUBLE_TOL)
//                {
//                    // On est avant ou sur le forward spot
//                    x0=0;/*zeroCurve->ForwardYield((ARM_Date)t,firstStart,
//                        K_COMP_PROP,K_ADJUSTED);*/
//                    if(isRateAdjMethod)
//                    {
//                        for(nodeIdx=0;nodeIdx<nbNode;nodeIdx++)
//                        {
//                            interCoef=(spotForward->Elt(nodeIdx)+spotShift)/
//                                        (spotForward0+spotShift);
//                            x = (x0+spotShift)*interCoef - spotShift;
//                            result->Elt(nodeIdx,zcDateIdx) = savAcc->
//                            Elt(nodeIdx,firstFwdIdx) * (1+delta*x);
//                        }
//                    }
//                    else
//                    {
//                        fwdShift=shift->Elt(firstFwdIdx);
//                        for(nodeIdx=0;nodeIdx<nbNode;nodeIdx++)
//                        {
//                            interCoef=(forward->Elt(nodeIdx,firstFwdIdx)+fwdShift)/
//                                      (forward0->Elt(firstFwdIdx)+fwdShift);
//                            x = (x0+fwdShift)*interCoef - fwdShift;
//                            result->Elt(nodeIdx,zcDateIdx) = savAcc->
//                            Elt(nodeIdx,firstFwdIdx) * (1+delta*x);
//                        }
//                    }
//                }
//                else if(isDriftAdjMethod)
//                {
//                    // Correction calee du drift du taux periode t->TnextFwdIdx
//                   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
//                        "Forward curve interpolation method not implemented");
//                }
//                else if(isRateAdjMethod)
//                {
//                    fwdShift=shift->Elt(firstFwdIdx);
//                    x0=0;/*zeroCurve->ForwardYield((ARM_Date)t,firstStart,
//                        K_COMP_PROP,K_ADJUSTED);*/
//                    xShift = (1-ratio)*spotShift + ratio*fwdShift;
//
//                    for(nodeIdx=0;nodeIdx<nbNode;++nodeIdx)
//                    {
//                        interCoef=
//                           (1-ratio)*log( (spotForward->Elt(nodeIdx)+spotShift)/
//                                         (spotForward0+spotShift) ) +
//                           ratio*log( (forward->Elt(nodeIdx,firstFwdIdx)+fwdShift)/
//                                      (forward0->Elt(firstFwdIdx)+fwdShift) );
//                        x = (x0+xShift)*exp(interCoef) - xShift;
//                        result->Elt(nodeIdx,zcDateIdx) = savAcc->
//                            Elt(nodeIdx,firstFwdIdx)*(1+delta*x);
//                    }
//                }
//                else
//                {
//                    // On interpole lineairement entre le forward spot et
//                    // le premier forward de diffusion connu
//                    for(nodeIdx=0;nodeIdx<nbNode;++nodeIdx)
//                    {
//                        x=spotForward->Elt(nodeIdx);
//                        result->Elt(nodeIdx,zcDateIdx) = savAcc->
//                            Elt(nodeIdx,firstFwdIdx)*(1+delta*
//                            (x + ratio*(forward->Elt(nodeIdx,firstFwdIdx)-x)));
//                    }
//                }
//            } // if un forward fixant a la slice existe
//
//            else
//            {
//                // Le tx fwd est prix egal au dernier connu
//                for(nodeIdx=0;nodeIdx<nbNode;++nodeIdx)
//                    result->Elt(nodeIdx,zcDateIdx)=savAcc->
//                        Elt(nodeIdx,firstFwdIdx)*(1+delta* forward->
//                        Elt(nodeIdx,firstFwdIdx));
//            }
//        } // if avant ou sur la start du premier forward de diffusion non fixe
//
//        else if(t <= ycFwdStartDate->Elt(fwdIdx) + K_NEW_DOUBLE_TOL)
//        {
//            // On est juste sur une date de Start de fwd
//            // i.e. Start(fwdIdx)-eps <= t <= Start(fwdIdx)+eps
//
//            for(nodeIdx=0;nodeIdx<nbNode;++nodeIdx)
//                result->Elt(nodeIdx,zcDateIdx)=savAcc->Elt(nodeIdx,fwdIdx);
//        }
//
//        else
//        {
//            // On est strictement entre deux forwards
//            // i.e. Start(fwdIdx)+eps < t < Start(nextFwdIdx)-eps
//
//            nextFwdIdx=fwdIdx+1;
//
//            startDate=ycFwdStartDate->Elt(fwdIdx);
//            nextStartDate=ycFwdStartDate->Elt(nextFwdIdx);
//
//            ratio=(t-startDate)/(nextStartDate - startDate);
//            delta=0.01*CountYears(ycFwdDayCount,t,nextStartDate);
//
//            if(isSavingAccountMethod)
//            {
//                // La partie stochastique du saving account reste egale
//                // a celle du prochain saving account connu
//                // Surtout ne pas interpoler entre les saving account encadrant
//                // car on va perdre la propriete martingale
//                zcT=0;/*zeroCurve->DiscountPrice((t-asOfDate)/K_YEAR_LEN);*/
//                zcNext=0;/*zeroCurve->DiscountPrice(
//                    (ycFwdStartDate->Elt(nextFwdIdx)-asOfDate)/K_YEAR_LEN);*/
//                interCoef=zcT/zcNext;
//                for(nodeIdx=0;nodeIdx<nbNode;++nodeIdx)
//                {
//                    result->Elt(nodeIdx,zcDateIdx)=interCoef*
//                        savAcc->Elt(nodeIdx,nextFwdIdx);
//                }
//            }
//            else if(isDriftAdjMethod)
//            {
//                // Correction calee du drift du taux periode t->TnextFwdIdx
//               throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
//                    "Forward curve interpolation method not implemented");
//            }
//            else if(isRateAdjMethod)
//            {
//                fwdShift=shift->Elt(fwdIdx);
//                double nextFwdShift=shift->Elt(nextFwdIdx);
//                x0=0;/*zeroCurve->ForwardYield((ARM_Date)t,
//                    (ARM_Date)nextStartDate,K_COMP_PROP,K_ADJUSTED);*/
//                xShift = (1-ratio)*fwdShift + ratio*nextFwdShift;
//
//                for(nodeIdx=0;nodeIdx<nbNode;++nodeIdx)
//                {
//                    interCoef=
//                       (1-ratio)*log( (forward->Elt(nodeIdx,fwdIdx)+fwdShift)/
//                                     (forward0->Elt(fwdIdx)+fwdShift) ) +
//                       ratio*log( (forward->Elt(nodeIdx,nextFwdIdx)+nextFwdShift)/
//                                  (forward0->Elt(nextFwdIdx)+nextFwdShift) );
//                    x = (x0+xShift)*exp(interCoef) - xShift;
//                    result->Elt(nodeIdx,zcDateIdx) = savAcc->
//                        Elt(nodeIdx,nextFwdIdx)*(1+delta*x);
//                }
//            }
//            else
//            {
//                // Le tx fwd est lineairement interpole entre
//                // les deux encardrant
//                for(nodeIdx=0;nodeIdx<nbNode;++nodeIdx)
//                {
//                    x=forward->Elt(nodeIdx,fwdIdx);
//                    result->Elt(nodeIdx,zcDateIdx)=savAcc->
//                        Elt(nodeIdx,nextFwdIdx)*(1+delta*
//                        (x + ratio*(forward->Elt(nodeIdx,nextFwdIdx)-x)));
//                }
//            }
//
//        } // if strictement entre 2 start de fwds de diffusion
//
//    } // for nb Zc a calculer
//
//    return result;
//}

/*--------------------------------------------------------------------------*/
/*  Calcul recursif des coefs d'un polynome de degre m passant par m+1 pts  */
/*--------------------------------------------------------------------------*/
ARM_Vector* PolynomialInterpolation(ARM_Vector* x,ARM_Vector* fx)
{
    int i,k,m,size;
    int nbValue=x->GetSize();

    ARM_Matrix *line = new ARM_Matrix(1,nbValue);
    for(k=0;k<nbValue;++k)
        line->Elt(0,k)=fx->Elt(k);

    ARM_Matrix *nextLine;

    double xi,xim,dx;
    for(m=1,size=nbValue-1;m<nbValue;++m,--size)
    {
        nextLine=new ARM_Matrix(m+1,size);
        for(i=0;i<size;++i)
        {
            xi=x->Elt(i);
            xim=x->Elt(i+m);
            dx=xi-xim;

            // k=0
            nextLine->Elt(0,i) = (line->Elt(0,i+1)*xi - line->Elt(0,i)*xim)/dx;

            // k=1 a m-1
            for(k=1;k<m;++k)
            {
                nextLine->Elt(k,i) = (line->Elt(k,i+1)*xi - line->Elt(k,i)*xim +
                    line->Elt(k-1,i) - line->Elt(k-1,i+1))/dx;
            }

            // k=m
            nextLine->Elt(m,i) = (line->Elt(m-1,i) - line->Elt(m-1,i+1))/dx;
        }
        delete line;
        line=nextLine;
    }

    ARM_Vector *result = new ARM_Vector(nbValue);
    for(k=0;k<nbValue;++k)
        result->Elt(k) = line->Elt(k,0);

    delete line;

    return result;
}



ARM_Vector* CalcSPLSecondDerivative(ARM_Vector* vect1, ARM_Vector* vect2)
{
    int i;

    double T_FACTOR = 0.9;

    double Ai;
    double Bi;
    double Ci;
    double Di;

    double Y2NullFlag = 0;

    double FirstFirstDeriv;
    double LastFirstDeriv;

    int sz = vect1->GetSize();

    ARM_Vector Io(sz);
    ARM_Vector S(sz);
    ARM_Vector SecondDeriv(sz);

    ARM_Vector SelectedX(sz);
    ARM_Vector SelectedY(sz);
 
    // Fill needed vectors
    int size = vect1->GetSize();
    for (i = 0; i < size; i++)
    {
        SelectedX[i] = (*vect1)[i];
        SelectedY[i] = (*vect2)[i];
    }

    // Compute the first derivatives
    FirstFirstDeriv = (SelectedY[1]-SelectedY[0])
                         /(SelectedX[1]-SelectedX[0])
                         -T_FACTOR*((SelectedY[2]-SelectedY[1])
                         /(SelectedX[2]-SelectedX[1])
                         -(SelectedY[1]-SelectedY[0])
                          /(SelectedX[1]-SelectedX[0]));

    LastFirstDeriv = (SelectedY[size-1]-SelectedY[size-2])
                        /(SelectedX[size-1]-SelectedX[size-2])
                        +T_FACTOR*((SelectedY[size-1]
                        -SelectedY[size-2])/(SelectedX[size-1]
                        -SelectedX[size-2])-(SelectedY[size-2]
                        -SelectedY[size-3])/(SelectedX[size-2]
                        -SelectedX[size-3]));

    if (Y2NullFlag )   
    {
       // "Natural spline"

       Io[0] = 0.0;
       S[0]  = 0.0;
    }
    else
    { 
       // "Clamped spline"

       Io[0] = -3.0/(SelectedX[1]-SelectedX[0])
                *(FirstFirstDeriv-((SelectedY[1]-SelectedY[0])
                                /(SelectedX[1]-SelectedX[0])));

       S[0] = -0.5 ;
    }

    // Calculate of S and I at each point
    for (i = 1; i < (size-1); i++)
    {
        Ai = SelectedX[i+1]-SelectedX[i];
        Bi = (SelectedX[i+1]-SelectedX[i-1])/3.0;
        Ci = SelectedX[i]-SelectedX[i-1];
        Di = (SelectedY[i+1]-SelectedY[i])/Ai
             -(SelectedY[i]-SelectedY[i-1])/Ci;

        Ai /= 6.0;
        Ci /= 6.0;

        Io[i] = Bi+Ci*S[i-1];
        S[i]  = -Ai/Io[i];
        Io[i] = (Di-Ci*Io[i-1])/Io[i];

    }

    // Calculate second deriv. at rank n
    if (Y2NullFlag )
    {
       SecondDeriv[size-1] = 0.0;
    }
    else
    {
       SecondDeriv[size-1] = 6.0/(SelectedX[size-1]
                                  -SelectedX[size-2])
                                  *(LastFirstDeriv-(SelectedY[size-1]
                                  -SelectedY[size-2])
                                   /(SelectedX[size-1]
                                   -SelectedX[size-2]))
                                   -Io[size-2];

       SecondDeriv[size-1] /= (2.0+S[size-2]);
    }

    // obtain all the second deriv.
    for (i = size-2; i >= 0; i--)
    {
        SecondDeriv[i] = Io[i]
                            +S[i]*SecondDeriv[i+1];
    }

    ARM_Vector* res = (ARM_Vector *) SecondDeriv.Clone();

    return(res);
}



/*----------------------------------------------------------------------------*

   This method implement either the normal spline interpolation or a specific
   2nd-derivative-keeping interpolation.
   According to the parameter keep2Der.

   N.B:
   ----

   The first time SplineInterpolateFunc has to be called
   with keep2Der = 1 

*----------------------------------------------------------------------------*/


void FindInXBounds(double InX,
                   ARM_Vector& VectX,
                   int& indInfX,
                   int& indSupX)
{
    int i;
    int size;
    int found = 0;
    double comparePrec = 1e-3;


    size = VectX.GetSize();

    for (i = 1; !found && i < size; i++)
    {
        indSupX = i ;

        if ( fabs(InX-VectX[indSupX]) <= comparePrec ) // Equality
        {
           indInfX = indSupX ;

           found = 1 ;
        }
        else
        {
           if ( InX < VectX[indSupX] )
           {
              indInfX = i-1;

              found = 1;
           }
        }
    }

    if (!found)
    {
       indSupX = size-1;
       indInfX = size-2;
    }
}



double SplineInterpolateFunc(ARM_Vector* vect1, ARM_Vector* vect2,
                             double inX,
                             ARM_Vector* SecondDerivCalc,
                             int keep2Der)
{
    int indInfX, indSupX;

    double res;

    double xi, yi, y2i;       // Points corresponding to indInfX of inX
    double xip1, yip1, y2ip1; // Points corresponding to indSupX of inX

    double dx, dy, y;

    double INIT_DX_VALUE = 1.0;

    ARM_Vector X = *vect1;
    ARM_Vector Y = *vect2;


    ARM_Vector SecondDeriv;

    if ( SecondDerivCalc == NULL )
    {
       ARM_Vector* res = CalcSPLSecondDerivative(vect1, vect2); 

       SecondDeriv = *res;

       delete res;
    }
    else
    {
       SecondDeriv = *SecondDerivCalc;
    }

    FindInXBounds(inX, X, indInfX, indSupX);

    if ( indInfX == indSupX )
    {
       res = Y[indInfX];
    }
    else
    {
       // Calculate (x,y,y2)

       xi  = X[indInfX];
       yi  = Y[indInfX];
       y2i = SecondDeriv[indInfX];

       xip1  = X[indSupX];
       yip1  = Y[indSupX];
       y2ip1 = SecondDeriv[indSupX];

       if (keep2Der)
       {
          // Clamped Spline: Specific Interpolation
          if ( inX < xi )
          {
             dx = INIT_DX_VALUE;

             y  = SplineInterpolateFunc(vect1, vect2, xi-dx/2.0,
                                        &SecondDeriv);

             dy = SplineInterpolateFunc(vect1, vect2, xi+dx/2.0,
                                        &SecondDeriv)-y;

             y = yi+dy/dx*(inX-xi)+0.5*y2i*(inX-xi)*(inX-xi);

             res = y ;

             return(res);
          }
          else if ( inX > xip1 ) // if above highest bound
          {
             dx = INIT_DX_VALUE;

             y  = SplineInterpolateFunc(vect1, vect2, xip1-dx/2.0,
                                        &SecondDeriv);

             dy = SplineInterpolateFunc(vect1, vect2, xip1+dx/2.0,
                                        &SecondDeriv)-y;

             y = yip1+dy/dx*(inX-xip1)+0.5*y2ip1*(inX-xip1)*(inX-xip1);

             res = y;

             return(res);
          }
       }

       // Normal Spline Interpolation
       double x;

       x = inX-xi;

       y = (y2ip1-y2i)/(6.0*(xip1-xi))*x*x*x;

       y = y+y2i/2.0*x*x;

       y = y+((yip1-yi)/(xip1-xi)-(xip1-xi)/6.0*(y2ip1+2.0*y2i))*x+yi;

       res = y;
    }

    return(res);
}



// Interpolation based on Lagrange Polynom N=3
double QuadraticLagrangeInterpol(double x,
                                 double x0, double y0,
                                 double x1, double y1,
                                 double x2, double y2)
{
    double Px; // The intrepolation result according to 
               // Lagrange's polynome of order N-1 (N equal 3)
    
    
    Px = (((x-x1)*(x-x2))/((x0-x1)*(x0-x2)))*y0
        +(((x-x0)*(x-x2))/((x1-x0)*(x1-x2)))*y1
        +(((x-x0)*(x-x1))/((x2-x0)*(x2-x1)))*y2;

    return(Px);
}



// Interpolation based on Lagrange Polynom N=4
double LagrangeInterpol(double x,
                        double x0, double y0,
                        double x1, double y1,
                        double x2, double y2,
                        double x3, double y3)
{
    double Px; // The intrepolation result according to 
               // Lagrange's polynome of order N-1 (N equal 4)
    
    
    Px =  (((x-x1)*(x-x2)*(x-x3))/((x0-x1)*(x0-x2)*(x0-x3)))*y0
         +(((x-x0)*(x-x2)*(x-x3))/((x1-x0)*(x1-x2)*(x1-x3)))*y1
         +(((x-x0)*(x-x1)*(x-x3))/((x1-x0)*(x1-x2)*(x1-x3)))*y2
         +(((x-x0)*(x-x1)*(x-x2))/((x3-x0)*(x3-x1)*(x3-x2)))*y3;

    return(Px);
}





/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
