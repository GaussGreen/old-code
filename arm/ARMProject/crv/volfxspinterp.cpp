/*
 * $Log: volfxspinterp.cpp,v $
 * Revision 1.7  2004/01/12 17:59:32  mab
 * if ( fabs(Xsup-Xinf) <= 1e-6 ) // Equality : NOW
 *
 * Revision 1.6  2004/01/12 17:57:19  mab
 * if ( fabs(Xsup-Xinf) <= 1e-2 ) // Equality
 * replaced by: if ( fabs(Xsup-Xinf) <= 1e-4 ) // Equality
 *
 * Revision 1.5  2003/11/19 11:14:01  mab
 * suppress const in InterpolInStrikeFwdTime
 * new call : res = SplineInterpolate(delta, !itsY2NullFlag);
 * instead of : res = SplineInterpolate(delta, itsY2NullFlag) ;
 *
 * Revision 1.4  2003/11/17 16:33:20  mab
 * Improvements in bounds detection
 *
 * Revision 1.3  2003/11/17 10:39:58  mab
 * Added : AddATMVolToSmile() in constructor
 *
 */

/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : volfxspinterp.cpp                                            */
/*                                                                            */
/* DESCRIPTION : Methods for the ARM_FXVolSmileInterpol class, a class        */
/*               dealing with spline FX Vol Interpolation                     */
/*                                                                            */
/* DATE        : Nov. 2003                                                    */
/*                                                                            */
/*----------------------------------------------------------------------------*/


#include "volfxspinterp.h"
#include "volcube.h"
#include "gaussian.h"


#define INIT_DX 1




/*----------------------------------------------------------------*/



void ARM_FXVolSmileInterpol::Init(void)
{
    itsFirstFirstDeriv = 0.0;

    itsLastFirstDeriv  = 0.0;

    itsY2NullFlag      = 1; // Natural spline
}



ARM_FXVolSmileInterpol::ARM_FXVolSmileInterpol(void)
{
    Init();
}



ARM_FXVolSmileInterpol::ARM_FXVolSmileInterpol(ARM_VolCurve* fxVol)
{
    Init();


    if ( fxVol->GetName() == ARM_VOL_CUBE )
    {
       ARM_VolCurve* ATMFxVol  = ((ARM_VolCube *) fxVol)->GetATMVol();

       ARM_VolLInterpol* smile = (ARM_VolLInterpol *)
                                  ((ARM_VolCube *) fxVol)->GetNthTenorVolMatrix(0);

       ARM_Matrix* smileMatrix = smile->GetVolatilities(); 

       ARM_Vector* matus     = ATMFxVol->GetExpiryTerms();
       ARM_Vector* deltaCols = smile->GetStrikes(); 
       ARM_Vector* ATMFxVolsVect = ATMFxVol->GetVolatilities()->GetColumn(0);
       

       itsATMVols   = *ATMFxVolsVect;
       delete ATMFxVolsVect;
       itsMatus     = *matus;
       itsDeltaCols = *deltaCols;

       itsSmile     = *smileMatrix;
    }
    else // Fx Vol. Matrix: we do the assumption that ATM vol is in col 0
    {
       ARM_Matrix* ATMSmileMatrix = fxVol->GetVolatilities(); 

       ARM_Vector* matus     = fxVol->GetExpiryTerms();
       ARM_Vector* deltaCols = ((ARM_VolLInterpol * ) fxVol)->GetStrikes(); 
       ARM_Vector* ATMFxVolsVect  = ATMSmileMatrix->GetColumn(0);
       

       itsATMVols   = *ATMFxVolsVect;
       delete ATMFxVolsVect;

       itsMatus     = *matus;

       // Copy the delta columns ordinates
       int colSize  = deltaCols->GetSize();
       itsDeltaCols.Resize(colSize-1);

       int i, j;

       // the first cell is reserved for ATM vol Column
       for (i = 1; i < colSize; i++)
       {
           itsDeltaCols[i-1] = deltaCols->Elt(i); 
       }

       // Copy the smile

       int nbMatus = matus->GetSize();

       itsSmile.Resize(nbMatus, colSize-1);

       for (i = 0; i < nbMatus; i++)
       {
           for (j = 0; j < (colSize-1); j++)
           { 
               itsSmile.Elt(i, j) = ATMSmileMatrix->Elt(i, j+1);
           } 
       }
    }

    itsCurSelectedDelta.Resize(itsDeltaCols.GetSize());
    itsCurSelectedVol.Resize(itsDeltaCols.GetSize());

    itsSecondDeriv.Resize(itsDeltaCols.GetSize());
    itsIo.Resize(itsDeltaCols.GetSize());
    itsS.Resize(itsDeltaCols.GetSize());

    AddATMVolToSmile();
}



void ARM_FXVolSmileInterpol::BitwiseCopy(const ARM_Object* srcVolCurve)
{
    ARM_FXVolSmileInterpol* vFxCurve = (ARM_FXVolSmileInterpol *) srcVolCurve; 


    itsMatus     = vFxCurve->itsMatus;
    itsDeltaCols = vFxCurve->itsDeltaCols;
    itsATMVols   = vFxCurve->itsATMVols;
    itsSmile     = vFxCurve->itsSmile;
    itsFirstFirstDeriv  = vFxCurve->itsFirstFirstDeriv;
    itsLastFirstDeriv   = vFxCurve->itsLastFirstDeriv;
    itsCurSelectedDelta = vFxCurve->itsCurSelectedDelta;
    itsCurSelectedVol   = vFxCurve->itsCurSelectedVol;
    itsSecondDeriv      = vFxCurve->itsSecondDeriv;
    itsIo               = vFxCurve->itsIo;
    itsS                = vFxCurve->itsS;
    itsY2NullFlag       = vFxCurve->itsY2NullFlag; 
}



/*----------------------------------------------------------------------------*/
/* Constructor(Copy)                                                          */
/*----------------------------------------------------------------------------*/

ARM_FXVolSmileInterpol::ARM_FXVolSmileInterpol(const ARM_FXVolSmileInterpol& volCurve) 
                       : ARM_Object(volCurve)
{
    Init();

    BitwiseCopy(&volCurve);
}



/*----------------------------------------------------------------------------*/
/* Assignment operator                                                        */
/*----------------------------------------------------------------------------*/

ARM_FXVolSmileInterpol& ARM_FXVolSmileInterpol::operator = 
                                           (const ARM_FXVolSmileInterpol& volCurve)
{
    (*this).ARM_Object::operator =(volCurve);

    BitwiseCopy(&volCurve);

    return (*this);
}



double ARM_FXVolSmileInterpol::InterpolDeltaIndex(const int indDate,
                                                  const double delta, 
                                                  const int y2Null)
{
    double    res ;

    // Test cas simple Interpolation sur la premiere ligne.
    itsY2NullFlag   = y2Null ;    


    // 2nd derivatives calculation
    CalculateSecondDerivative(indDate);

    res = SplineInterpolate(delta, !itsY2NullFlag);

    return(res);
}



/*----------------------------------------------------------------------------*

   This method implement either the normal spline interpolation or a specific
   2nd-derivative-keeping interpolation.
   According to the parameter keep2Der.

*----------------------------------------------------------------------------*/

double ARM_FXVolSmileInterpol::SplineInterpolate(const double delta,
                                                 const int keep2Der)
{
    int indInfDelta, indSupDelta;

    double res;

    double xi, yi, y2i;       // Points corresponding to indInfDelta of delta
    double xip1, yip1, y2ip1; // Points corresponding to indSupDelta of delta 

    double dx, dy, y;


    FindDeltaBounds(indInfDelta, indSupDelta, delta);        

    if ( indInfDelta == indSupDelta ) 
    {
       res = itsCurSelectedVol[indInfDelta] ;
    }
    else
    {
       // Calculate (x,y,y2)

       xi  = itsCurSelectedDelta[indInfDelta];
       yi  = itsCurSelectedVol[indInfDelta];
       y2i = itsSecondDeriv[indInfDelta];
                
       xip1  = itsCurSelectedDelta[indSupDelta];
       yip1  = itsCurSelectedVol[indSupDelta];
       y2ip1 = itsSecondDeriv[indSupDelta];

       if (keep2Der)
       {
          // Clamped Spline: Specific Interpolation
          if ( delta < xi )
          {
             dx = INIT_DX;

             y = SplineInterpolate(xi-dx/2.0);

             dy = SplineInterpolate(xi+dx/2.0)-y;

             y = yi+dy/dx*(delta-xi)+0.5*y2i*(delta-xi)*(delta-xi);

             res = y ;

             return(res);
          }
          else if ( delta > xip1 ) // if above highest bound
          {
             dx = INIT_DX;

             y = SplineInterpolate(xip1-dx/2.0);

             dy = SplineInterpolate(xip1+dx/2.0)-y;

             y = yip1+dy/dx*(delta-xip1)+0.5*y2ip1*(delta-xip1)*(delta-xip1);

             res = y;

             return(res);
          }
       }

       // Normal Spline Interpolation
       double x;

       x = delta-xi;

       y = (y2ip1-y2i)/(6.0*(xip1-xi))*x*x*x;

       y = y+y2i/2.0*x*x;

       y = y+((yip1-yi)/(xip1-xi)-(xip1-xi)/6.0*(y2ip1+2.0*y2i))*x+yi;

       res = y;
    }

    return(res);        
}



void ARM_FXVolSmileInterpol::FindDeltaBounds(int& indInfDelta,
                                             int& indSupDelta,
                                             const double delta)
{
    int i;
    int size; 
    int found = 0;
    double comparePrec = 1e-3;


    size = itsDeltaCols.GetSize();

    for (i = 1; !found && i < size; i++)
    {
        indSupDelta = i ;

        if ( fabs(delta-itsDeltaCols[indSupDelta]) <= comparePrec ) // Equality
        {
           indInfDelta = indSupDelta ;

           found = 1 ;
        }
        else
        {
           if ( delta < itsDeltaCols[indSupDelta] )
           {
              indInfDelta = i-1;

              found = 1;
           }
        }
    }

    if (!found)
    {
       indSupDelta = size-1;
       indInfDelta = size-2;
    }
}



double ARM_FXVolSmileInterpol::InterpolInStrikeFwdTime(double Forward, 
                                                       double Strike,
                                                       double matu, 
                                                       double Precision,
                                                       double sigmaATMF,
                                                       int y2Null)
{
    double sigma, prevSigma, delta ;

    int nbIterMax = 100; 
    int nbIter;

    sigma  = sigmaATMF/100.0;
    nbIter = 0;

    do
    {
        prevSigma = sigma ;

        // Delta Fwd
        delta = GetBSDelta(Forward, Strike, sigma, matu)*100.0;

        sigma = InterpolDeltaMatu(matu, delta, y2Null)/100.0;

        nbIter++;
    }
    while (( fabs(sigma-prevSigma) > Precision )
           &&
           ( nbIter < nbIterMax)
          );

    return(sigma*100.0);
}



double ARM_FXVolSmileInterpol::InterpolDeltaMatu(const double matu,
                                                 const double delta,
                                                 const int y2Null )
{
    double res ;
    double deltaInf, deltaSup;

    int _indInfDate, _indSupDate;

    // Retrieve bounds (inf,sup) for the maturity
    FindMatuBounds(_indInfDate, _indSupDate, matu);

    // Retrieve the delta on index _indInfDate 
    deltaInf = InterpolDeltaIndex(_indInfDate, delta, y2Null) ;

    
    // Retrieve the delta on index _indSupDate 
    deltaSup = InterpolDeltaIndex(_indSupDate, delta, y2Null) ;


    // Interpolate lineary 
    res = LinInterpol(itsMatus[_indInfDate],
                      itsMatus[_indSupDate],
                      deltaInf, deltaSup,
                      matu);

    return(res);
}



double ARM_FXVolSmileInterpol::LinInterpol(const double Xinf,
                                           const double Xsup, 
                                           const double Yinf,
                                           const double Ysup,
                                           const double X)
{
    double res;

    if ( fabs(Xsup-Xinf) <= 1e-6 ) // Equality
       res = Yinf ;
    else
       res = ((Xsup-X)*Yinf+(X-Xinf)*Ysup)/(Xsup-Xinf);

    return(res);
}



void ARM_FXVolSmileInterpol::FindMatuBounds(int& indInfDate,
                                            int& indSupDate,
                                            const double matu)        
{
    int i;
    int found = 0;
    double comparePrec = 1e-3;



    int size = itsMatus.GetSize();

    for (i = 1; !found && i < size ; i++)
    {
        indSupDate = i ;

        if ( fabs(matu-itsMatus[indSupDate]) <= comparePrec ) // Equality
        {
           indInfDate = indSupDate ;
           found = 1 ;
        }
        else
        {
           if ( matu < itsMatus[indSupDate] )
           {
              indInfDate = i-1;
              found = 1 ;
           }
        }
    }

    if (!found)
    {
       indSupDate = size-1;
       indInfDate = size-2;
    }
}



void ARM_FXVolSmileInterpol::CalculateSecondDerivative(const int indDate)
{
    int i;


    double Ai;
    double Bi;
    double Ci;
    double Di;


    // Fill needed vectors 
    int size = itsDeltaCols.GetSize();
    for (i = 0; i < size; i++)
    {
        itsCurSelectedDelta[i] = itsDeltaCols[i];
        itsCurSelectedVol[i]   = itsSmile.Elt(indDate, i);
    }                          

    // Compute the first derivatives
    itsFirstFirstDeriv = (itsCurSelectedVol[1]-itsCurSelectedVol[0])
                         /(itsCurSelectedDelta[1]-itsCurSelectedDelta[0])
                         -FACTOR_T*((itsCurSelectedVol[2]-itsCurSelectedVol[1])
                         /(itsCurSelectedDelta[2]-itsCurSelectedDelta[1])
                         -(itsCurSelectedVol[1]-itsCurSelectedVol[0])
                          /(itsCurSelectedDelta[1]-itsCurSelectedDelta[0]));    

    itsLastFirstDeriv = (itsCurSelectedVol[size-1]-itsCurSelectedVol[size-2])
                        /(itsCurSelectedDelta[size-1]-itsCurSelectedDelta[size-2])
                        +FACTOR_T*((itsCurSelectedVol[size-1]
                        -itsCurSelectedVol[size-2])/(itsCurSelectedDelta[size-1]
                        -itsCurSelectedDelta[size-2])-(itsCurSelectedVol[size-2]
                        -itsCurSelectedVol[size-3])/(itsCurSelectedDelta[size-2]
                        -itsCurSelectedDelta[size-3]));

    if (itsY2NullFlag )    
    { 
       // "Natural spline"

       itsIo[0] = 0.0;
       itsS[0]  = 0.0;
    }
    else
    {  
       // "Clamped spline"

       itsIo[0] = -3.0/(itsCurSelectedDelta[1]-itsCurSelectedDelta[0])
                *(itsFirstFirstDeriv-((itsCurSelectedVol[1]-itsCurSelectedVol[0])
                                /(itsCurSelectedDelta[1]-itsCurSelectedDelta[0])));

       itsS[0] = -0.5 ;
    }

    // Calculate of S and I at each point
    for (i = 1; i < (size-1); i++)
    {
        Ai = itsCurSelectedDelta[i+1]-itsCurSelectedDelta[i];
        Bi = (itsCurSelectedDelta[i+1]-itsCurSelectedDelta[i-1])/3.0;
        Ci = itsCurSelectedDelta[i]-itsCurSelectedDelta[i-1];
        Di = (itsCurSelectedVol[i+1]-itsCurSelectedVol[i])/Ai
             -(itsCurSelectedVol[i]-itsCurSelectedVol[i-1])/Ci;

        Ai /= 6.0;
        Ci /= 6.0;

        itsIo[i] = Bi+Ci*itsS[i-1];
        itsS[i]  = -Ai/itsIo[i];
        itsIo[i] = (Di-Ci*itsIo[i-1])/itsIo[i];

    }

    // Calculate second deriv. at rank n
    if (itsY2NullFlag )
    {
       itsSecondDeriv[size-1] = 0.0;
    }
    else
    {
       itsSecondDeriv[size-1] = 6.0/(itsCurSelectedDelta[size-1]
                                  -itsCurSelectedDelta[size-2])
                                  *(itsLastFirstDeriv-(itsCurSelectedVol[size-1]
                                  -itsCurSelectedVol[size-2])
                                   /(itsCurSelectedDelta[size-1]
                                   -itsCurSelectedDelta[size-2]))
                                   -itsIo[size-2];

       itsSecondDeriv[size-1] /= (2.0+itsS[size-2]);
    }

    // obtain all the second deriv. 
    for (i = size-2; i >= 0; i--)
    {
        itsSecondDeriv[i] = itsIo[i]
                            +itsS[i]*itsSecondDeriv[i+1];
    }
}


void ARM_FXVolSmileInterpol::AddATMVolToSmile(void)
{
    int i, j;

    int size1 = itsMatus.GetSize();
    int size2 = itsDeltaCols.GetSize();

    for (i = 0; i < size1; i++)
    {
        for (j = 0 ; j < size2; j++)
        {
            itsSmile.Elt(i, j) = itsSmile.Elt(i, j)+itsATMVols[i];
        }
    }
}



ARM_FXVolSmileInterpol::~ARM_FXVolSmileInterpol()
{
}



double ARM_FXVolSmileInterpol::GetBSDelta(const double fwd,
                                          const double strike, 
                                          const double vol,
                                          const double T,
                                          int CallPut)
{
    double delta ;

    if ( T <= 0.0 )
       delta = 0.0 ;
    else
    {
       // Optimization

       double vsqrt = vol*sqrt(T);

       delta = CallPut*cdfNormal(CallPut*(log(fwd/strike)/vsqrt+0.5*vsqrt));
    }

    return(delta);
} 



/*----------------------------------------------------------------------------*/
/*---- End Of File ----*/
