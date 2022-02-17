/*
 * $Log: volcube.h,v $
 * Revision 1.22  2004/02/11 14:36:36  jpriaudel
 * BumpVolatility added
 *
 * Revision 1.21  2003/11/25 21:09:45  jpriaudel
 *  constructor added
 *
 * Revision 1.20  2003/11/14 11:20:53  mab
 * Take in account Splined FX VOL
 *
 * Revision 1.19  2003/09/25 10:06:45  ykhlif
 * correction de SetATMVol().
 *
 * Revision 1.18  2003/09/09 11:42:31  ebenhamou
 * move sparsevolcube to its own file
 *
 * Revision 1.17  2003/08/29 07:08:52  ebenhamou
 * added validation when filling vol in sparse vol cube
 *
 * Revision 1.15  2003/08/26 11:33:33  ebenhamou
 * added sparse vol cube
 *
 * Revision 1.13  2003/08/18 09:50:37  mab
 * improvement in the principal method :
 * virtual double VolatilityFunction(double m1, double K, double m2)
 * virtual double VolatilityFunction(double m1, double m2) :
 * No default parameter but 2 methods!
 *
 * Revision 1.12  2003/04/03 08:57:07  mab
 * Correction of View prototype
 *
 * Revision 1.11  2003/02/11 14:49:37  mab
 * Added : void ParallelShift(double value);
 *
 * Revision 1.10  2002/03/25 16:07:11  nicolasm
 * Ecriture de la methode Print
 *
 * Revision 1.9  2001/10/03 16:08:09  smysona
 * SetSTickyStrike
 *
 * Revision 1.8  2001/08/09 18:12:03  smysona
 * Correction bitwiseCopy
 *
 * Revision 1.6  2001/08/06 09:59:38  smysona
 * Ajout mode ATM
 *
 * Revision 1.5  2001/07/23 17:11:00  smysona
 *  Correction clone
 *
 * Revision 1.4  2001/07/23 12:35:12  mab
 * Rajout de view
 *
 * Revision 1.3  2001/07/20 12:52:09  mab
 * Rajout de : View()
 *
 * Revision 1.2  2001/07/18 13:55:28  mab
 * dos2unix
 *
 */

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*  Class to handle a Vol cube, ie a Vol depending on                         */
/*  three parameters                                                          */
/*  - moneyness                                                               */
/*  - expiry of the option                                                    */
/*  - tenor of the underluying                                                */
/*                                                                            */
/*  The cube is essentially composed by several layers of                     */
/*  vol matrices, one for each underluying.                                   */
/*                                                                            */
/*----------------------------------------------------------------------------*/


#ifndef _VOLCUBE_H
#define _VOLCUBE_H


#include <vector>

#include "armglob.h"
#include "volint.h"
//#include "frmutils.h"


/*!
 * A volcube is a vector of volcurve indexed by itsUnderlyings
 * usually itsUnderlying is the tenor of the underlying rate
 * hence for a one year rate, the value should be one
 * 
 * We can also define an ATM vol with the volCurve itsATMVol
 * but this is not compulsory. 
 * 
 * if you want to use ATM vol, then itsATMref = true
 * the design is that the ATM Vol gives the overall
 * level of vol and then the ARM_volcurve gives the smile 
 * on top of it
 *
 * if itsStickyDelta=true, then uses moneyness
 * otherwise itsStickyDelta=false, uses by strike
 * if using strike, then the matrix itsStrikeLevels has to be filled
 * with the row corresponding to the itsUnderlying
 * and the column to the strike
 * hence something like
 *
 *                    tenor\strike
 *                                 0.01   0.03   0.05
 *                        1y     ( *      *      *    )
 *    itsStrikeLevels =   5y     ( *      *      *    )
 *                        10y    ( *      *      *    )
 *
 */

class ARM_VolCube : public ARM_VolLInterpol         
{
    private :
    
        vector < ARM_VolCurve* > * itsVols;
    
        ARM_VolCurve*              itsATMVol;
    
        ARM_Vector*                itsUnderlyings;
    
        bool                       itsATMref;
    
        bool                       itsStickyDelta;
    
        ARM_Matrix*                itsStrikeLevels; // ATM strike reference for
        
        // stticky strike
        // underlying*expiry
       
        void Init(void);
    
        double VolatilityFunction(double expiry, 
                                  double moneyness, 
                                  double underlying);
    
        double VolatilityFunction(double expiry, 
                                  double underlying);
    
    public :
    
        ARM_VolCube(void);
    
        ARM_VolCube(vector<ARM_VolCurve*>* inVols, ARM_Vector* underlyings, 
                    const ARM_Date& lastKnownDate = ARM_Date() );
    
        ARM_VolCube(ARM_VolCurve** inVols, int size, ARM_Vector* tenors);

        ARM_VolCube(ARM_VolCurve* atmVol, ARM_VolCurve** inVols, 
                    int size, ARM_Vector* tenors,
                    int volType = K_ATMF_VOL,
					int checkCcy = 1);
    
        ARM_VolCube(const ARM_VolCube& volCube);
    
        virtual ~ARM_VolCube(void);

	   inline virtual ARM_Matrix* GetVolatilities() const
        {
           return(itsATMVol->GetVolatilities());
        }
    

virtual double CalcNumericalObjectSignature(void);

        double VolatilityFunctionByStrike(double expiry, double moneyness, 
                                          double underlying = 0.0);

		double ComputeCorrelByExpiry(double aExpiry, double aTenor1, double aTenor2);

		double ComputeSmileOnly(double expiry, double moneyness, double underlying);

//		double ComputeHyperCorrel(double aOptMat, double aTenor1,double aTenor2, double aStrike = 0.0);
    
    
        /***********************/
        /*    ARM Stuff        */
        /***********************/
    
        ARM_VolCube& operator = (const ARM_VolCube &);
    
        void BitwiseCopy(const ARM_Object* srcObject);
    
        void Copy(const ARM_Object* VolIn);
    
        virtual ARM_Object* Clone(void);
    
        void View(char* id = NULL, FILE* ficOut = NULL);
    
        void ParallelShift(double value); 
    
        ARM_VolCurve* GetNthTenorVolMatrix(int i);

		virtual ARM_VolCurve*	GetVolCurve(string& aTenor);
  
        virtual void BumpVolatility(double value,
                            int nthLine = 0,
                            int nthCol = 0,
                            int isCumul = K_NO,
						    int isAbsolute = K_YES);

        virtual void BumpSmile(double value,
							   double tenor,
							   int nthLine = 0,
							   int nthCol = 0,
							   int isCumul = K_NO,
							   int isAbsolute = K_YES);

        void ParallelShiftAll(double value);

        inline ARM_Vector* GetUnderLyings(void)
        {
            return itsUnderlyings;
        }
    
        inline vector < ARM_VolCurve* >* GetVols(void)
        {
            return itsVols;
        }
    
        void SetInterpType(int interpType)
	    {
		    itsATMVol->SetInterpType(interpType);
	    }

        inline ARM_VolCurve* GetATMVol(void)
        {
            return itsATMVol;
        }
    
        inline bool GetATMref(void)
        {
            return itsATMref;
        }
    
        inline bool GetStickyDelta(void)
        {
            return itsStickyDelta;
        }
    
        inline ARM_Matrix* GetStrikeLevels(void)
        {
            return itsStrikeLevels;
        }
    
        inline void SetUnderLyings(ARM_Vector* in)
        {
            itsUnderlyings = in;
        }
    
        inline void SetVols(vector < ARM_VolCurve* >* in)
        {
            itsVols = in;
        }
    
        inline void SetATMVol(ARM_VolCurve* inVol)
        {
            if (itsATMVol)
            {
               delete itsATMVol;
               
               itsATMVol = NULL;
            }

            if (inVol)
               itsATMVol = (ARM_VolCurve *) inVol->Clone();
        }
    
        inline void SetATMref(bool in)
        {
            itsATMref = in;
        }
    
        inline void SetStickyDelta(bool in)
        {
            itsStickyDelta = in;
        }
    
        inline void SetStrikeLevels(ARM_Matrix* in)
        {
            itsStrikeLevels = in;
        }
    
        void SwitchSwopt2StickyStrike(void);
};



#endif
/*----------------------------------------------------------------------------*/
/*---- End Of File ----*/