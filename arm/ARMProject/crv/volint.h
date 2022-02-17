/*
 * $Log: volint.h,v $
 * Revision 1.18  2004/03/26 16:49:52  mab
 * ApproxaimatedFwd now in : fromto.h
 *
 * Revision 1.17  2003/11/20 15:25:00  ebenhamou
 * remove this
 *
 * Revision 1.16  2003/11/14 11:26:33  mab
 * Added : int volType = K_ATMF_VOL in constructor
 *
 * Revision 1.15  2003/10/22 16:39:27  jpriaudel
 * modif pour Olivier
 *
 * Revision 1.14  2003/10/07 15:09:06  ebenhamou
 *  added const accessor for const correctness
 *
 * Revision 1.13  2003/09/25 10:35:42  mab
 * A flag added to : ConvertToBSVol()
 *
 * Revision 1.12  2003/09/23 09:33:21  mab
 * Added: ConvertToBSVol(..)
 *
 * Revision 1.11  2003/09/22 12:39:07  mab
 * added: const in constructors parameters
 *
 * Revision 1.10  2003/09/18 17:11:02  ebenhamou
 * change to const
 *
 * Revision 1.9  2003/09/16 15:57:38  mab
 * Added: AP addins :  ConvertToNormalVol()
 * double ImpliedLogNorIRVol()
 *
 * Revision 1.8  2003/08/18 09:49:15  mab
 * improvement in the principal method :
 * virtual double VolatilityFunction(double m1, double K, double m2)
 * virtual double VolatilityFunction(double m1, double m2) :
 * No default parameter but 2 methods!
 *
 * Revision 1.7  2003/02/11 14:59:55  mab
 * Added : void UpdateCol(ARM_Vector* pCol, double tenor);
 * ARM_Vector* GetCol(double tenor);
 *
 * Revision 1.6  2002/09/17 16:36:33  mab
 * Added : ARM_Vector* GetStrikes(void)
 * method
 *
 * Revision 1.5  2001/04/18 16:34:20  abizid
 * Ajout du View
 *
 * Revision 1.4  2001/04/03 11:59:14  nicolasm
 * Print.
 *
 * Revision 1.3  1999/09/07 08:16:18  nicolasm
 * Ajout champ volType dans le constructeur principal
 *
 * Revision 1.2  1999/04/14 13:28:30  mab
 * Rajout Log pour RCS
 *
 */


/*----------------------------------------------------------------------------*
  volint.h

    
  Header for the ARM_VolLInterpol class, a class for computing a ARM_VolCurve
 
*----------------------------------------------------------------------------*/ 
#ifndef _VOLINT_H
#define _VOLINT_H




#include "dates.h"
#include "volcurv.h"
#include "linalg.h"
#include "interpol.h"



#define K_MIN_NUM_YEAR_TERMS          3


class ARM_VolFlat;

class ARM_VolLInterpol : public ARM_VolCurve 
{
    private:
        
        ARM_Vector* itsStrikes;  // Strikes  (colomns of the matrix) 

        // Methods

	public:

        // here m2 is not relevant
        virtual double VolatilityFunction(double m1, double k, double m2);

        virtual double VolatilityFunction(double m1, double m2);

    private:

        virtual void Init(void);

    public:

        ARM_VolLInterpol(void);

        ARM_VolLInterpol(const ARM_Date& asOf, ARM_Vector* yearTerms, 
                         ARM_Vector* strikes, ARM_Matrix* volatilities,
                         int strikeType = K_STK_TYPE_PRICE,
                         int volType = K_ATMF_VOL,
                         ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);

        // Constructor for at the money vol curve 
        // (used eg for cap & floor ATM vol) 
        // -> itsStrikes contains 1 elt

        ARM_VolLInterpol(const ARM_Date& asOf, ARM_Vector* yearTerms, 
                         ARM_Vector* volatilities, 
                         int strikeType = K_STK_TYPE_PRICE,
                         int volType = K_ATMF_VOL,
                         ARM_Currency* ccy = ARM_DEFAULT_CURRENCY);

        ARM_VolLInterpol(ARM_VolFlat* volFlat, 
                         ARM_Vector* expiries = NULL,
                         ARM_Vector* undTenors = NULL);

        ARM_VolLInterpol(const ARM_VolLInterpol &);
    
        virtual ~ARM_VolLInterpol(void);

        ARM_VolLInterpol& operator = (const ARM_VolLInterpol &);

		virtual double CalcNumericalObjectSignature(void);

		void BitwiseCopy(const ARM_Object* srcVollint);
 
        void Copy(const ARM_Object* vollint);
 
        virtual ARM_Object* Clone(void);

        void View(char* id = NULL, FILE* fOut = NULL);

        void Set(int szyt, double* yt, int szstk, double* stk, 
                 double* vol, char* date, int strikeType=K_STK_TYPE_PRICE);
 
        void SetStrikes(ARM_Vector *strikes);

        ARM_Vector* GetStrikes(void);

		/// const version for const correctness!
        ARM_Vector* GetStrikes() const;

        void UpdateCol(ARM_Vector* pCol, double tenor);

        ARM_Vector* GetCol(double tenor);

        // Conversion de la vol forward Black-Scholes en vol Normale:
        // Rq: utilise ImpliedLogNorIRVol
        ARM_VolLInterpol* ConvertToNormalVol(ARM_ZeroCurve* YieldCurve,
                                             int IsSwaptionVol = 1, 
                                             int InPct = 1);

        // Interpolation dans la matrice de vol absolue et conversion
        // en vol forward Black-Scholes
        double ImpliedLogNorIRVol(double Mat, double Tenor,
                                  ARM_ZeroCurve* YieldCurve,
                                  int IsSwaptionVol = 1,
                                  int InPct = 1);

        ARM_VolLInterpol* ConvertToBSVol(ARM_ZeroCurve* YieldCurve,
                                         int IsSwaptionVol = 1,
                                         int InPct = 1,
                                         int Post100 = 0);
};

#endif
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
