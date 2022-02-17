/*
 * icm_option.h
 * Author : FBenatig
 * September 2004
 *

/*----------------------------------------------------------------------------*
 
    icm_option.h
 
    Implements the ICM_Option class, a class to deal with options on CDS, CDS Index, Tranches
         
    deals with option priced with analytic valuation formulae in the B&S framework
 
*----------------------------------------------------------------------------*/
#ifndef _ICMOPTION_H
#define _ICMOPTION_H

#include "ICMKernel/inst/icm_security.h"
#include "ICMKernel/pricer/icm_pricer.h"

class ICM_Option : public ICM_Security
{
    private:

        double				itsStrike;		// option strike                                
        int                 itsOptionType;  // option type K_Call=1 or K_Put=-1
		qDEF_MAT			itsKoType;		// KO Type KO = 1 , No KO with Acceleration 0, No KO without Acceleration -1
		//ARM_Date			itsExpiryDate;	// expiry date of the option in ICM_Security itsendDateNA
		string			itsUnderlyingMaturityTerm ; // in Terms
		ARM_Date		itsUnderlyingMaturityDate ; // in Date
		double			itsNotional;	// notional * quantity (-1 /1)

    public:

        ICM_Option(void)
        {
            Init();
        }
		ICM_Option( const string&  underMaturity,
					const ARM_Date& underMaturityDate,
					   const ARM_Date&	ExpiryDate,
					   const string&	ccy,
					   const qCDS_ADJ cds_adj,
					   bool endAdj,
					   double strike,
					   int optionType,
				       qDEF_MAT KO,
					   double Notional);

		ICM_Option(const ARM_Date&  underMaturity,
					   const ARM_Date&	ExpiryDate,
					   double strike,
					   int optionType,
				       qDEF_MAT KO );

        ICM_Option(const ICM_Option &);
        ICM_Option& operator = (const ICM_Option& option);                            

		~ICM_Option(void){}

		void View(char* id = NULL, FILE* fOut = NULL);
		void BitwiseCopy(const ARM_Object* srcOption);
		void Copy(const ARM_Object* srcOption);
		ARM_Object* Clone(void);
       
		
		void Init();
		
		int IsCall(void) {return(itsOptionType == K_CALL);}
        int IsPut(void) {return(itsOptionType == K_PUT);}
        
        void SetStrike(double strike) { itsStrike = strike ; }
		double GetStrike(void) const { return(itsStrike); }
		// just Renamed
        void SetExpiry(const ARM_Date& Date) { this->SetEndDateNA(Date); }
		const ARM_Date& GetExpiry(void) const  { return(this->GetEndDateNA()); }

		void SetKoType(qDEF_MAT KO) {itsKoType = KO ;}
		qDEF_MAT GetKoType(void) const {return itsKoType ;}
		
		void SetOptionType(int optiontype){ itsOptionType = optiontype;}
		int GetOptionType(void){ return itsOptionType;}

		void SetUnderMaturityDate( const ARM_Date&  date){	itsUnderlyingMaturityDate = date;}
		const ARM_Date& GetUnderMaturityDate(void) { return itsUnderlyingMaturityDate; }
		
		void SetUnderMaturityTerm( const string&  date){	itsUnderlyingMaturityTerm = date;}
		const string& GetUnderMaturityTerm(void) { return itsUnderlyingMaturityTerm; }

		void SetNotional( double q) {itsNotional = q;}
		const double GetNotional() const {return itsNotional ;}
protected :
		 void Set( const ARM_Date&  undermaturity,	 
				 const ARM_Date& ExpiryDate,
				 double Strike,
				 int optiontype,
				 qDEF_MAT KO,
				 double Notional);
};


#endif
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/






