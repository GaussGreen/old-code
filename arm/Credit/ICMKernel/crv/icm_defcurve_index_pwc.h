#error no longer part of the project 

#ifndef _ICM_DEFCURVE_INDEX_H
#define _ICM_DEFCURVE_INDEX_H

#include "ICMKernel\crv\icm_constant_piecewise.h"

/*********************************************************************************/
/*! \class  ICM_DefcurveIndex icm_defcurve_index_pwc.h.h "icm_defcurve_index_pwc.h.h"
 *  \author d. pouponneau vs jp riaudel
 *	\version 1.0
 *	\date   1 january 2006
 *	\file   icm_defcurve_index_pwc.h.h
 *	\brief  Default Curve */
/***********************************************************************************/

class ICM_DefcurveIndex : public ICM_Constant_Piecewise
{
    private:
		vector<double>	itsRefSpread;
		ICM_Pricer*     itsSltPricer;	// associated during calibration. 

    public:
		ICM_DefcurveIndex(const ICM_DefcurveIndex&);  
		ICM_DefcurveIndex(const ARM_Date& asOf,
						 const vector<string>& terms,
						 ARM_Vector* rates,
						 ARM_Vector* Refrates,
						 double& Recovery,
						 ARM_ZeroCurve* zc,
						 int intRule,
						 int adjStartDate,
						 qCDS_ADJ adj /* = qCredit_Default */ ,  	
						 const std::string& ccy  ,
						 const std::string& label,
						 bool issummitcurve /*= true */,
						 const ARM_VolCurve* VolCurve /*= NULL */,
						 long PayFreq /*= K_QUARTERLY */,
						 qDEFCURVE_CALIB_ALGO calibAlgo,
						 const std::string& calibrationData,
						 int Lag)
		{
			Init();

			Set (asOf,terms,rates,Refrates,Recovery,zc,intRule,adjStartDate,adj,ccy,
				label,issummitcurve,VolCurve,PayFreq,calibAlgo,calibrationData,Lag);
		}

		void Set (const ARM_Date& asOf,
				 const vector<string>& terms,
                 ARM_Vector* rates,
				 ARM_Vector* Refrates,
				 double Recovery,
				 ARM_ZeroCurve* zc,
				 int intRule,
				 int adjStartDate,
				 qCDS_ADJ adj /* = qCredit_Default */ ,  	
                 const std::string& ccy  ,
				 const std::string& label/* = NULL */ ,
				 bool issummitcurve /* true */ ,
				 const ARM_VolCurve* VolCurve /* = NULL */,
				 long PayFreq /* = K_QUARTERLY */,
				 qDEFCURVE_CALIB_ALGO calibAlgo,
				 const std::string& calibData,
				 int Lag)
		{
			
			itsRefSpread.resize(Refrates->GetSize()+1);
			for (int i=1;i<Refrates->size()+1;i++) {itsRefSpread[i]=Refrates->Elt(i-1);}

			itsRefSpread[0] = Refrates->Elt(0);

			ICM_Constant_Piecewise::Set (asOf,terms,rates,Recovery,zc,intRule,adjStartDate,adj,ccy,
								label,issummitcurve,VolCurve,PayFreq,calibAlgo,calibData,Lag,ICM_Parameters());

		}

        ICM_DefcurveIndex(void) { Init();}

		void Init()
		{
			SetName(ICM_DEFCURVE_INDEX);
			itsRefSpread.clear();
			itsSltPricer=NULL;
		}

		void SetRefSpreads(ARM_Vector* RefSpread) 
		{itsRefSpread.resize(RefSpread->GetSize());
		for (int i=0;i<itsRefSpread.size();i++) {itsRefSpread[i]=RefSpread->Elt(i);}}

		void SetRefSpreads(const vector<double>&	RefSpread) {itsRefSpread=RefSpread;}
		// vector<double>& GetRefSpread() {return itsRefSpread;}

		//void View(char* id = NULL, FILE* ficOut = NULL ){}

// 		void BitwiseCopy(const ARM_Object* src)
// 		{	ICM_DefcurveIndex* dcindex = (ICM_DefcurveIndex*) src;
// 			itsRefSpread = dcindex->itsRefSpread;	}

//         void Copy(const ARM_Object* srcZc)
// 		{ICM_Constant_Piecewise::Copy(srcZc);
// 		 BitwiseCopy(srcZc);}

        virtual ARM_Object* Clone(void) ;
//         { ICM_Constant_Piecewise* theClone = new ICM_DefcurveIndex();
//           theClone->Copy(this);
//           return(theClone);}

	private:
		virtual void Calibrate();
	public:
		double Evaluate(const double& x);

		// inline double GetRefSpread(int i) { return itsRefSpread[i];}

		virtual ICM_DefaultCurve* GenerateShiftCurve(const vector<string>& terms, 
												  const ARM_Vector& epsilon ,
												  qSENSITIVITY_TYPE mode ) const ;

/** JLA Removed
	private:
		virtual ICM_DefaultCurve* xGenerateShiftCurve(double epsilon ,
												  qSENSITIVITY_TYPE mode ) const;
												  **/ 
	private:
		
		ICM_DefcurveIndex& operator=(const ICM_DefcurveIndex&); // NA
};


#endif /*---- End of file ----*/
