
#ifndef _ICMOPTION_RESTRIKABLE_H
#define _ICMOPTION_RESTRIKABLE_H

#include "ICMKernel/inst/icm_option.h"
#include "ICMKernel/pricer/icm_pricer.h"

class ICM_Option_Tranche : public ICM_Option
{
    private:

		double itsRehauss ;
		int itsExerciceType ;
		double itsObservationFreq ;
		int itsDiffCDO ;
		int itsIsCMSpread ;
		double itsCMSpreadMatu ;
		double itsInitSpread ;
		ICM_Security* itsUnderlying ; // always ICM_Mezz for restrikable

    public:

        ICM_Option_Tranche(void)
        {
            Init();
        }
		void Init(void);

		/*ICM_Option_Tranche( const ARM_Date& underMaturity,
						const ARM_Date& ExpiryDate,
						double strike,
						int optionType,
						ICM_Security* Underlying,
						double rehauss,
						double Freq,
						int ExecType,
						int DiffCDO)
		{
			Init();
			Set(underMaturity, ExpiryDate,strike, optionType, Underlying, rehauss, Freq, ExecType, DiffCDO);
		}*/

		ICM_Option_Tranche( const ARM_Date& TrigerStartDate,
						const ARM_Date& ExpiryDate,
						double strike,
						double InitSpread,
						int optionType,
						ICM_Security* Underlying,
						double rehauss,
						double Freq,
						int ExecType,
						int DiffCDO,
						int IsCMSpread,
						double CMSpreadMatu)
		{
			Init();
			Set(TrigerStartDate, ExpiryDate,strike,InitSpread, optionType, Underlying, rehauss, Freq, ExecType, DiffCDO,IsCMSpread,CMSpreadMatu);
		}

                          
		/*void Set( const ARM_Date& underMaturity,
						const ARM_Date& ExpiryDate,
						double strike,
						int optionType,
						ICM_Security* Underlying,
						double rehauss,
						double Freq,
						int ExecType,
						int DiffCDO);*/
		void Set( const ARM_Date& TrigerStartDate,
						const ARM_Date& ExpiryDate,
						double strike,
						double InitSpread,
						int optionType,
						ICM_Security* Underlying,
						double rehauss,
						double Freq,
						int ExecType,
						int DiffCDO,
						int IsCMSpread,
						double CMSpreadMatu);


		~ICM_Option_Tranche(void){}

		// Setter/Getter
		void SetUnderlying(ICM_Security* Underlying){ itsUnderlying = Underlying;}
		ICM_Security* GetUnderlying(){ return itsUnderlying ;}

		void SetRehauss(double Rehauss){ itsRehauss = Rehauss;}
		double GetRehauss(){ return itsRehauss ;}
       
		double GetTriggerFreq(){ return itsObservationFreq ;}
		void SetTriggerFreq(double TriggerFreq){ itsObservationFreq = TriggerFreq;}
       
		int GetDiffCDO(){ return itsDiffCDO ;}

		void SetDiffCDO(int DiffCDO){ itsDiffCDO = DiffCDO;}

		double GetInitSpread(){return itsInitSpread;}
		void SetInitSpread( double InitSpread) { itsInitSpread = InitSpread;}

		int IsCMSpread() {return itsIsCMSpread;}
		double GetCMSpreadMatu(){return itsCMSpreadMatu;}

		double SetCMSpreadMatu(double matu)
		{
			itsIsCMSpread = 1 ;
			itsCMSpreadMatu = matu;
		}
		void View(char* id = NULL, FILE* ficOut = NULL);
		

};


#endif
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/






