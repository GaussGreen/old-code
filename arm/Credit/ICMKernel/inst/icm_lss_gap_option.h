#ifndef _ICM_LSS_GAP_OPTION_H
#define _ICM_LSS_GAP_OPTION_H

#include "armglob.h"
#include "security.h"
#include "ICMKernel/inst/icm_mez.h"
#include "ICMKernel/pricer/icm_pricer.h"

class ICM_Lss_Gap_Option : public ICM_Security
{
    private:

		ICM_Mez*		its_Underlying_Cdo; 
		ICM_QMatrix<double>*	its_Spread_Trigger;
		ICM_QMatrix<double>*	its_Default_Trigger;
		ICM_QMatrix<double>*	its_MTM_Trigger;

		double			its_MTM_Single_condition;

    public:

        ICM_Lss_Gap_Option(void) {Init();}

		ICM_Lss_Gap_Option(ICM_Mez* cdo,
						   ICM_QMatrix<double>* spread_triggers,
						   ICM_QMatrix<double>* default_triggers,
						   ICM_QMatrix<double>* mtm_triggers,
						   double	single_cond) 
		{
			Init();
			Set(cdo,spread_triggers,default_triggers,mtm_triggers,single_cond);	
		}

		void Set(ICM_Mez* cdo,
						   ICM_QMatrix<double>* spread_triggers,
						   ICM_QMatrix<double>* default_triggers,
						   ICM_QMatrix<double>* mtm_triggers,
						   double	single_cond) 
		{
				if (cdo)
				{its_Underlying_Cdo = (ICM_Mez*) cdo->Clone();}

				if (spread_triggers)
				{its_Spread_Trigger = (ICM_QMatrix<double>*) spread_triggers->Clone();}

				if (default_triggers)
				{its_Default_Trigger = (ICM_QMatrix<double>*) default_triggers->Clone();}

				if (mtm_triggers)
				{its_MTM_Trigger = (ICM_QMatrix<double>*) mtm_triggers->Clone();}

				its_MTM_Single_condition=single_cond;	
		}

		void Init()
		{
		 SetName(ICM_LSS_GAP_OPTION);
		 its_Underlying_Cdo=NULL;
		 its_Spread_Trigger=NULL;
		 its_Default_Trigger=NULL;
		 its_MTM_Trigger=NULL;
		 its_MTM_Single_condition=-1.;	
		}

		~ICM_Lss_Gap_Option(void)
		{if (its_Underlying_Cdo)
			 delete its_Underlying_Cdo;
		 its_Underlying_Cdo=NULL;
		 if (its_Spread_Trigger)
			 delete its_Spread_Trigger;
		 its_Spread_Trigger=NULL;
		 if (its_Default_Trigger)
			 delete its_Default_Trigger;
		 its_Default_Trigger=NULL;
		 if (its_MTM_Trigger)
			 delete its_MTM_Trigger;
		 its_MTM_Trigger=NULL;

		 its_MTM_Single_condition=-1.;		
		};

		void View(char* id = NULL, FILE* ficOut = NULL)
		{
		    FILE* fOut;
			char  fOutName[200];

		    if ( ficOut == NULL )
			{
			ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
	
		   (void) unlink(fOutName);

			fOut = fopen(fOutName, "w"); 
			}
			else
			{
			fOut = ficOut;
			} 

			fprintf(fOut, "\t\t\t ----------------- Spread Triggers ----------------- \n");
			its_Spread_Trigger->View(id, fOut);
			fprintf(fOut, "\t\t\t ----------------- Default Triggers ----------------- \n");
			its_Default_Trigger->View(id, fOut);
			fprintf(fOut, "\t\t\t ----------------- Mtm Triggers ----------------- \n");
			its_MTM_Trigger->View(id, fOut);

			if ( ficOut == NULL )
			{
			fclose(fOut);
			}
		}

		void BitwiseCopy(const ARM_Object* object)
		{
			    ICM_Lss_Gap_Option* option = (ICM_Lss_Gap_Option*) object;
				if (its_Underlying_Cdo)
				{its_Underlying_Cdo = (ICM_Mez*) option->its_Underlying_Cdo->Clone();}

				if (its_Spread_Trigger)
				{its_Spread_Trigger = (ICM_QMatrix<double>*) option->its_Spread_Trigger->Clone();}

				if (its_Default_Trigger)
				{its_Default_Trigger = (ICM_QMatrix<double>*) option->its_Default_Trigger->Clone();}

				if (its_MTM_Trigger)
				{its_MTM_Trigger = (ICM_QMatrix<double>*) option->its_MTM_Trigger->Clone();}

				its_MTM_Single_condition=option->its_MTM_Single_condition;	
		}

		void Copy(const ARM_Object* srcOption)
		{
			   ICM_Security::Copy(srcOption);
			   BitwiseCopy(srcOption);
		}

		ARM_Object* Clone(void)
		{
			ICM_Lss_Gap_Option* theClone = new ICM_Lss_Gap_Option();
			theClone->Copy(this);
			return(theClone);
		}
		
		ICM_Mez* GetUnderlyingCdo() {return	its_Underlying_Cdo;} 
		ICM_QMatrix<double>* GetSpreadTrigger() {return	its_Spread_Trigger;}
		ICM_QMatrix<double>* GetDefaultTrigger() {return its_Default_Trigger;}
		ICM_QMatrix<double>* GetMTMTrigger() {return its_MTM_Trigger;}
		double GetMTMSinglecondition() {return its_MTM_Single_condition;}

};


#endif
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/






