
#ifndef _ICM_CREDIT_MANAGER_H
#define _ICM_CREDIT_MANAGER_H


/*********************************************************************************/
/*! \class  ICM_Credit_Manager icm_credit_manager.h "icm_credit_manager.h"
 *  \author 
 *	\version 1.0
 *	\date   Sept 15 2004
 *	\file   icm_credit_manager.h
 *	\brief  Credit Manager */
/***********************************************************************************/

#include "ARMKernel\inst\security.h"

class ICM_Credit_Manager : public ARM_Object
{

	// ----------------------------------------------------------
	// DATA
	// ----------------------------------------------------------

    private:

		bool			itsActivateFlag;

		ARM_Date		itsValDate;		// AsOfDate
		
		ARM_Date		itsMaxCDSDate;	// MaxCDSDate --> will have to take PayLag into account

		// PreComputes the ZC Values for CDS Calibration
		ARM_Vector*		itsZCValuesForCalibration;

		double			itsElapsed_Time;	// to monitor the time

	// ----------------------------------------------------------
	// ----------------------------------------------------------
    
	public:

        ICM_Credit_Manager(void) { Init();}

        void BitwiseCopy(const ARM_Object* src);

        void Copy(const ARM_Object* srcZc);
 
        ARM_Object* Clone(void)
        {
            ICM_Credit_Manager* theClone = new ICM_Credit_Manager();
 
             theClone->Copy(this);
 
            return(theClone);
        }


		~ICM_Credit_Manager(void)
		{
			if (itsZCValuesForCalibration)
				delete itsZCValuesForCalibration;
			itsZCValuesForCalibration = NULL;
		}

        void Init(void);

		void Reset();

	// ----------------------------------------------------------
	// SPECIFIC CLASS IMPLEMENTATION
	// ----------------------------------------------------------

	public :
		
		inline bool IsActivated(void) { return itsActivateFlag;}
		inline void SetActivateFlag(bool ActivateFlag) { itsActivateFlag = ActivateFlag;}

		inline void	SetValDate(const ARM_Date& ValDate) {itsValDate = ValDate; }
		inline ARM_Date GetValDate(void) { return itsValDate;}

		inline void	SetMaxCDSDate(const ARM_Date& ValDate) {itsMaxCDSDate = ValDate; }
		inline ARM_Date GetMaxCDSDate(void) { return itsMaxCDSDate;}

		inline void SetZCValuesForCalibration(ARM_Vector* zc) 
		{
			if (itsZCValuesForCalibration) 
				delete itsZCValuesForCalibration;
			itsZCValuesForCalibration = zc; 
		}

		inline ARM_Vector* GetZCValuesForCalibration(void) { return itsZCValuesForCalibration;}

		inline void	SetElasped_Time(double elapsed_time) {itsElapsed_Time = elapsed_time; }
		inline double GetElasped_Time(void) { return itsElapsed_Time;}

		double	DiscountPrice(int i);

		void Display();

};		



#endif /*---- End of file ----*/
