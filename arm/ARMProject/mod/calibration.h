/*
 $Log: calibration.h,v $
 Revision 1.6  2003/02/11 17:16:29  mab
 Improvements

 Revision 1.5  2002/11/21 16:40:59  mab
 Added const in Constructor definition

 Revision 1.4  2002/11/21 16:37:17  mab
 SetName put in Init

 Revision 1.3  2001/01/23 09:08:30  nicolasm
 Suppression de delete itsModel
 et ajout if (itsZeroCurve) dans le destructeur

 Revision 1.2  2000/10/12 16:26:58  mab
 Ajout fonction SetZeroCurve

 Revision 1.1  1999/11/10 08:54:59  nicolasm
 Initial revision

 */

/*----------------------------------------------------------------------------*/
#ifndef _CALIB_H
#define _CALIB_H



#include "armglob.h"
#include "model.h"




class ARM_Calibration :  public ARM_Object
{
    private :

        ARM_Date itsAsOf;


        ARM_Model* itsModel;

        ARM_ZeroCurve* itsZeroCurve;

       
        int itsCalibrationMode; // K_GLOBAL_CALIB or K_BOOTSTRAP_CALIB

    public :


        ARM_Calibration(void) 
        {
            Init();
        }

        ARM_Calibration(ARM_Model* model)
        {
            Init();

            itsModel = model;
        }

        ARM_Calibration(ARM_Date& asOf, ARM_ZeroCurve* zc);
  

        ARM_Calibration(const ARM_Calibration& calib) : ARM_Object(calib)
        {
            Init();

            BitwiseCopy(&calib);
        }

       ~ARM_Calibration(void) 
        {
            if (itsZeroCurve)
               delete itsZeroCurve;
        }

        ARM_Calibration& operator = (const ARM_Calibration &calib)
        {
            (*this).ARM_Object::operator = (calib);

            BitwiseCopy(&calib);

            return(*this);
        }

        ARM_Model* GetModel(void)
        {
            return(itsModel);
        }

        int GetCalibrationMode(void)
        {
            return(itsCalibrationMode);
        }

        void SetGetCalibrationMode(int mode) 
        {
            itsCalibrationMode = mode;
        }
 
        void SetZeroCurve(ARM_ZeroCurve* zc)
        {
            itsZeroCurve = (ARM_ZeroCurve *) zc->Clone();
        }

        ARM_ZeroCurve* GetZeroCurve(void)
        {
            return(itsZeroCurve);
        }

        ARM_Date& GetAsOfDate(void)
        {
            return(itsAsOf);
        }

        // ARM Stuff

        void Init(void)
        {
            SetName(ARM_CALIBRATION);

            itsModel = NULL;

            itsZeroCurve = NULL;

            itsCalibrationMode = K_GLOBAL_CALIB;
        }

        void BitwiseCopy(const ARM_Object *ocalib)
        {
            ARM_Calibration* calib = (ARM_Calibration *) ocalib;


            if (calib->itsModel)
				itsModel = (ARM_Model *) calib->itsModel->Clone();    

            if (calib->itsZeroCurve)
				itsZeroCurve = (ARM_ZeroCurve *) calib->itsZeroCurve->Clone();

            itsCalibrationMode = calib->itsCalibrationMode;
        }

        void Copy(const ARM_Object* ocalib)
        {
            ARM_Object::Copy(ocalib);

            BitwiseCopy(ocalib);
        }

        ARM_Object* Clone(void)
        {
            ARM_Calibration* theClone = new ARM_Calibration();


            theClone->Copy(this);

            return(theClone);
        }
};



#endif
/*----------------------------------------------------------------------------*/
/*---- End Of File ----*/
