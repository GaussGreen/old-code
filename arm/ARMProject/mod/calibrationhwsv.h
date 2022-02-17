/*
 $Log: calibrationhwsv.h,v $
 Revision 1.16  2003/09/08 11:57:01  mab
 Added Init() Method in default constructor

 Revision 1.15  2003/09/02 16:28:14  sgasquet
  Ajout constructeur prenant le modele de marche en input

 Revision 1.14  2003/07/04 07:23:08  ebenhamou
 move the destructor to the cpp file since we only know the declaration of the ARM_Portfolio and not its definition

 Revision 1.12  2002/11/29 10:21:28  sgasquet
 Correction oubli ";"

 Revision 1.11  2002/11/27 16:30:39  mab
 Ajout du SetName dans la methode Init

 Revision 1.10  2002/11/21 12:28:01  mab
 Init and const added in copy constructor

 Revision 1.9  2002/01/16 14:48:20  sgasquet
  Ajout infos pour calage stickys

 Revision 1.7  2001/11/20 08:58:30  vberger
 correction de l'ambiguitë sur les constructeurs

 Revision 1.6  2001/10/24 09:03:17  vberger
 crëation d'une struture d'input de calage
 modification de la structure d'output
 ajout de menbres tels que volmin, volmax, amin, amax

 Revision 1.5  2001/06/25 11:54:10  vberger
 modif de calib pour la calibration des MCPRCS

 Revision 1.4  2001/05/07 15:54:56  vberger
 ajout de void Calibrate(char * msg,ARM_Calibration_PORTTYPE porttype,ARM_CalibrationHWSV_OUTPUT& Output)

 Revision 1.3  2000/10/03 17:19:50  nicolasm
 Ajout de msg comme parametre a Calibrate

 Revision 1.2  2000/08/03 12:22:57  nicolasm
 Ajout indicateurs

 Revision 1.1  1999/11/10 08:55:21  nicolasm
 Initial revision

 */

/*-------------------------------------------------------------------*/
#ifndef _CALIBHWSV_H
#define _CALIBHWSV_H



#include "volcurv.h"
#include "zerocurv.h"
#include "calibration.h"
#include "x1fhwsigvar.h"
#include "sigmadates.h"
#include "bsmodel.h"


#define MAXNBHWSIGMAPLOT 100
#define MAXNBSECURITYINPF 200




typedef struct _ARM_CalibrationHWSV_INPUT
{
    // HW Input parameters
    ARM_DFHWSigVarTree_CalibrationAction itsCalibrationAction;
    ARM_Calibration_PORTTYPE itsPfType;
    double itsMeanRevMin;
    double itsMeanRevMax;
    double itsVolMin;
    double itsVolMax;

} 
ARM_CalibrationHWSV_INPUT;


typedef struct _ARM_CalibrationHWSV_OUTPUT
{
    // HW Output parameters

    double meanrev;
    int HWVolSize ;
    ARM_Date Asof ;
    ARM_Calibration_CALAGEQUALITY HWVectorStability;
    ARM_Date datesched[MAXNBHWSIGMAPLOT];
    double volsched[MAXNBHWSIGMAPLOT];

    // PF Data : description of the security in pf 

    int PFSize ;
    ARM_Date SecurityExeDate[MAXNBSECURITYINPF];
    ARM_Date SecurityStartDate[MAXNBSECURITYINPF];
    ARM_Date SecurityEndDate[MAXNBSECURITYINPF];
    double SecurityStrike[MAXNBSECURITYINPF];

    // fitting quality in price and BS vols

    double InputVol[MAXNBSECURITYINPF];
    double OutputVol[MAXNBSECURITYINPF];
    double ErrorVol[MAXNBSECURITYINPF];
    double InputPrice[MAXNBSECURITYINPF];
    double OutputPrice[MAXNBSECURITYINPF];
    double ErrorPrice[MAXNBSECURITYINPF];
    double calagequality;

} 
ARM_CalibrationHWSV_OUTPUT;





class ARM_CalibrationHWSV :  public   ARM_Calibration
{
    private :
   
         ARM_VolCurve* itsVolCurve;

         ARM_Model* itsMktModelToFit;

         ARM_Security* itsSecurity;

         double itsamin;

         double itsamax;

         double itsvolmin;

         double itsvolmax;

         double itsIndicPrix; // indicateur de qualite de calib sur prix
        
         double itsIndicVol;  // indicateur de qualite de calib sur Vol

         ARM_Vector* itsSigmaDates;

         ARM_Portfolio* itsPf;

         static ARM_ConfigHWSVCalibration* itsConfigCalib;

    public :
    
        ARM_CalibrationHWSV(ARM_Date& asOf, ARM_Model* model, ARM_Security* sec,
                            double amin, double amax, 
                            double volmin, double volmax, 
                            ARM_Vector* dates = NULL,
                            ARM_Portfolio* pf = NULL);

        ARM_CalibrationHWSV(ARM_Date& asOf, ARM_ZeroCurve* zc, 
                            ARM_VolCurve* vol, ARM_Security* sec,
                            double amin = 0.0001, double amax = 0.5, 
                            ARM_Vector* dates = NULL, 
                            ARM_Portfolio* pf = NULL);

        ARM_CalibrationHWSV(ARM_Date& asOf, ARM_ZeroCurve* zc, 
                            ARM_VolCurve* vol, ARM_Security* sec,
                            double amin, double amax, 
                            double volmin, double volmax, 
                            ARM_Vector* dates = NULL,
                            ARM_Portfolio* pf = NULL,
                            ARM_VolCurve* RhoSwopt = NULL,
                            ARM_VolCurve* NuSwopt = NULL,
                            ARM_VolCurve* BetaSwopt = NULL,
							int isAlphaOrSigma = 1);

        
        ARM_CalibrationHWSV(void) 
        {
            Init();
        }

       ~ARM_CalibrationHWSV(void);

        ARM_CalibrationHWSV(const ARM_CalibrationHWSV& calib) 
                           : ARM_Calibration(calib)
        {
            Init();

            BitwiseCopy(&calib);
        }


        ARM_CalibrationHWSV &operator = (const ARM_CalibrationHWSV &calib)
        {
            (*this).ARM_Calibration::operator = (calib);

            BitwiseCopy(&calib);

            return(*this);
        }

        void Init(void)
        {
            SetName(ARM_HWSVCALIBRATION);

            itsVolCurve = NULL;

            itsSecurity = NULL;

            itsMktModelToFit = NULL;

            itsamin = 0.001;

            itsamax = 0.5;

            itsvolmin = 0.001;

            itsvolmax = 10.0;

            itsSigmaDates = NULL;

            itsPf = NULL;

            SetName(ARM_HWSVCALIBRATION);
        }
    
        ARM_GYCSigVarModel* Calibrate(char* msg);

        void Calibrate(char* msg,
                       ARM_Calibration_PORTTYPE porttype, 
                       ARM_CalibrationHWSV_OUTPUT& Output);
    

        ARM_GYCSigVarModel* Calibrate(char* msg, ARM_Calibration_TYPE porttype,
                                      ARM_CalibrationHWSV_OUTPUT& Output);

  
        // ****** Services *****

        void BitwiseCopy(const ARM_Object* ocalib)
        {
            ARM_CalibrationHWSV* calib = (ARM_CalibrationHWSV *) ocalib;


            itsVolCurve = calib->itsVolCurve;

            itsSecurity = calib->itsSecurity;

            itsMktModelToFit = calib->itsMktModelToFit;

            itsamin = calib->itsamin;

            itsamax = calib->itsamax;

            itsvolmin = calib->itsvolmin;

            itsvolmax = calib->itsvolmax;

            itsSigmaDates = (ARM_Vector *) calib->itsSigmaDates->Clone();

            itsPf = calib->itsPf;
        }


        void Copy(const ARM_Object *ocalib)
        {
            ARM_Calibration::Copy(ocalib);

            BitwiseCopy(ocalib);
        }


        ARM_Object *Clone(void)
        {
            ARM_CalibrationHWSV *theClone = new ARM_CalibrationHWSV();

            theClone->Copy(this);

            return theClone;
        }

        double GetIndicPrix(void)
        {
            return itsIndicPrix;
        }

        double GetIndicVol(void)
        {
            return itsIndicVol;
        }

 static void SetConfigCalib(ARM_ConfigHWSVCalibration* config)
        {
            ARM_CalibrationHWSV::itsConfigCalib = config;
        }

 static ARM_ConfigHWSVCalibration *GetConfigCalib(void)
        {
            return itsConfigCalib;
        }
};



#endif
/*---------------------------------------------------------------*/
/*---- End Of File ----*/
