/*
 $Log: calibrationfrm.h,v $
 Revision 1.5  2003/09/08 12:36:26  mab
 Init() Method Added in default constructor+Formatting

 Revision 1.4  2003/09/02 16:06:38  sgasquet
 Ajout constructeur prenant le modele de marche en input

 Revision 1.3  2002/01/17 14:22:52  sgasquet
 Ajout parametres de calage

 Revision 1.2  2002/01/16 14:48:42  sgasquet
 Ajout infos pour calage stickys

 */

/*----------------------------------------------------------------------*/
#ifndef _CALIBFRM_H
#define _CALIBFRM_H


#include "volcurv.h"
#include "zerocurv.h"
#include "calibration.h"
#include "frmana.h"
#include "sigmadates.h"
#include "bsmodel.h"
#include "calibrationhwsv.h"


class ARM_CalibrationFRM :  public   ARM_Calibration
{
private :

    ARM_Security* itsSecurity;
    ARM_Model* itsMktModelToFit;
    ARM_VolCurve* itsVolCurve;

    double itsIndicPrix; // indicateur de qualite de calib sur prix
    double itsIndicVol;  // indicateur de qualite de calib sur Vol

    int itsMfine ;
    int itsShape_type ;
    double itsDecay ;
    double itsSlope ;
    double itsAsymptote;
    int itsNbFactor;

    ARM_Vector* itsCorrelatedIndexes;
    ARM_Vector* itsIndexes;
    ARM_Matrix* itsCorrelations;
    ARM_Portfolio *itsPf;
    double itsprecisionVol;
    double itsminVols;
    double itsmaxVols;
    long itsnbMaxIter;


public :
                            
    ARM_CalibrationFRM(ARM_Date &asOf, ARM_ZeroCurve *zc,  
                       ARM_VolCurve *VolCurve,
                       ARM_Security *sec,
                       int mfine = 3, int shape_type = K_ROW , 
                       double decay=0.0, double slope=0.0, 
                       double asymptote=.0, int NbFactor = 1, 
                       double precisionVol=0.0001, double minVols=0.001,
                       double maxVols=0.5, long nbMaxIter=500,
                       ARM_Vector* CorrelatedIndexes = NULL,  
                       ARM_Vector* Indexes = NULL, 
                       ARM_Matrix* Correlations = NULL,
                       ARM_Portfolio *itsPf = NULL);

    ARM_CalibrationFRM(ARM_Date &asOf, ARM_Model *itsMktModelToFit,
                       ARM_Security *sec,
                       int mfine = 3, int shape_type = K_ROW , 
                       double decay=0.0, double slope=0.0, 
                       double asymptote=.0, int NbFactor = 1, 
                       double precisionVol=0.0001, double minVols=0.001,
                       double maxVols=0.5, long nbMaxIter=500,
                       ARM_Vector* CorrelatedIndexes = NULL,  
                       ARM_Vector* Indexes = NULL, 
                       ARM_Matrix* Correlations = NULL,
                       ARM_Portfolio* itsPf = NULL);


    ARM_CalibrationFRM(void) 
    {
        Init();
    }


    ~ARM_CalibrationFRM(void);


    ARM_CalibrationFRM(const ARM_CalibrationFRM& calib) : ARM_Calibration(calib)
    {
        Init();

        BitwiseCopy(&calib);
    }


    ARM_CalibrationFRM &operator = (const ARM_CalibrationFRM &calib)
    {
        (*this).ARM_Calibration::operator = (calib);

        BitwiseCopy(&calib);

        return (*this);
    }

    
    ARM_FrmAna* Calibrate(char* msg,
                          ARM_Calibration_TYPE caltype,
                          ARM_CalibrationHWSV_OUTPUT& Output);



    // ****** Services de Copie *****
    void Init(void)
    {
        SetName(ARM_FRMCALIBRATION);

        itsSecurity = NULL;
        itsMktModelToFit = NULL;
        itsVolCurve = NULL;
        itsIndicPrix = 0.0 ;
        itsIndicVol = 0.0 ;
        itsMfine = 3;
        itsShape_type = K_ROW ;
        itsDecay = 0.0;
        itsSlope = 0.0;
        itsAsymptote = 0.0;
        itsNbFactor = 1;
        itsCorrelatedIndexes = NULL;
        itsIndexes = NULL;
        itsCorrelations = NULL;
        itsPf = NULL;
    }


    void BitwiseCopy(const ARM_Object *ocalib)
    {
        ARM_CalibrationFRM* calib = (ARM_CalibrationFRM *) ocalib;

        itsSecurity = calib->itsSecurity;
        itsMktModelToFit = calib->itsMktModelToFit;
        itsVolCurve = calib->itsVolCurve;
        itsIndicPrix = calib->itsIndicPrix;
        itsIndicVol = calib->itsIndicVol;
        itsMfine = calib->itsMfine;
        itsShape_type = calib->itsShape_type ;
        itsDecay = calib->itsDecay;
        itsSlope = calib->itsSlope;
        itsAsymptote = calib->itsAsymptote;
        itsNbFactor = calib->itsNbFactor;
        
        if (calib->itsCorrelatedIndexes)
           itsCorrelatedIndexes = (ARM_Vector *) calib->itsCorrelatedIndexes->Clone();

        if (calib->itsIndexes)
           itsIndexes =  (ARM_Vector *) calib->itsIndexes->Clone();

        if (calib->itsCorrelations)
          itsCorrelations =  (ARM_Matrix *) calib->itsCorrelations->Clone();

        itsPf = calib->itsPf;
    }


    void Copy(const ARM_Object* ocalib)
    {
        ARM_Calibration::Copy(ocalib);

        BitwiseCopy(ocalib);
    }


    ARM_Object* Clone(void)
    {
        ARM_CalibrationFRM* theClone = new ARM_CalibrationFRM();

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

};



#endif
/*----------------------------------------------------------------------*/
/*---- End Of File ----*/
