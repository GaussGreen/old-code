/*
 * $Log: autocal.h,v $
 * Revision 1.3  2002/07/02 10:25:02  mab
 * Yan Li add-on
 *
 */



#ifndef _AUTOCAL_H
#define _AUTOCAL_H



/*----------------------------------------------------------------------------*/







class ARM_Frm_AutoCalibrator : public ARM_FRM_InterpCalibrator
{
     public :

         ARM_VolCurve* itsSwpBlackVols;  // Market input (sticky delta)
         ARM_VolCurve* itsCapBlackVols;

         ARM_VolCurve* itsSwpVolsForComputations; // Used for valuation 
                                                  // in the tree
                  // (Market vols may have been switched to sticky strike)

         ARM_VolCurve*       itsSwapSmile;
         ARM_VolCurve*       itsCapSmile;

         int                 itsMode;

         ARM_Portfolio*      itsCalibPortfolio;
         ARM_Security*       itsPricedSecurity;

         double              itsMinVol;
         double              itsMaxVol;
         double              itsAccuracy;
         long                itsNbMaxStep;

         double*             itsMatCurve;
         double*             itsValCurve;
         int                 itsCurveSize;

         ARM_Matrix*         itsMatCorr;

         ARM_Vector*         itsSwp2CalVols;
         ARM_Vector*         itsCap2CalVols;
         ARM_Vector*         itsCalVols;


     public :

         ARM_Frm_AutoCalibrator(void)
         {
             Init();
         }

         // 1 vol type constructor : IDX or IRG
         ARM_Frm_AutoCalibrator(int type, ARM_VolCurve* VolCurve,
                                ARM_VolCurve* Smile,
                                ARM_Vector* CorrelatedIndexes= NULL, 
                                ARM_Vector* indexes = NULL,
                                ARM_Matrix* correlations = NULL);

         // 2 vol type constructor
         ARM_Frm_AutoCalibrator(ARM_VolCurve* VolIdxCurve,
                                ARM_VolCurve* IdxSmile,
                                ARM_VolCurve* VolIrgCurve,
                                ARM_VolCurve* IrgSmile,
                                ARM_Vector* CorrelatedIndexes= NULL, 
                                ARM_Vector* indexes = NULL,
                                ARM_Matrix* correlations = NULL);

        ~ARM_Frm_AutoCalibrator(void);


         void SetVolOptimParams(double minVol = 1e-8,
                                double maxVol = 1, 
                                double accuracy = 1e-8,
                                long maxIters = 1e8);

         void Init(void)
         {
             ARM_FRM_InterpCalibrator::Init();

             itsSwpBlackVols = NULL;
             itsCapBlackVols = NULL;
             itsSwapSmile = NULL;
             itsCapSmile = NULL;
             itsMode = AUTO_MODE::NONE;
             itsCalibPortfolio = NULL;

             SetCalibratingModel(NULL);

             itsMatCurve = NULL;
             itsValCurve = NULL;
             itsMatCorr  = NULL;


             itsMinVol    = 1e-8;
             itsMaxVol    = 2.0;
             itsAccuracy  = 1e-8;
             itsNbMaxStep = 100;

             killCalibratingModel = true;

             itsSwp2CalVols = NULL;
             itsCap2CalVols = NULL;
             itsCalVols     = NULL;
         }

         virtual void CreateCalibPortfolio(void);

         virtual void CreateAnaModel(void);

         double BootstrapAnaModel(int mode = 0);

         void GenerateVolMatrix(int vol = VOL_TYPE::IDX);


         bool IsAutoCalibrator()
         {
             return true;
         }

         inline ARM_Security* GetPricedSecurity(void)
         {
             return itsPricedSecurity;
         }

         inline void SetPricedSecurity(ARM_Security *sec)
         {
             itsPricedSecurity = sec;
         }

         inline ARM_Vector* GetCap2CalVols(void) {return itsCap2CalVols;}
         inline ARM_Vector* GetSwp2CalVols(void) {return itsSwp2CalVols;}
         inline ARM_Vector* GetCalVols(void)     {return itsCalVols;}

         inline void SetCap2CalVols(ARM_Vector* in) {itsCap2CalVols = in;}
         inline void SetSwp2CalVols(ARM_Vector* in) {itsSwp2CalVols = in;}
         inline void SetCalVols(ARM_Vector* in)     {itsCalVols = in;}

         inline ARM_VolCurve*    GetSwpBlackVols(void) {return itsSwpBlackVols;}
         inline ARM_VolCurve*    GetCapBlackVols(void) {return itsCapBlackVols;}
         inline ARM_VolCurve*    GetSwapSmile(void)    {return itsSwapSmile;}
         inline ARM_VolCurve*    GetCapSmile(void)     {return itsCapSmile;}

         inline int GetMode(void)                      {return itsMode;}
         inline int GetCurveSize(void)                 {return itsCurveSize;}

         inline ARM_Portfolio* GetCalibPf(void)  {return itsCalibPortfolio;}

         inline void SetSwpBlackVols(ARM_VolCurve* m)    {itsSwpBlackVols = m;}
         inline void SetCapBlackVols(ARM_VolCurve* m)    {itsCapBlackVols = m;}
         inline void SetSwapSmile(ARM_VolCurve* m)       {itsSwapSmile = m;}
         inline void SetCapSmile(ARM_VolCurve* m)        {itsCapSmile = m;}

         inline void Setmode(int m)               {itsMode = m;}
         inline void SetCurveSize(int n)          {itsCurveSize = n;}

         inline void SetCalibPf(ARM_Portfolio* p) {itsCalibPortfolio = p;}


         inline double* GetMatCurve(void)            {return itsMatCurve;}
         inline double* GetValCurve(void)            {return itsValCurve;}

         inline void SetMatCurve(double* c)          {itsMatCurve = c;}
         inline void SetValCurve(double* c)          {itsValCurve = c;}

         inline ARM_Matrix* GetMatCorr(void)          {return itsMatCorr;}
         inline void        SetMatCorr(ARM_Matrix* m) {itsMatCorr = m;}



         void Output(void);

         void EndCalibration(void);

         void PrepareCalibration(void);    

         void SwitchSwopt2StickyStrike(void);

};


#endif

/*----------------------------------------------------------------------------*/
/*---- End Of File ----*/
