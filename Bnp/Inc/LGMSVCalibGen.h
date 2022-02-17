#ifndef LGMSVCALIBGENH
#define LGMSVCALIBGENH

#include "DiagCalibDLM.h"
#include "DiagCalibGen.h"
#include "LGMSVClosedformApprox.h"
#include "cpdcalib.h"
#include "lgmsvpde.h"
#include "srt_h_all.h"

Err GetTargetVol_LGMSV(void *Inst, void *GlobalConst, void *Model,

                       CALIBGEN_PARAMS CalibConsts, double *target);

Err GetFirstGuessVol_LGMSV(void *Model, void *GlobalParam, int vol_index,
                           double target, double *vol1);

Err GetSecondGuessVol_LGMSV(void *Model, void *GlobalParam, int vol_index,
                            double vol1, double price1, double target,
                            double *vol2);

Err BumpVol_LGMSV(void *Model, void *GlobalParam, int vol_index, double vol);

Err SetVol_LGMSV(void *Model, CALIBGEN_PARAMS CalibConsts, void *GlobalParam,
                 int vol_index, double vol);

Err GetLimitAndLastVol_LGMSV(void *Model, CALIBGEN_PARAMS CalibParams,
                             void *GlobalParam,

                             int vol_index, double *last_vol,
                             double *limit_down, double *limit_up);

Err ExtrapolVol_LGMSV(void *Model, void *GlobalParam, int last_vol_index);

Err UpdateParamsAfterVol_LGMSV(void *Inst_, void *InstParam, void *GlobalParam,
                               void *Model, CALIBGEN_PARAMS CalibConsts);

Err PriceInstVol_LGMSV(void *Inst_, void *InstParam, void *GlobalParam,
                       void *Model, double *InstPrice);

Err PriceInstVolSmile_LGMSV(void *Inst_, void *InstParam, void *GlobalParam,
                            void *Model, double *InstPrice);

Err GetTargetRho_LGMSV(void *Inst_, void *GlobalConst, void *Model,

                       CALIBGEN_PARAMS CalibConsts, double *target);

Err GetFirstGuessRho_LGMSV(void *Model, void *GlobalConst, int index_smile,
                           double target, double *param1);

Err GetSecondGuessRho_LGMSV(void *Model, void *GlobalConst, int index_smile,
                            double rho1, double price1, double target,
                            double *rho2);

Err GetLimitAndLastRho_LGMSV(void *Model, CALIBGEN_PARAMS CalibConsts,
                             void *GlobalConst, int index_smile,
                             double *last_rho, double *limit_down,
                             double *limit_up);

Err BumpRho_LGMSV(void *Model, void *GlobalParam, int smile_param, double rho);

Err SetRho_LGMSV(void *Model, CALIBGEN_PARAMS CalibConsts, void *GlobalParam,
                 int smile_param, double rho);

Err ExtrapolRho_LGMSV(void *Model, void *GlobalParam, int last_smile_index);

Err UpdateParamsAfterRho_LGMSV(void *Inst_, void *InstParam, void *GlobalParam,
                               void *Model, CALIBGEN_PARAMS CalibConsts);

Err PriceInstRho_LGMSV(void *Inst_, void *InstParam, void *GlobalParam,
                       void *Model, double *InstPrice);

Err GetTargetAlpha_LGMSV(void *Inst_, void *GlobalConst, void *Model,

                         CALIBGEN_PARAMS CalibConsts, double *target);

Err GetFirstGuessAlpha_LGMSV(void *Model, void *GlobalConst, int index_smile,
                             double target, double *param1);

Err GetSecondGuessAlpha_LGMSV(void *Model, void *GlobalConst, int index_smile,
                              double alpha1, double price1, double target,
                              double *alpha2);

Err GetLimitAndLastAlpha_LGMSV(void *Model, CALIBGEN_PARAMS CalibConsts,
                               void *GlobalConst, int index_smile,
                               double *last_alpha, double *limit_down,
                               double *limit_up);

Err UpdateParamsAfterAlpha_LGMSV(void *Inst_, void *InstParam,
                                 void *GlobalParam, void *Model,
                                 CALIBGEN_PARAMS CalibConsts);

Err PriceInstAlpha_LGMSV(void *Inst_, void *InstParam, void *GlobalParam,
                         void *Model, double *InstPrice);

Err BumpAlpha_LGMSV(void *Model, void *GlobalParam, int smile_param,
                    double alpha);

Err SetAlpha_LGMSV(void *Model, CALIBGEN_PARAMS CalibConsts, void *GlobalParam,
                   int smile_param, double alpha);

Err ExtrapolAlpha_LGMSV(void *Model, void *GlobalParam, int last_smile_index);

Err GetTargetLambda_LGMSV(void *Inst_, void *GlobalConst, void *Model,

                          CALIBGEN_PARAMS CalibConsts, double *target);

Err GetFirstGuessLambda_LGMSV(void *Model, void *GlobalConst, int index_lambda,
                              double target, double *param1);

Err GetSecondGuessLambda_LGMSV(void *Model, void *GlobalConst, int index_lambda,
                               double lambda1, double price1, double target,
                               double *lambda2);

Err GetLimitAndLastLambda_LGMSV(void *Model, CALIBGEN_PARAMS CalibConsts,
                                void *GlobalConst, int index_lambda,
                                double *last_lambda, double *limit_down,
                                double *limit_up);

Err BumpLambda_LGMSV(void *Model, void *GlobalParam, int index_lambda,
                     double lambda);

Err SetLambda_LGMSV(void *Model, CALIBGEN_PARAMS CalibConsts, void *GlobalParam,
                    int index_lambda, double lambda);

Err ExtrapolLambda_LGMSV(void *Model, void *GlobalParam, int last_smile_index);

Err UpdateParamsAfterLambda_LGMSV(void *Inst_, void *InstParam,
                                  void *GlobalParam, void *Model,
                                  CALIBGEN_PARAMS CalibConsts);

Err PriceInstLambda_LGMSV(void *Inst_, void *InstParam, void *GlobalParam,
                          void *Model, double *InstPrice);

Err GetTargetRho2_LGMSV(void *Inst_, void *GlobalConst, void *Model,

                        CALIBGEN_PARAMS CalibConsts, double *target);

Err GetFirstGuessRho2_LGMSV(void *Model, void *GlobalConst, int index_rho2,
                            double target, double *param1);

Err GetSecondGuessRho2_LGMSV(void *Model, void *GlobalConst, int index_rho2,
                             double rho21, double price1, double target,
                             double *rho22);

Err GetLimitAndLastRho2_LGMSV(void *Model, CALIBGEN_PARAMS CalibConsts,
                              void *GlobalConst, int index_rho2,
                              double *last_rho2, double *limit_down,
                              double *limit_up);

Err BumpRho2_LGMSV(void *Model, void *GlobalParam, int index_rho2, double rho2);

Err SetRho2_LGMSV(void *Model, CALIBGEN_PARAMS CalibConsts, void *GlobalParam,
                  int index_rho2, double rho2);

Err ExtrapolRho2_LGMSV(void *Model, void *GlobalParam, int last_smile_index);

Err UpdateParamsAfterRho2_LGMSV(void *Inst_, void *InstParam, void *GlobalParam,
                                void *Model, CALIBGEN_PARAMS CalibConsts);

Err PriceInstRho2_LGMSV(void *Inst_, void *InstParam, void *GlobalParam,
                        void *Model, double *InstPrice);

Err GetTargetFlatAlpha_LGMSV(void *Inst_, void *GlobalConst, void *Model,

                             CALIBGEN_PARAMS CalibConsts, double *target);

Err GetFirstGuessFlatAlpha_LGMSV(void *Model, void *GlobalConst,
                                 int index_flat_alpha, double target,
                                 double *param1);

Err GetSecondGuessFlatAlpha_LGMSV(void *Model, void *GlobalConst,
                                  int index_flat_alpha, double flat_alpha1,
                                  double price1, double target,
                                  double *flat_alpha2);

Err GetLimitAndLastFlatAlpha_LGMSV(void *Model, CALIBGEN_PARAMS CalibConsts,
                                   void *GlobalConst, int index_flat_alpha,
                                   double *last_flat_alpha, double *limit_down,
                                   double *limit_up);

Err BumpFlatAlpha_LGMSV(void *Model, void *GlobalParam, int index_flat_alpha,
                        double flat_alpha);

Err SetFlatAlpha_LGMSV(void *Model, CALIBGEN_PARAMS CalibConsts,
                       void *GlobalParam, int index_flat_alpha,
                       double flat_alpha);

Err ExtrapolFlatAlpha_LGMSV(void *Model, void *GlobalParam,
                            int last_smile_index);

Err UpdateParamsAfterFlatAlpha_LGMSV(void *Inst_, void *InstParam,
                                     void *GlobalParam, void *Model,
                                     CALIBGEN_PARAMS CalibConsts);

Err PriceInstFlatAlpha_LGMSV(void *Inst_, void *InstParam, void *GlobalParam,
                             void *Model, double *InstPrice);

Err GetTargetFlatRho_LGMSV(void *Inst_, void *GlobalConst, void *Model,

                           CALIBGEN_PARAMS CalibConsts, double *target);

Err GetFirstGuessFlatRho_LGMSV(void *Model, void *GlobalConst,
                               int index_flat_rho, double target,
                               double *param1);

Err GetSecondGuessFlatRho_LGMSV(void *Model, void *GlobalConst,
                                int index_flat_rho, double flat_rho1,
                                double price1, double target,
                                double *flat_rho2);

Err GetLimitAndLastFlatRho_LGMSV(void *Model, CALIBGEN_PARAMS CalibConsts,
                                 void *GlobalConst, int index_flat_rho,
                                 double *last_flat_rho, double *limit_down,
                                 double *limit_up);

Err BumpFlatRho_LGMSV(void *Model, void *GlobalParam, int index_flat_rho,
                      double flat_rho);

Err SetFlatRho_LGMSV(void *Model, CALIBGEN_PARAMS CalibConsts,
                     void *GlobalParam, int index_flat_rho, double flat_rho);

Err ExtrapolFlatRho_LGMSV(void *Model, void *GlobalParam, int last_smile_index);

Err UpdateParamsAfterFlatRho_LGMSV(void *Inst_, void *InstParam,
                                   void *GlobalParam, void *Model,
                                   CALIBGEN_PARAMS CalibConsts);

Err PriceInstFlatRho_LGMSV(void *Inst_, void *InstParam, void *GlobalParam,
                           void *Model, double *InstPrice);
#endif