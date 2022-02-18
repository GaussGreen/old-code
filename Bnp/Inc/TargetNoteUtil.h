#ifndef __COBRA_WRAP__
#ifndef _TARGET_NOTE_UTIL_
#define _TARGET_NOTE_UTIL_

#include "TargetNoteProdStruct.h"

init_TARN_Struct(TARN_Struct* tarn);
copy_TARN_Struct(TARN_Struct* in_tarn, TARN_Struct* out_tarn);
free_TARN_Struct(TARN_Struct* tarn);

/*
init_TARN_All_Struct(
                                                TARN_Market_Struct* tarn_market_struct,
                                                TARN_Deal_Struct* tarn_deal_struct,
                                                TARN_Model_Struct* tarn_model_struct,
                                                TARN_Calibration_Struct* tarn_calibration_struct,
                                                TARN_Pricing_Struct* tarn_pricing_struct,
                                                TARN_Output_Struct* tarn_output_struct
                                        );

free_TARN_All_Struct(
                                                TARN_Market_Struct* tarn_market_struct,
                                                TARN_Deal_Struct* tarn_deal_struct,
                                                TARN_Model_Struct* tarn_model_struct,
                                                TARN_Calibration_Struct* tarn_calibration_struct,
                                                TARN_Pricing_Struct* tarn_pricing_struct,
                                                TARN_Output_Struct* tarn_output_struct
                                        );
*/

init_TARN_Market_Struct(TARN_Market_Struct* tarn_market_struct);
copy_TARN_Market_Struct(
    TARN_Market_Struct* in_market_struct, TARN_Market_Struct* out_market_struct);
free_TARN_Market_Struct(TARN_Market_Struct* tarn_market_struct);

init_TARN_Deal_Struct(TARN_Deal_Struct* tarn_deal_struct);
copy_TARN_Deal_Struct(TARN_Deal_Struct* in_deal_struct, TARN_Deal_Struct* out_deal_struct);
free_TARN_Deal_Struct(TARN_Deal_Struct* tarn_deal_struct);

init_TARN_Model_Struct(TARN_Model_Struct* tarn_model_struct);
copy_TARN_Model_Struct(TARN_Model_Struct* in_model_struct, TARN_Model_Struct* out_model_struct);
free_TARN_Model_Struct(TARN_Model_Struct* tarn_model_struct);

init_TARN_Calibration_Struct(TARN_Calibration_Struct* tarn_calibration_struct);
copy_TARN_Calibration_Struct(
    TARN_Calibration_Struct* in_calibration_struct,
    TARN_Calibration_Struct* out_calibration_struct);
free_TARN_Calibration_Struct(TARN_Calibration_Struct* tarn_calibration_struct);

init_TARN_Pricing_Struct(TARN_Pricing_Struct* tarn_pricing_struct);
copy_TARN_Pricing_Struct(
    TARN_Pricing_Struct* in_pricing_struct, TARN_Pricing_Struct* out_pricing_struct);
free_TARN_Pricing_Struct(TARN_Pricing_Struct* tarn_pricing_struct);

init_TARN_Output_Struct(TARN_Output_Struct* tarn_output_struct);
copy_TARN_Output_Struct(
    TARN_Output_Struct* in_output_struct, TARN_Output_Struct* out_output_struct);
free_TARN_Output_Struct(TARN_Output_Struct* tarn_output_struct);

char* TARN_calcHistory(TARN_Struct* tarn, TARN_AUX* aux);

#endif /* _TARGET_NOTE_UTIL_ */
#endif