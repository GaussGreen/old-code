/* ===============================================================================

   FILENAME	:	srt_h_und_utils.h
	
   PURPOSE:     Utility Functions to deal with underlying for the outside world
   
   =============================================================================== */

#ifndef SRT_H_UND_UTILS_H
#define SRT_H_UND_UTILS_H


SRT_Boolean srt_f_IsUnderlyingDefined(char *und_name);

SRT_Boolean srt_f_IsUnderlyingInterestRate(char *und_name);

SRT_Boolean srt_f_IsUnderlyingIrOneFactor(char *und_name);

SRT_Boolean srt_f_IsUnderlyingIrTwoFactor(char *und_name);

void srt_f_GetUnderlyingModelName(char *und_name, char *mdl_name);

void srt_f_GetUnderlyingYieldCurveName(char *und_name, char *yc_name);

void srt_f_GetUnderlyingCurrency(char *und_name, char *ccy_name);

#endif