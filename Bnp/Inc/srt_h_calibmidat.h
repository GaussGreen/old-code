/****************************************************************/
/* Calboot with swaptions for Midat, and so on					*/
/****************************************************************/
Err SrtGrfnSwaptionCalboot (long			*plStartDates,
							long			*plEndDates,
							String			szFrequency,
							String			szBasis,
							String			szRefRateCode,
							double 			*dStrike,
							String			szRecPay,
							long			lNumDates,
							String			szModelName,
							/* LGM Parameters */
							double			dTau,
							SRT_Boolean			bUseTwoFactor,
							double			dAlpha,
							double			dGamma,
							double			dRho,
							/* GRFN Parameters for underlying */
							String			*pszGrfnParamNames,
							String			*pszGrfnParamValues,
							int				iNumGrfnParams, 
							String			szYieldCurveName,
							Err				(*pfGetVol)(double dStart, double dEnd, double dStrike, 
														double dForward, double dSpread, double *pdBsVol),
							String			szVolType,

							/* Outputs from this calibration */	
							String          szUndName);