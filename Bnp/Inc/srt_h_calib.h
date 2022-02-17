/* ========================================================================

  FILENAME:  srt_h_calib.h	


  PURPOSE:   The core for any interest rate calibration.
             
  ========================================================================== */

#ifndef SRT_H_CALIB_H
#define SRT_H_CALIB_H


Err srt_f_calib_stoch_vol(
		SrtGrfnParam     *grfnparam,
		SrtMdlType        mdl_type,
		SrtMdlDim		  mdl_dim,
		SwapDP            *sdp,                    
		double            *strike,
		double            *bond_strike,
		StructType 		  *opt_type,
		SrtReceiverType   *rec_pay,
		double            *price,
		double            *vega,
		String            *ref_rate_code,
		int               num_instruments,
		double            *fra_mat,
		double            **corr_matrix,
		long              num_tenors,
		SrtUndPtr         undptr,
		SrtCalibParam     *calib_param,
		double            **sigDatesValues,
		long              numSigs,
		long              numSigCols,
		double            **tauDatesValues,
		long              numTaus,
		long              numTauCols,
		double            *chisq);

Err srt_f_calib_main(
		SrtGrfnParam      *psGrfnParam,
		SrtMdlType        eModelType,
		SrtMdlDim		    eModelDim,
		SwapDP            *psSwapDp,                    
		double            *pdStrike,
		double            *pdBondStrike,
		StructType 		  *peProductType,
		SrtReceiverType   *peRecPay,
		double            *pdMktPrice,
		double            *pdMktVega,
		String           *pszRefRateCode,
		long                lNumInstruments,
		double            *pdFraMaturities,
		double          **ppdCorrelationMatrix,
		long               plNumTenors,
		SrtUndPtr           sUndPtr,
		SrtCalibParam     *psCalibParams,
		double          **ppdSigmaValues,
		long                lNumSigmas,
		long                lNumSigmaCols,
		double          **ppdTauValues,
		long                lNumTaus,
		long                lNumTauCols,
		double			  *pdChiSquare);

Err srt_f_calib_core(
		SrtGrfnParam      *psGrfnParam,
		SrtMdlType        eModelType,
		SrtMdlDim		    eModelDim,
		SwapDP            *psSwapDp,                    
		double            *pdStrike,
		double            *pdBondStrike,
		StructType 		  *peProductType,
		SrtReceiverType   *peRecPay,
		double            *pdMktPrice,
		double            *pdMktVega,
		String           *pszRefRateCode,
		int	                lNumInstruments,
		double            *pdFraMaturities,
		double          **ppdCorrelationMatrix,
		long               plNumTenors,
		SrtUndPtr           sUndPtr,
		SrtCalibParam     *psCalibParams,
		double          **ppdSigmaValues,
		long                lNumSigmas,
		long                lNumSigmaCols,
		double          **ppdTauValues,
		long                lNumTaus,
		long                lNumTauCols,
		double			  *pdChiSquare);


Err srt_f_bootstrap_calibrate(
		SrtUndPtr          undptr,
		SrtGrfnParam       *grfnparam,
		SwapDP             *sdp,
		double             *pdStrikeRate,
		double             *bnd_strike,
		int	               *type,
		SrtReceiverType	   *recpay,
		String             *pszRefRateCode,
		double             *price,
		double			   *ATMprice,
		double             tau,
		double             alpha,
		double             gamma,
		double             rho,
		int                *num_inp);


Err  srt_f_calib_fx_core(char		**pszFXCalibStrings,
						 char		**pszFXCalibValueStrings,
						 int		iNumCalibParamss,
					
						 int		iNumParams,

						 int		iNumLocVols,
						 double		**ppdFXVolCrv,
							  
   					     int		iNumInstrs,
					     long		*plFXOptMDates,
					     double		*pdFXOptMktVols,
					     double		*pdWeights,
					  
					     char		*szFXUndName,				  
					     double		*ChiSquare);




#endif

