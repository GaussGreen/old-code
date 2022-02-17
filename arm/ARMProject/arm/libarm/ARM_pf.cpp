#include "ARM_interglob.h"
#include "ARM_pf.h"






long ARM_PF(const VECTOR<long>& insts,
            const VECTOR<double>& coeffs,
            const VECTOR<double>& marketPrices,
			const VECTOR<double>& precisions,
            ARM_result& result, long objId)
{
    long idx = 0;
    long ARM_REQUEST_ID;
    long ARM_REQUEST_NBPAR = 5;
    CCString stringObjectId;

    

    try
    {
        /*--- parameters checking ---*/
        if(!((insts.size () == coeffs.size ()) && (coeffs.size () == marketPrices.size ())))
        {
            result.setMsg ("ARM_ERR: instruments, coefficients and prices array must have same size");
            return ARM_KO;
        }

		if ( (precisions.size () > 0) && (precisions.size () != insts.size ()) )
        {
            result.setMsg ("ARM_ERR: instruments, coefficients, precisions and prices array must have same size");
            return ARM_KO;
        }

        ARM_CORBA_init ();

        if (CORBA_OBJECT_INTERFACE)
        {
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

            stringObjectId = GetLastCurCellEnvValue ();

            if ( objId != -1 )
            {
                ARM_REQUEST_ID = RPC_SETPF;
                reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
                REQ_SetLong (reqIn, objId, 0);
                idx = 1;
            }
            else
            {
                ARM_REQUEST_ID = RPC_CREATEPF;
                reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
                idx = 0;
            }

            REQ_SetLong (reqIn, insts.size (), idx);
            REQ_SetLongVector (reqIn, insts, ++idx);
            REQ_SetDoubleVector (reqIn, coeffs, ++idx);
            REQ_SetDoubleVector (reqIn, marketPrices, ++idx);
            REQ_SetDoubleVector (reqIn, precisions, ++idx);
            
            /*--- CORBA Server Call ---*/

            CORBA_OBJECT_INTERFACE->Send (*reqIn, reqOut);

            result.set (reqOut);

            /*--- requests freeing ---*/
            REQ_Delete (reqOut);
            REQ_Delete (reqIn);

            ARM_RESULT();
        }
    }

    catch(const CORBA::Exception &e)
    {
        CORBA_ERR();
    }

    return ARM_KO;
}



long ARM_MKTPRICE (long idPF, ARM_result& result)
{
    long ARM_REQUEST_ID;
    long ARM_REQUEST_NBPAR = 1;
    long idx = 0;

    

    try
    {
        ARM_CORBA_init ();

        if (CORBA_OBJECT_INTERFACE)
        {
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

            ARM_REQUEST_ID = RPC_MKTPRICE;

            reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
            
            REQ_SetLong (reqIn, idPF, idx);
            
            /*--- CORBA Server Call ---*/

               CORBA_OBJECT_INTERFACE->Send (*reqIn, reqOut);

            result.set (reqOut);
        
            /*--- requests freeing ---*/
            REQ_Delete (reqOut);
            REQ_Delete (reqIn);
            
              ARM_RESULT();
        }
    }
    
    catch(const CORBA::Exception &e)
    {
        CORBA_ERR();
    }

    return ARM_KO;
}



long ARM_PFSPLESTIMATE (long pfId,
                        long splineId,
                        double settlement,
                        ARM_result& result, long objId)
{
    long idx = 0;
    long ARM_REQUEST_ID;
    long ARM_REQUEST_NBPAR = 3;
    CCString stringObjectId;
            
    try
    {
        ARM_CORBA_init ();

        if(CORBA_OBJECT_INTERFACE)
        {
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

            stringObjectId = GetLastCurCellEnvValue ();

            if (objId != -1)
            {
                ARM_REQUEST_ID = RPC_SETPFSPLESTIMATED;
                reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
                REQ_SetLong (reqIn, objId, 0);
                idx = 1;
            }
            else
            {
                ARM_REQUEST_ID = RPC_PFSPLESTIMATE;
                reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
                idx = 0;
            }

            REQ_SetLong (reqIn, pfId, idx);
            REQ_SetLong (reqIn, splineId, ++idx);
            REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(settlement), ++idx);
            
            /*--- CORBA Server Call ---*/

               CORBA_OBJECT_INTERFACE->Send (*reqIn, reqOut);

            result.set (reqOut);

            /*--- requests freeing ---*/
            REQ_Delete (reqOut);
            REQ_Delete (reqIn);

            ARM_RESULT();
        }
    }

    catch(const CORBA::Exception &e)
    {
        CORBA_ERR();
    }

    return ARM_KO;
}



long ARM_PFVSKESTIMATE (long pfId,
                        double settlement,
                        double a,
                        ARM_result& result, long objId)
{
    long idx = 0;
    long ARM_REQUEST_ID;
    long ARM_REQUEST_NBPAR = 3;
    CCString stringObjectId;

    

    try
    {
        ARM_CORBA_init ();

        if(CORBA_OBJECT_INTERFACE)
        {
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

            stringObjectId = GetLastCurCellEnvValue ();

            if ( objId != -1 )
            {
                ARM_REQUEST_ID = RPC_SETPFVSKESTIMATED;
                reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
                REQ_SetLong (reqIn, objId, 0);
                idx = 1;
            }
            else
            {
                ARM_REQUEST_ID = RPC_PFVSKESTIMATE;
                reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
                idx = 0;
            }

            REQ_SetLong (reqIn, pfId, idx);
            REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(settlement), ++idx);
            REQ_SetDouble (reqIn, a, ++idx);
            
            /*--- CORBA Server Call ---*/

               CORBA_OBJECT_INTERFACE->Send (*reqIn, reqOut);

            result.set (reqOut);

            /*--- requests freeing ---*/
            REQ_Delete (reqOut);
            REQ_Delete (reqIn);

            ARM_RESULT();
        }
    }

    catch(const CORBA::Exception &e)
    {
        CORBA_ERR();
    }

    return ARM_KO;
}



long ARM_PFGYCFIT (long pfId,
                   double curveId,
                   double settlement,
                   ARM_result& result, long objId)
{
    long idx = 0;
    long ARM_REQUEST_ID;
    long ARM_REQUEST_NBPAR = 3;
    CCString stringObjectId;
            
    try
    {
        ARM_CORBA_init ();

        if (CORBA_OBJECT_INTERFACE)
        {
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

            stringObjectId = GetLastCurCellEnvValue ();

            if ( objId != -1 )
            {
                ARM_REQUEST_ID = RPC_SETPFHWXESTIMATED;
                reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
                REQ_SetLong (reqIn, objId, 0);
                idx = 1;
            }
            else
            {
                ARM_REQUEST_ID = RPC_PFHWXESTIMATE;
                reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
                idx = 0;
            }

            REQ_SetLong (reqIn, pfId, idx);
            REQ_SetLong (reqIn, curveId, ++idx);
            REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(settlement), ++idx);
            
            /*--- CORBA Server Call ---*/

               CORBA_OBJECT_INTERFACE->Send (*reqIn, reqOut);

            result.set (reqOut);

            /*--- requests freeing ---*/
            REQ_Delete (reqOut);
            REQ_Delete (reqIn);

            ARM_RESULT();
        }
    }

    catch(const CORBA::Exception &e)
    {
        CORBA_ERR();
    }

    return ARM_KO;
}



long ARM_PFTHEOPRICE (long pfId,
                      long modelId,
                      double settlement,
                      ARM_result& result)
{
    long ARM_REQUEST_ID;
    long ARM_REQUEST_NBPAR = 3;
    long idx = 0;
        
    try
    {
        ARM_CORBA_init ();

        if (CORBA_OBJECT_INTERFACE)
        {
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

            ARM_REQUEST_ID = RPC_PFTHEOPRICE;

            reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
            
            REQ_SetLong (reqIn, pfId, idx);
            REQ_SetLong (reqIn, modelId, ++idx);
            REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(settlement), ++idx);
            
            /*--- CORBA Server Call ---*/

               CORBA_OBJECT_INTERFACE->Send (*reqIn, reqOut);

            result.set (reqOut);
        
            /*--- requests freeing ---*/
            REQ_Delete (reqOut);
            REQ_Delete (reqIn);
            
              ARM_RESULT();
        }
    }
    
    catch(const CORBA::Exception &e)
    {
        CORBA_ERR();
    }

    return ARM_KO;
}



long ARM_PFGYCSIGVARFIT (long pfId,
                         long curveId,
                         const VECTOR<double>& matCurve,
                         double in_min_meanRev,
                         double in_max_meanRev,
                         double in_min_vol,
                         double in_max_vol,
                         double in_precision_meanRev,
                         double in_precision_vol,
                         long nbMaxIter,
                         ARM_result& result, long objId)
{
    long idx = 0;
    long ARM_REQUEST_ID;
    long ARM_REQUEST_NBPAR = 11;
    CCString stringObjectId;
            
    try
    {
        ARM_CORBA_init ();

        if (CORBA_OBJECT_INTERFACE)
        {
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

            stringObjectId = GetLastCurCellEnvValue ();

            if ( objId != -1 )
            {
                ARM_REQUEST_ID = RPC_SETPFHWSIGVARESTIMATED;
                reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
                REQ_SetLong (reqIn, objId, 0);
                idx = 1;
            }
            else
            {
                ARM_REQUEST_ID = RPC_PFHWSIGVARESTIMATE;
                reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
                idx = 0;
            }

            VECTOR<CCString> matCurve_str;
            for (int i = 0; i < matCurve.size (); i++)
            {
                matCurve_str.push_back(XLDATE2ARMDATE(matCurve[i]));
            }

            /*MSG_printf_message (MSG_INFO, "pfId = %ld\n", pfId);
            MSG_printf_message (MSG_INFO, "curveId = %ld\n", curveId);
            MSG_printf_message (MSG_INFO, "size = %ld\n", matCurve.size ());
            for (i = 0; i < matCurve.size (); i++)
            {
                MSG_printf_message (MSG_INFO, "date[%d] = %s\n", i, (const char*)matCurve_str[i]);
            }
            MSG_printf_message (MSG_INFO, "in_min_meanRev = %lf\n", in_min_meanRev);
            MSG_printf_message (MSG_INFO, "in_max_meanRev = %lf\n", in_max_meanRev);
            MSG_printf_message (MSG_INFO, "in_min_vol = %lf\n", in_min_vol);
            MSG_printf_message (MSG_INFO, "in_max_vol = %lf\n", in_max_vol);
            MSG_printf_message (MSG_INFO, "in_precision_meanRev = %lf\n", in_precision_meanRev);
            MSG_printf_message (MSG_INFO, "in_precision_vol = %lf\n", in_precision_vol);
            MSG_printf_message (MSG_INFO, "nbMaxIter = %ld\n", nbMaxIter);*/

            REQ_SetLong (reqIn, pfId, idx);
            REQ_SetLong (reqIn, curveId, ++idx);
            REQ_SetLong (reqIn, matCurve.size (), ++idx);
            REQ_SetStringVector (reqIn, matCurve_str, ++idx);
            REQ_SetDouble (reqIn, in_min_meanRev, ++idx);
            REQ_SetDouble (reqIn, in_max_meanRev, ++idx);
            REQ_SetDouble (reqIn, in_min_vol, ++idx);
            REQ_SetDouble (reqIn, in_max_vol, ++idx);
            REQ_SetDouble (reqIn, in_precision_meanRev, ++idx);
            REQ_SetDouble (reqIn, in_precision_vol, ++idx);
            REQ_SetLong (reqIn, nbMaxIter, ++idx);
                        
            /*--- CORBA Server Call ---*/

               CORBA_OBJECT_INTERFACE->Send (*reqIn, reqOut);

            result.set (reqOut);

            /*--- requests freeing ---*/
            REQ_Delete (reqOut);
            REQ_Delete (reqIn);

            ARM_RESULT();
        }
    }

    catch(const CORBA::Exception &e)
    {
        CORBA_ERR();
    }

    return ARM_KO;
}



long ARM_PFMODFIT (const CCString& modName,
                   long pfId,
                   double settlement,
                   long zcId,
                   VECTOR<double>& vList,
                   VECTOR<double>& fList,
                   long step,
                   double horizon,
                   long nag_algo,
                   ARM_result& result,
                   long objId)
{
    long idx = 0;
    long ARM_REQUEST_ID;
    long ARM_REQUEST_NBPAR = 10;
    CCString stringObjectId;

    

    try
    {
        /*--- parameters checking ---*/
        if ( vList.size () != fList.size ()) 
        {
            result.setMsg ("ARM_ERR: values and flags array must have same size");
            return ARM_KO;
        }

        ARM_CORBA_init ();

        if (CORBA_OBJECT_INTERFACE)
        {
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

            stringObjectId = GetLastCurCellEnvValue ();

            if ( objId != -1 )
            {
                ARM_REQUEST_ID = RPC_SETPFMODESTIMATED;
                reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
                REQ_SetLong (reqIn, objId, 0);
                idx = 1;
            }
            else
            {
                ARM_REQUEST_ID = RPC_PFMODESTIMATE;
                reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
                idx = 0;
            }

            long tmp_size = vList.size ();
            if(tmp_size == 0)
            {
                vList.push_back (0.0);
                fList.push_back (0.0);
            }

            VECTOR<long> fListTmp;
            for (int i = 0; i < fList.size (); i++)
            {
                fListTmp.push_back ((long)fList[i]);
            }

            REQ_SetString (reqIn, (const char*)modName, idx);
            REQ_SetLong (reqIn, pfId, ++idx);
            REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(settlement), ++idx);
            REQ_SetLong (reqIn, zcId, ++idx);
            REQ_SetLong (reqIn, tmp_size, ++idx);
            REQ_SetDoubleVector (reqIn, vList, ++idx);
            REQ_SetLongVector (reqIn, fListTmp, ++idx);
            REQ_SetLong (reqIn, step, ++idx);
            REQ_SetString (reqIn, (const char*)XLDATE2ARMDATE(horizon), ++idx);
            REQ_SetLong (reqIn, nag_algo, ++idx);
            
            /*--- CORBA Server Call ---*/

               CORBA_OBJECT_INTERFACE->Send (*reqIn, reqOut);

            result.set (reqOut);

            /*--- requests freeing ---*/
            REQ_Delete (reqOut);
            REQ_Delete (reqIn);

            ARM_RESULT();
        }
    }

    catch(const CORBA::Exception &e)
    {
        CORBA_ERR();
    }

    return ARM_KO;
}



long ARM_PFGYCSIGVARGLOBFIT(long pfId, 
                            long curveId,
                            const VECTOR<double>& matCurve,
                            double start_meanRev,
                            const VECTOR<double>& start_vol,
                            ARM_result& result,
                            long objId)
{
    long idx = 0;
    long ARM_REQUEST_ID;
    long ARM_REQUEST_NBPAR = 7;
    CCString stringObjectId;
            
    try
    {
        /*--- parameters checking ---*/

        if ( matCurve.size () != start_vol.size () )
        {
           result.setMsg ("ARM_ERR: Sigma plots array and maturities array must have same size");
           return ARM_KO;
        }

        ARM_CORBA_init ();

        if (CORBA_OBJECT_INTERFACE)
        {
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

            stringObjectId = GetLastCurCellEnvValue ();

            if ( objId != -1 )
            {
               ARM_REQUEST_ID = RPC_SETPFHWSIGVARGLOBMIN;
               reqIn = REQ_NewInitializedRequest(ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
               REQ_SetLong (reqIn, objId, 0);
               idx = 1;
            }
            else
            {
               ARM_REQUEST_ID = RPC_PFHWSIGVARGLOBMIN;
               reqIn = REQ_NewInitializedRequest(ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
               idx = 0;
            }

            VECTOR<CCString> matCurve_str;
            for (int i = 0; i < matCurve.size (); i++)
            {
                matCurve_str.push_back(XLDATE2ARMDATE(matCurve[i]));
            }

            REQ_SetLong(reqIn, pfId, idx);
            REQ_SetLong(reqIn, curveId, ++idx);

            REQ_SetLong(reqIn, matCurve_str.size(), ++idx);
            REQ_SetStringVector(reqIn, matCurve_str, ++idx);

            REQ_SetDouble(reqIn, start_meanRev, ++idx);

            REQ_SetLong(reqIn, start_vol.size(), ++idx);
            REQ_SetDoubleVector(reqIn, start_vol, ++idx);
            
            
            /*--- CORBA Server Call ---*/

               CORBA_OBJECT_INTERFACE->Send(*reqIn, reqOut);

            result.set(reqOut);

            /*--- requests freeing ---*/
            REQ_Delete(reqOut);
            REQ_Delete(reqIn);

            ARM_RESULT();
        }
    }

    catch(const CORBA::Exception &e)
    {
        CORBA_ERR();
    }

    return ARM_KO;
}



long ARM_PFGYCSIGVARPENALFIT(long pfId, 
                             long curveId,
                             const VECTOR<double>& matCurve,
                             double start_meanRev,
                             const VECTOR<double>& start_vol_vect,
                             const VECTOR<double>& penal_vect,
                             const VECTOR<double>& coeff_vect,
                             double theAccuracy,
                             ARM_result& result,
                             long objId)
{
    long idx = 0;
    long ARM_REQUEST_ID;
    long ARM_REQUEST_NBPAR = 10;
    CCString stringObjectId;
            
    try
    {
        /*--- parameters checking ---*/

        if ( matCurve.size () != start_vol_vect.size () )
        {
           result.setMsg ("ARM_ERR: Sizes of Maturities and Volatilities Vectors must be equal");
    
           return(ARM_KO);
        }

        if ( (penal_vect.size ()-1) != start_vol_vect.size () )
        {
           result.setMsg ("ARM_ERR: Sizes of Penalty and Volatilities Vectors must be equal");
    
           return(ARM_KO);
        }

        if ( penal_vect.size () != coeff_vect.size () )
        {
           result.setMsg ("ARM_ERR: Sizes of Penalty and Coeffs. Vectors must be equal");
    
           return(ARM_KO);
        }

        ARM_CORBA_init ();

        if (CORBA_OBJECT_INTERFACE)
        {
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

            stringObjectId = GetLastCurCellEnvValue ();

            if ( objId != -1 )
            {
               ARM_REQUEST_ID = RPC_SETPFHWSIGVARPENALFIT;
               reqIn = REQ_NewInitializedRequest(ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
               REQ_SetLong (reqIn, objId, 0);
               idx = 1;
            }
            else
            {
               ARM_REQUEST_ID = RPC_PFHWSIGVARPENALFIT;
               reqIn = REQ_NewInitializedRequest(ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
               idx = 0;
            }

            VECTOR<CCString> matCurve_str;

            for (int i = 0; i < matCurve.size (); i++)
            {
                matCurve_str.push_back(XLDATE2ARMDATE(matCurve[i]));
            }

            REQ_SetLong(reqIn, pfId, idx);
            REQ_SetLong(reqIn, curveId, ++idx);

            REQ_SetLong(reqIn, matCurve_str.size(), ++idx);
            REQ_SetStringVector(reqIn, matCurve_str, ++idx);

            REQ_SetDouble(reqIn, start_meanRev, ++idx);

            REQ_SetLong(reqIn, start_vol_vect.size(), ++idx);
            REQ_SetDoubleVector(reqIn, start_vol_vect, ++idx);

       
            REQ_SetDoubleVector(reqIn, penal_vect, ++idx);
            REQ_SetDoubleVector(reqIn, coeff_vect, ++idx);

            REQ_SetDouble(reqIn, theAccuracy, ++idx);
    
            
            /*--- CORBA Server Call ---*/

               CORBA_OBJECT_INTERFACE->Send(*reqIn, reqOut);

            result.set(reqOut);

            /*--- requests freeing ---*/
            REQ_Delete(reqOut);
            REQ_Delete(reqIn);

            ARM_RESULT();
        }
    }

    catch(const CORBA::Exception &e)
    {
        CORBA_ERR();
    }

    return ARM_KO;
}



long ARM_PFLOGDECVOLFIT(long pfId, 
                        long curveId,
                        const VECTOR<double>& proba,
                        double theAccuracy,
                        long shapeType,
                        double decay,
                        double slope,
                        double asymptote,
                        const VECTOR<double>& matCurve,
                        const VECTOR<double>& volinit_vect,
                        const VECTOR<double>& coeff_vect,
                        ARM_result& result,
                        long objId)
{
    long idx = 0;
    long ARM_REQUEST_ID;
    long ARM_REQUEST_NBPAR = 13;
    CCString stringObjectId;
            
    try
    {
        /*--- parameters checking ---*/

        if (shapeType)
        {
            if (volinit_vect.size () != proba.size ())
            {
                result.setMsg ("ARM_ERR: Sizes of Initial Vols and Proba Vectors must be equal");
                
                return(ARM_KO);
            }
        }
        else
        {
            if (volinit_vect.size () != matCurve.size ())
            {
                result.setMsg ("ARM_ERR: Sizes of Initial Vols and Maturities Vectors must be equal");
                
                return(ARM_KO);
            }
        }

        if ( volinit_vect.size () != coeff_vect.size () )
        {
            result.setMsg ("ARM_ERR: Sizes of Initial Vols and Coeffs. Vectors must be equal");
    
            return(ARM_KO);
        }

        ARM_CORBA_init ();

        if (CORBA_OBJECT_INTERFACE)
        {
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

            stringObjectId = GetLastCurCellEnvValue ();

            if ( objId != -1 )
            {
                ARM_REQUEST_ID = RPC_SETPFLOGDECVOLFIT;
                reqIn = REQ_NewInitializedRequest(ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
                REQ_SetLong (reqIn, objId, 0);
                idx = 1;
            }
            else
            {
                ARM_REQUEST_ID = RPC_PFLOGDECVOLFIT;
                reqIn = REQ_NewInitializedRequest(ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
                idx = 0;
            }

            VECTOR<CCString> matCurve_str;

            long realMatCurveSize = 0;

            if ( matCurve.size () != 0 )
            {
                realMatCurveSize = matCurve.size ();

                for (int i = 0; i < matCurve.size (); i++)
                {
                    matCurve_str.push_back(XLDATE2ARMDATE(matCurve[i]));
                }
            }
            else
            {
                matCurve_str.push_back("NULLDATE");  
            }

            REQ_SetLong(reqIn, pfId, idx);
            REQ_SetLong(reqIn, curveId, ++idx);

            REQ_SetLong(reqIn, realMatCurveSize, ++idx);

            REQ_SetStringVector(reqIn, matCurve_str, ++idx);

            REQ_SetLong(reqIn, proba.size(), ++idx);
            REQ_SetDoubleVector(reqIn, proba, ++idx);

            REQ_SetDouble(reqIn, theAccuracy, ++idx);

            REQ_SetLong(reqIn, shapeType, ++idx);
            REQ_SetDouble(reqIn, decay, ++idx);
            REQ_SetDouble(reqIn, slope, ++idx);
            REQ_SetDouble(reqIn, asymptote, ++idx);

            REQ_SetDoubleVector(reqIn, volinit_vect, ++idx);
            REQ_SetDoubleVector(reqIn, coeff_vect, ++idx);


            /*--- CORBA Server Call ---*/

            CORBA_OBJECT_INTERFACE->Send(*reqIn, reqOut);

            result.set(reqOut);

            /*--- requests freeing ---*/
            REQ_Delete(reqOut);
            REQ_Delete(reqIn);

            ARM_RESULT();
        }
    }

    catch(const CORBA::Exception &e)
    {
        CORBA_ERR();
    }

    return ARM_KO;
}



long ARM_PFHW2FPENALFIT(long pfId, 
                        long curveId,
                        const VECTOR<double>& initVect,
                        const VECTOR<double>& penal_vect,
                        const VECTOR<double>& coeff_vect,
                        double theAccuracy,
                        ARM_result& result,
                        long objId)
{
    long idx = 0;
    long ARM_REQUEST_ID;
    long ARM_REQUEST_NBPAR = 7;
    CCString stringObjectId;
            
    try
    {
        /*--- parameters checking ---*/

        if ( penal_vect.size () != 5 )
        {
           result.setMsg ("ARM_ERR: Penalties Vector size must be equal to 5");
    
           return(ARM_KO);
        }

        if ( penal_vect.size () != coeff_vect.size () )
        {
           result.setMsg ("ARM_ERR: Sizes of Penalties and weights. Vectors must be equal");
    
           return(ARM_KO);
        }


        ARM_CORBA_init ();

        if (CORBA_OBJECT_INTERFACE)
        {
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

            stringObjectId = GetLastCurCellEnvValue ();

            if ( objId != -1 )
            {
               ARM_REQUEST_ID = RPC_SETPFHW2FPENALFIT;
               reqIn = REQ_NewInitializedRequest(ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
               REQ_SetLong (reqIn, objId, 0);
               idx = 1;
            }
            else
            {
               ARM_REQUEST_ID = RPC_PFHW2FPENALFIT;
               reqIn = REQ_NewInitializedRequest(ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
               idx = 0;
            }

            REQ_SetLong(reqIn, pfId, idx);
            REQ_SetLong(reqIn, curveId, ++idx);
            
            REQ_SetDouble(reqIn, theAccuracy, ++idx);

            REQ_SetLong(reqIn, initVect.size(), ++idx);
            REQ_SetDoubleVector(reqIn, initVect, ++idx);

            REQ_SetDoubleVector(reqIn, penal_vect, ++idx);
            REQ_SetDoubleVector(reqIn, coeff_vect, ++idx);
    
            
            /*--- CORBA Server Call ---*/

               CORBA_OBJECT_INTERFACE->Send(*reqIn, reqOut);

            result.set(reqOut);

            /*--- requests freeing ---*/
            REQ_Delete(reqOut);
            REQ_Delete(reqIn);

            ARM_RESULT();
        }
    }

    catch(const CORBA::Exception &e)
    {
        CORBA_ERR();
    }

    return(ARM_KO);
}



long ARM_PFGYCTWOFACTSIGVARFIT(long pfId,
                               long curveId,
                               const VECTOR<double>& matCurve,
                               double in_a,
                               double in_b,
                               double in_sigmaRatio,
                               double in_rho,
                               double in_min_vol,
                               double in_max_vol,
                               long   in_precision_vol,
                               long nbMaxIter,
                               ARM_result& result, long objId)
{
    long idx = 0;
    long ARM_REQUEST_ID;
    long ARM_REQUEST_NBPAR = 12;
    CCString stringObjectId;
    


    try
    {
        ARM_CORBA_init();

        if (CORBA_OBJECT_INTERFACE)
        {
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

            stringObjectId = GetLastCurCellEnvValue();

            if ( objId != -1 )
            {
               ARM_REQUEST_ID = RPC_SETPFHW2FSIGVARESTIMATED;
               reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
               REQ_SetLong (reqIn, objId, 0);
               idx = 1;
            }
            else
            {
               ARM_REQUEST_ID = RPC_PFHW2FSIGVARESTIMATE;
               reqIn = REQ_NewInitializedRequest (ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
               idx = 0;
            }

            VECTOR<CCString> matCurve_str;

            for (int i = 0; i < matCurve.size (); i++)
            {
                matCurve_str.push_back(XLDATE2ARMDATE(matCurve[i]));
            }

            REQ_SetLong (reqIn, pfId, idx);
            REQ_SetLong (reqIn, curveId, ++idx);
            REQ_SetLong (reqIn, matCurve.size (), ++idx);
            REQ_SetStringVector (reqIn, matCurve_str, ++idx);
            REQ_SetDouble (reqIn, in_a, ++idx);
            REQ_SetDouble (reqIn, in_b, ++idx);
            REQ_SetDouble (reqIn, in_sigmaRatio, ++idx);
            REQ_SetDouble (reqIn, in_rho, ++idx);
            REQ_SetDouble (reqIn, in_min_vol, ++idx);
            REQ_SetDouble (reqIn, in_max_vol, ++idx);
            REQ_SetDouble (reqIn, in_precision_vol, ++idx);
            REQ_SetLong (reqIn, nbMaxIter, ++idx);
                        
            /*--- CORBA Server Call ---*/

               CORBA_OBJECT_INTERFACE->Send (*reqIn, reqOut);

            result.set (reqOut);

            /*--- requests freeing ---*/
            REQ_Delete (reqOut);
            REQ_Delete (reqIn);

            ARM_RESULT();
        }
    }

    catch(const CORBA::Exception &e)
    {
        CORBA_ERR();
    }

    return(ARM_KO);
}



long ARM_INSTLOGDECVOLFIT(long secId, 
						       long curveId,
						       long volBSId,
                               const VECTOR<double>& proba,
                               long LiborType,
							   double strike,
							   double theAccuracy,
							   long shapeType,
							   double decay,
							   double slope,
							   double asymptote,
						       const VECTOR<double>& matCurve,
                               const VECTOR<double>& volinit_vect,
                               const VECTOR<double>& coeff_vect,
                               ARM_result& result,
				               long objId)
{
    long idx = 0;
    long ARM_REQUEST_ID;
    long ARM_REQUEST_NBPAR = 16;
    CCString stringObjectId;
            
    try
    {
        /*--- parameters checking ---*/

        if (shapeType)
        {
            if (volinit_vect.size () != proba.size ())
            {
                result.setMsg ("ARM_ERR: Sizes of Initial Vols and Proba Vectors must be equal");
                
                return(ARM_KO);
            }
        }
        else
        {
            if (volinit_vect.size () != matCurve.size ())
            {
                result.setMsg ("ARM_ERR: Sizes of Initial Vols and Maturities Vectors must be equal");
                
                return(ARM_KO);
            }
        }

        if ( volinit_vect.size () != coeff_vect.size () )
        {
            result.setMsg ("ARM_ERR: Sizes of Initial Vols and Coeffs. Vectors must be equal");
    
            return(ARM_KO);
        }

        ARM_CORBA_init ();

        if (CORBA_OBJECT_INTERFACE)
        {
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

            stringObjectId = GetLastCurCellEnvValue ();

            if ( objId != -1 )
            {
                ARM_REQUEST_ID = RPC_SETINSTLOGDECVOLFIT;
                reqIn = REQ_NewInitializedRequest(ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
                REQ_SetLong (reqIn, objId, 0);
                idx = 1;
            }
            else
            {
                ARM_REQUEST_ID = RPC_INSTLOGDECVOLFIT;
                reqIn = REQ_NewInitializedRequest(ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
                idx = 0;
            }

            VECTOR<CCString> matCurve_str;

            long realMatCurveSize = 0;

            if ( matCurve.size () != 0 )
            {
                realMatCurveSize = matCurve.size ();

                for (int i = 0; i < matCurve.size (); i++)
                {
                    matCurve_str.push_back(XLDATE2ARMDATE(matCurve[i]));
                }
            }
            else
            {
                matCurve_str.push_back("NULLDATE");  
            }

            REQ_SetLong(reqIn, secId, idx);
            REQ_SetLong(reqIn, curveId, ++idx);
            REQ_SetLong(reqIn, volBSId, ++idx);
            REQ_SetLong(reqIn, LiborType, ++idx);
            REQ_SetDouble(reqIn, strike, ++idx);

            REQ_SetLong(reqIn, realMatCurveSize, ++idx);

            REQ_SetStringVector(reqIn, matCurve_str, ++idx);

            REQ_SetLong(reqIn, proba.size(), ++idx);
            REQ_SetDoubleVector(reqIn, proba, ++idx);


            REQ_SetDouble(reqIn, theAccuracy, ++idx);

            REQ_SetLong(reqIn, shapeType, ++idx);
            REQ_SetDouble(reqIn, decay, ++idx);
            REQ_SetDouble(reqIn, slope, ++idx);
            REQ_SetDouble(reqIn, asymptote, ++idx);

            REQ_SetDoubleVector(reqIn, volinit_vect, ++idx);
            REQ_SetDoubleVector(reqIn, coeff_vect, ++idx);


            /*--- CORBA Server Call ---*/

            CORBA_OBJECT_INTERFACE->Send(*reqIn, reqOut);

            result.set(reqOut);

            /*--- requests freeing ---*/
            REQ_Delete(reqOut);
            REQ_Delete(reqIn);

            ARM_RESULT();
        }
    }

    catch(const CORBA::Exception &e)
    {
        CORBA_ERR();
    }

    return ARM_KO;
}



long ARM_PFINSTLOGDECVOLFIT(long pfId,
							long secId, 
							long curveId,
							const VECTOR<double>& proba,
							long UsePFResetDates,
							double theAccuracy,
							long shapeType,
							double decay,
							double slope,
							double asymptote,
                            long VolBSId,
							const VECTOR<double>& matCurve,
							const VECTOR<double>& volinit_vect,
							const VECTOR<double>& coeff_vect,
							ARM_result& result,
							long objId)
{
    long idx = 0;
    long ARM_REQUEST_ID;
    long ARM_REQUEST_NBPAR = 16;
    CCString stringObjectId;
            
    try
    {
        /*--- parameters checking ---*/

        if (shapeType)
        {
            if (volinit_vect.size () != proba.size ())
            {
                result.setMsg ("ARM_ERR: Sizes of Initial Vols and Proba Vectors must be equal");
                
                return(ARM_KO);
            }
        }
        else
        {
            if (volinit_vect.size () != matCurve.size ())
            {
                result.setMsg ("ARM_ERR: Sizes of Initial Vols and Maturities Vectors must be equal");
                
                return(ARM_KO);
            }
        }

        if ( volinit_vect.size () != coeff_vect.size () )
        {
            result.setMsg ("ARM_ERR: Sizes of Initial Vols and Coeffs. Vectors must be equal");
    
            return(ARM_KO);
        }

        ARM_CORBA_init ();

        if (CORBA_OBJECT_INTERFACE)
        {
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqOut = NULL;
            ARM_CorbaRequestModule::ARM_CORBA_REQUEST* reqIn = NULL;

            stringObjectId = GetLastCurCellEnvValue ();

            if ( objId != -1 )
            {
                ARM_REQUEST_ID = RPC_SETPFINSTLOGDECVOLFIT;
                reqIn = REQ_NewInitializedRequest(ARM_REQUEST_ID, ARM_REQUEST_NBPAR + 1);
                REQ_SetLong (reqIn, objId, 0);
                idx = 1;
            }
            else
            {
                ARM_REQUEST_ID = RPC_PFINSTLOGDECVOLFIT;
                reqIn = REQ_NewInitializedRequest(ARM_REQUEST_ID, ARM_REQUEST_NBPAR);
                idx = 0;
            }

            VECTOR<CCString> matCurve_str;

            long realMatCurveSize = 0;

            if ( matCurve.size () != 0 )
            {
                realMatCurveSize = matCurve.size ();

                for (int i = 0; i < matCurve.size (); i++)
                {
                    matCurve_str.push_back(XLDATE2ARMDATE(matCurve[i]));
                }
            }
            else
            {
                matCurve_str.push_back("NULLDATE");  
            }

            REQ_SetLong(reqIn, pfId, idx);
            REQ_SetLong(reqIn, secId, ++idx);
            REQ_SetLong(reqIn, curveId, ++idx);

            REQ_SetLong(reqIn, realMatCurveSize, ++idx);

            REQ_SetStringVector(reqIn, matCurve_str, ++idx);

            REQ_SetLong(reqIn, proba.size(), ++idx);
            REQ_SetDoubleVector(reqIn, proba, ++idx);

            REQ_SetLong(reqIn, UsePFResetDates, ++idx);

            REQ_SetDouble(reqIn, theAccuracy, ++idx);

            REQ_SetLong(reqIn, shapeType, ++idx);
            REQ_SetDouble(reqIn, decay, ++idx);
            REQ_SetDouble(reqIn, slope, ++idx);
            REQ_SetDouble(reqIn, asymptote, ++idx);

            REQ_SetDoubleVector(reqIn, volinit_vect, ++idx);
            REQ_SetDoubleVector(reqIn, coeff_vect, ++idx);

            if (VolBSId != ARM_KO)
                REQ_SetLong(reqIn, VolBSId, ++idx);
            else
                REQ_SetLong(reqIn, 0, ++idx);

            /*--- CORBA Server Call ---*/

            CORBA_OBJECT_INTERFACE->Send(*reqIn, reqOut);

            result.set(reqOut);

            /*--- requests freeing ---*/
            REQ_Delete(reqOut);
            REQ_Delete(reqIn);

            ARM_RESULT();
        }
    }

    catch(const CORBA::Exception &e)
    {
        CORBA_ERR();
    }

    return ARM_KO;
}


/*---- End Of File ----*/
// EOF %M%