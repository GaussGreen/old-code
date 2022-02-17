#pragma warning(disable : 4541)

#include "firstToBeIncluded.h"
#include <libCCtools++\CCstring.h>

#include "ARM_local_persistent.h"

#include <ARM\libarm\ARM_result.h>
#include <ARM\local_xlarm\ARM_local_interglob.h>

#include <inst\irindex.h>
#include <inst\multiidx.h>

long ARMLOCAL_LIBOR (long liborTypeId,
					 bool ccyIsObject,
					 const CCString& ccyName,
					 long resetFreqId,
					 long payFreqId,
                     long daycount,
					 long intRuleId,
					 ARM_result& result,
					 long objId)
{
	long irindexId;

	ARM_IRIndex* createdIrIndex=NULL;
	ARM_IRIndex* irIndex=NULL;
	ARM_Currency* ccy=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		if(ccyName == "DEFAULT")
		{
			ccy = ARM_DEFAULT_CURRENCY;
		}
		else if (ccyIsObject)
		{
			long ccyId = LocalGetNumObjectId (ccyName);
			ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0)
			{
				result.setMsg ("ARM_ERR: Currency is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			ccy = new ARM_Currency((const char*) ccyName);
		}

		createdIrIndex = new ARM_IRIndex((ARM_INDEX_TYPE) liborTypeId, resetFreqId, 
										 payFreqId, ccy,daycount);
		createdIrIndex->SetIntRule(intRuleId);


        if (!(ccyIsObject)
			&&
			!(ccyName == "DEFAULT"))
		{
		   delete ccy;
		   ccy = NULL;
		}

		if (createdIrIndex == NULL)
		{
			result.setMsg ("ARM_ERR: Libor is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			irindexId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdIrIndex);

			if (irindexId == RET_KO)
			{
				if (createdIrIndex)
					delete createdIrIndex;
				createdIrIndex = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(irindexId);

			return ARM_OK;
		}
		else
		{
			irIndex = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(irIndex, ARM_IRINDEX) == 1)
			{
				if (irIndex)
				{
					delete irIndex;
					irIndex = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdIrIndex, objId);

				return ARM_OK;
			}

			else
			{
				if (createdIrIndex)
					delete createdIrIndex;
				createdIrIndex = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}				
		}	
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdIrIndex)
			delete createdIrIndex;
		createdIrIndex = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_IRINDEX (long dayCountId,
					   long payFreqId,
					   double maturity,
					   long compMethId,
					   long fwdRuleId,
					   long resetTimingId,
					   long resetGap,
					   long payTimingId,
					   long payGap,
					   bool ccyIsObject,
					   const CCString& ccyName,
					   long indexType,
					   long decompFreq,
					   long intRuleId,
					   long resetFreqId,
					   ARM_result& result,
					   long objId)
{
	long irindexId;

	ARM_IRIndex* createdIrIndex=NULL;
	ARM_IRIndex* irIndex=NULL;
	ARM_Currency* ccy=NULL;
	ARM_INDEX_TYPE armIndexType;

	int intRule = 1; // TMP

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	switch(indexType)
	{
		case K_FIXED:
		{
			armIndexType = IDXFIXED;
		};
		break;

		case K_LIBOR1M:
		{
			armIndexType = LIBOR1M;
		};
		break;

		case K_LIBOR2M:
		{
			armIndexType = LIBOR2M;
		};
		break;

		case K_LIBOR3M:
		{
			armIndexType = LIBOR3M;
		};
		break;

		case K_LIBOR6M:
		{
			armIndexType = LIBOR6M;
		};
		break;

		case K_LIBOR1Y:
		{
			armIndexType = LIBOR1Y;
		};
		break;

		case K_EURIBOR1M:
		{
			armIndexType = EURIBOR1M;
		};
		break;

		case K_EURIBOR2M:
		{
			armIndexType = EURIBOR2M;
		};
		break;

		case K_EURIBOR3M:
		{
			armIndexType = EURIBOR3M;
		};
		break;

		case K_EURIBOR6M:
		{
			armIndexType = EURIBOR6M;
		};
		break;

		case K_EURIBOR1Y:
		{
			armIndexType = EURIBOR1Y;
		};
		break;

		case K_PIBOR1M:
		{
			armIndexType = PIBOR1M;
		};
		break;

		case K_PIBOR2M:
		{
			armIndexType = PIBOR2M;
		};
		break;

		case K_PIBOR3M:
		{
			armIndexType = PIBOR3M;
		};
		break;

		case K_PIBOR6M:
		{
			armIndexType = PIBOR6M;
		};
		break;

		case K_PIBOR1Y:
		{
			armIndexType = PIBOR1Y;
		};
		break;

		case K_CMT2:
		{
			armIndexType = CMT2;
		};
		break;

		case K_CMT5:
		{
			armIndexType = CMT5;
		};
		break;

		case K_CMT10:
		{
			armIndexType = CMT10;
		};
		break;

		case K_CMT15:
		{
			armIndexType = CMT15;
		};
		break;

		case K_CMT20:
		{
			armIndexType = CMT20;
		};
		break;

		case K_CMT30:
		{
			armIndexType = CMT30;
		};
		break;

		case K_CMS1:
		{
			armIndexType = CMS1;
		};
		break;

		case K_CMS2:
		{
			armIndexType = CMS2;
		};
		break;

		case K_CMS3:
		{
			armIndexType = CMS3;
		};
		break;

		case K_CMS4:
		{
			armIndexType = CMS4;
		};
		break;

		case K_CMS5:
		{
			armIndexType = CMS5;
		};
		break;

		case K_CMS6:
		{
			armIndexType = CMS6;
		};
		break;

		case K_CMS7:
		{
			armIndexType = CMS7;
		};
		break;

		case K_CMS8:
		{
			armIndexType = CMS8;
		};
		break;

		case K_CMS9:
		{
			armIndexType = CMS9;
		};
		break;

		case K_CMS10:
		{
			armIndexType = CMS10;
		};
		break;

		case K_CMS11:
		{
			armIndexType = CMS11;
		};
		break;

		case K_CMS12:
		{
			armIndexType = CMS12;
		};
		break;

		case K_CMS13:
		{
			armIndexType = CMS13;
		};
		break;

		case K_CMS14:
		{
			armIndexType = CMS14;
		};
		break;

		case K_CMS15:
		{
			armIndexType = CMS15;
		};
		break;

		case K_CMS16:
		{
			armIndexType = CMS16;
		};
		break;

		case K_CMS17:
		{
			armIndexType = CMS17;
		};
		break;

		case K_CMS18:
		{
			armIndexType = CMS18;
		};
		break;

		case K_CMS19:
		{
			armIndexType = CMS19;
		};
		break;

		case K_CMS20:
		{
			armIndexType = CMS20;
		};
		break;

		case K_CMS21:
		{
			armIndexType = CMS21;
		};
		break;

		case K_CMS22:
		{
			armIndexType = CMS22;
		};
		break;

		case K_CMS23:
		{
			armIndexType = CMS23;
		};
		break;

		case K_CMS24:
		{
			armIndexType = CMS24;
		};
		break;

		case K_CMS25:
		{
			armIndexType = CMS25;
		};
		break;

		case K_CMS26:
		{
			armIndexType = CMS26;
		};
		break;

		case K_CMS27:
		{
			armIndexType = CMS27;
		};
		break;

		case K_CMS28:
		{
			armIndexType = CMS28;
		};
		break;

		case K_CMS29:
		{
			armIndexType = CMS29;
		};
		break;

		case K_CMS30:
		{
			armIndexType = CMS30;
		};
		break;

		case K_TEC5:
		{
			armIndexType = TEC5;
		};
		break;

		case K_TEC10:
		{
			armIndexType = TEC10;
		};
		break;

		case K_T4M:
		{
			armIndexType = T4M;
		};
		break;

		case K_T4M_FIXED:
		{
			armIndexType = T4M_FIXED;
		};
		break;

		case K_TAM:
		{
			armIndexType = TAM;
		};
		break;

		case K_TAG:
		{
			armIndexType = TAG;
		};
		break;

        case K_EONIA: // identical to K_TMP
		{
			armIndexType = EONIA; // TMP;
		};
		break;

		default:
		{
			armIndexType = EURIBOR6M;
		}
	}

	CCString msg ("");

	try
	{
		if(ccyName == "DEFAULT")
		{
			ccy = ARM_DEFAULT_CURRENCY;
		}
		else if (ccyIsObject)
		{
			long ccyId = LocalGetNumObjectId (ccyName);
			ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0)
			{
				result.setMsg ("ARM_ERR: Currency is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			ccy = new ARM_Currency((const char*) ccyName);
		}

		if (armIndexType == IDXFIXED)
		{
			createdIrIndex = new ARM_IRIndex(ccy->GetCcyName(),dayCountId);

			if (payFreqId != K_DEF_FREQ)
			{
				createdIrIndex->SetPayFrequency(payFreqId);
			}

			if (resetFreqId != K_DEF_FREQ)
			{
				createdIrIndex->SetResetFrequency(resetFreqId);
			}

			createdIrIndex->SetResetTiming(resetTimingId);
			createdIrIndex->SetPayTiming(payTimingId);

			createdIrIndex->SetCompMeth(compMethId);

			createdIrIndex->SetFwdRule(fwdRuleId);

			if (resetGap != 10000)
				createdIrIndex->SetResetGap(resetGap);

			if (payGap != 10000)
				createdIrIndex->SetPayGap(payGap);

			createdIrIndex->SetDecompFreq(decompFreq);

			// case of float rate maturity == reset frequency (Libor rate)

			if ( maturity == 1.0/createdIrIndex->GetResetFrequency() )
			{
				createdIrIndex->SetYearTerm(maturity);
				createdIrIndex->SetTerm(createdIrIndex->GetResetFrequency());
			}
			else if (maturity > 1.0/createdIrIndex->GetResetFrequency() )
			{
				// case of float rate maturity > reset frequency 
				// (6m Libor reset every 3m)

				createdIrIndex->SetYearTerm(maturity);
				createdIrIndex->SetTerm(int(1.0/maturity));
			}
			else if (maturity > 1.0)
			{
				// case of float rate maturity > 1 year (CMS, CMT)

				createdIrIndex->SetYearTerm(maturity);
				createdIrIndex->SetTerm(K_ANNUAL);
			}

			createdIrIndex->SetIntRule(intRuleId);
		}
		else
		{

			createdIrIndex = new ARM_IRIndex(armIndexType,
											 resetFreqId,
											 payFreqId,
											 ccy);

			if ( !(( armIndexType >= K_CMS1) && ( armIndexType <= K_CMS30 ))
                && createdIrIndex->IsLiborIndex() )
			{
				if (dayCountId != -1)
				{
					createdIrIndex->SetDayCount(dayCountId);
				}

				createdIrIndex->SetResetTiming(resetTimingId);
				createdIrIndex->SetPayTiming(payTimingId);

				createdIrIndex->SetCompMeth(compMethId);

				createdIrIndex->SetFwdRule(fwdRuleId);

				if (resetGap != 10000)
					createdIrIndex->SetResetGap(resetGap);

				if (payGap != 10000)
					createdIrIndex->SetPayGap(payGap);

				createdIrIndex->SetDecompFreq(decompFreq);
				
				// case of float rate maturity == reset frequency (Libor rate)

				if ( maturity == 1.0/createdIrIndex->GetResetFrequency() )
				{
					createdIrIndex->SetYearTerm(maturity);
					createdIrIndex->SetTerm(createdIrIndex->GetResetFrequency());
				}
				else if (maturity > 1.0/createdIrIndex->GetResetFrequency() )
				{
					// case of float rate maturity > reset frequency 
					// (6m Libor reset every 3m)

					createdIrIndex->SetYearTerm(maturity);
					createdIrIndex->SetTerm(int(1.0/maturity));
				}
				else if (maturity > 1.0)
				{
					// case of float rate maturity > 1 year (CMS, CMT)

					createdIrIndex->SetYearTerm(maturity);
					createdIrIndex->SetTerm(K_ANNUAL);
				}

				createdIrIndex->SetIndexStyle(IN_ARREARS);

				if ( (createdIrIndex->GetResetTiming() == K_ADVANCE)
					&& (createdIrIndex->GetPayTiming() == K_ARREARS)
					&& (createdIrIndex->GetResetFrequency() == createdIrIndex->GetPayFrequency())
				   )
				{
					createdIrIndex->SetIndexStyle(VANILLA);
				}

				createdIrIndex->SetIntRule(intRuleId);


			}
			else
			{
				if ( dayCountId == -1 )
				{
					dayCountId = createdIrIndex->GetDayCount();
				}

				if ( resetFreqId == K_DEF_FREQ )
				{
					resetFreqId = createdIrIndex->GetResetFrequency();
				}

				if ( payFreqId == K_DEF_FREQ )
				{
					payFreqId = createdIrIndex->GetPayFrequency();
				}

				if ( resetGap == 10000 )
				{
					resetGap = createdIrIndex->GetResetGap();
				}

				if ( payGap == 10000 )
				{
					payGap = createdIrIndex->GetPayGap();
				}
				
				ARM_INDEX_TYPE Libor_Type = ccy->GetVanillaIndexType();

				createdIrIndex->Set(armIndexType, Libor_Type, resetFreqId, compMethId, ccy);

				createdIrIndex->Set(dayCountId, resetFreqId, payFreqId,
									maturity, compMethId, fwdRuleId, intRuleId,
									resetTimingId, resetGap, payTimingId, payGap,
									ccy, armIndexType, decompFreq);

                //createdIrIndex->SetIndexType(armIndexType);

				


            }
		}

        if (!(ccyIsObject)
			&&
			!(ccyName == "DEFAULT"))
		{
		   delete ccy;
		   ccy = NULL;
		}

		if (createdIrIndex == NULL)
		{
			result.setMsg ("ARM_ERR: IrIndex is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			irindexId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdIrIndex);

			if (irindexId == RET_KO)
			{
				if (createdIrIndex)
					delete createdIrIndex;
				createdIrIndex = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(irindexId);

			return ARM_OK;
		}
		else
		{
			irIndex = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(irIndex, ARM_IRINDEX) == 1)
			{
				if (irIndex)
				{
					delete irIndex;
					irIndex = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdIrIndex, objId);

				return ARM_OK;
			}
			else
			{
				if (createdIrIndex)
					delete createdIrIndex;
				createdIrIndex = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdIrIndex)
			delete createdIrIndex;
		createdIrIndex = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_FixedIndex(long dayCountId,
						 const CCString& ccy,
						 ARM_result& result,
						 long objId)
{
	long irindexId;

	ARM_IRIndex* createdIrIndex=NULL;
	ARM_IRIndex* irIndex=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
	    char* ccyChar = ccy;
        createdIrIndex = new ARM_IRIndex(ccyChar, dayCountId);
        delete ccyChar;
      

		if (createdIrIndex == NULL)
		{
			result.setMsg ("ARM_ERR: Fixed Index is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			irindexId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdIrIndex);

			if (irindexId == RET_KO)
			{
				if (createdIrIndex)
					delete createdIrIndex;
				createdIrIndex = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(irindexId);

			return ARM_OK;
		}
		else
		{
			irIndex = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(irIndex, ARM_IRINDEX) == 1)
			{
				if (irIndex)
				{
					delete irIndex;
					irIndex = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdIrIndex, objId);

				return ARM_OK;
			}

			else
			{
				if (createdIrIndex)
					delete createdIrIndex;
				createdIrIndex = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}				
		}	
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdIrIndex)
			delete createdIrIndex;
		createdIrIndex = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_IRINDEX_MONEY_MARKET (const CCString& mmTerm,
									const CCString& ccy,
									ARM_result& result,
									long objId)
{
	long irindexId;

	ARM_IRIndex* createdIrIndex = NULL;
	ARM_IRIndex* irIndex = NULL;

	ARM_Currency* ccyName = NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	char* sCCY;
	char* sTerm;

	try
	{
		sCCY = ccy;
		sTerm = mmTerm;

		if (sCCY != NULL)
			ccyName = new ARM_Currency(sCCY);

		createdIrIndex = new ARM_IRIndex(sTerm, ccyName);

		if (ccyName)
			delete ccyName;
		ccyName = NULL;

		if (sCCY)
			delete sCCY;
		sCCY = NULL;

		if (sTerm)
			delete sTerm;
		sTerm = NULL;

		if (createdIrIndex == NULL)
		{
			result.setMsg ("ARM_ERR: IR Index is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			irindexId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdIrIndex);

			if (irindexId == RET_KO)
			{
				if (createdIrIndex)
					delete createdIrIndex;
				createdIrIndex = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(irindexId);

			return ARM_OK;
		}
		else
		{
			irIndex = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(irIndex, ARM_IRINDEX) == 1)
			{
				if (irIndex)
				{
					delete irIndex;
					irIndex = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdIrIndex, objId);

				return ARM_OK;
			}

			else
			{
				if (createdIrIndex)
					delete createdIrIndex;
				createdIrIndex = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}				
		}	
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdIrIndex)
			delete createdIrIndex;
		createdIrIndex = NULL;

		if (ccyName)
			delete ccyName;
		ccyName = NULL;

		if (sCCY)
			delete sCCY;
		sCCY = NULL;

		if (sTerm)
			delete sTerm;
		sTerm = NULL;

		ARM_RESULT();
	}
}



long ARMLOCAL_CMS(long CMSType,
				  long liborTypeId,
				  bool ccyIsObject,
				  const CCString& ccyName,
				  ARM_result& result,
				  long objId)
{
	long irindexId;

	ARM_IRIndex* createdIrIndex=NULL;
	ARM_IRIndex* irIndex=NULL;
	ARM_Currency* ccy=NULL;

	if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
	{
		result.setMsg ("ARM_ERR: Pb with accessing objects");
		return ARM_KO;
	}

	CCString msg ("");

	try
	{
		if(ccyName == "DEFAULT")
		{
			ccy = ARM_DEFAULT_CURRENCY;
		}
		else if (ccyIsObject)
		{
			long ccyId = LocalGetNumObjectId (ccyName);
			ccy = (ARM_Currency *) LOCAL_PERSISTENT_OBJECTS->GetObject(ccyId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(ccy, ARM_CURRENCY) == 0)
			{
				result.setMsg ("ARM_ERR: Currency is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			ccy = new ARM_Currency((const char*) ccyName);
		}
		
		createdIrIndex = new ARM_IRIndex((ARM_INDEX_TYPE) CMSType,
										 (ARM_INDEX_TYPE) liborTypeId,
										 K_DEF_FREQ,
										 K_ANNUAL,
										 ccy);

		if (createdIrIndex == NULL)
		{
			result.setMsg ("ARM_ERR: Index is null");
			return ARM_KO;
		}

		if (objId == -1)
		{
			CREATE_GLOBAL_OBJECT();
			
			irindexId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdIrIndex);

			if (irindexId == RET_KO)
			{
				if (createdIrIndex)
					delete createdIrIndex;
				createdIrIndex = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(irindexId);

			return ARM_OK;
		}
		else
		{
			irIndex = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(irIndex, ARM_IRINDEX) == 1)
			{
				if (irIndex)
				{
					delete irIndex;
					irIndex = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)createdIrIndex, objId);

				result.setLong(objId);

				return ARM_OK;
			}

			else
			{
				if (createdIrIndex)
					delete createdIrIndex;
				createdIrIndex = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}				
		}	
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		if (createdIrIndex)
			delete createdIrIndex;
		createdIrIndex = NULL;

		ARM_RESULT();
	}
}


long ARMLOCAL_MultiIrindex(VECTOR<CCString>& irindexVect,
						   VECTOR<double>& weightVect,
						   ARM_result& result, 
						   long objId)
{
	long MultiIrindexId;

	int irindexVec_size = irindexVect.size ();
	int weightVect_size = weightVect.size ();

	//verify the size
	if( irindexVec_size!=weightVect_size || irindexVec_size < 2 )
	{
		result.setMsg ("ARM_ERR: check your matrix dimension");
		return ARM_KO;		
	}

	ARM_Vector* vWeight = new ARM_Vector(weightVect_size);
	ARM_IRIndex* Idx[20];
	for(int i=0; i<20; i++)
	{
		Idx[i] = NULL;
	}

	ARM_MultiIndex* multiIdx= NULL;
	ARM_MultiIndex* multiIdxOrg = NULL;

	CCString msg ("");

	try
	{
		char vIrindex[200][20];
		
		long vIrindexIds[200];

		if (CHECK_GLOBAL_OBJECT(LOCAL_PERSISTENT_OBJECTS) == 0)
		{
			result.setMsg ("ARM_ERR: Pb with accessing objects");
			return ARM_KO;
		}
		int i;
		for(i=0; i<irindexVec_size; i++)
		{
			sprintf(vIrindex[i],(const char*)irindexVect[i]);
			vWeight->Elt(i) = weightVect[i];
			vIrindexIds[i] = LocalGetNumObjectId(irindexVect[i]);
			Idx[i] = new ARM_IRIndex();
			Idx[i] = (ARM_IRIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject((long)vIrindexIds[i])->Clone();
			
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(Idx[i], ARM_IRINDEX) == 0)
			{
				result.setMsg ("ARM_ERR: Index is not of a good type");
				return ARM_KO;
			}
		}
		
		multiIdx = new ARM_MultiIndex(irindexVec_size, Idx, vWeight);
		
		for( i=0; i<20; i++)
		{
			if(Idx[i])
				delete Idx[i];
			Idx[i] = NULL;
		}

		if(vWeight = NULL)
			delete vWeight;
		vWeight = NULL;

		if ( objId == -1 )
		{
			CREATE_GLOBAL_OBJECT();

			MultiIrindexId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object *) multiIdx);

			if ( MultiIrindexId == RET_KO )
			{
				if (multiIdx)
					delete multiIdx;
				multiIdx = NULL;

				result.setMsg("ARM_ERR: Pb with inserting object");				
				
                return(ARM_KO);
			}

			result.setLong(MultiIrindexId);

			return ARM_OK;			
		}
		else
		{
			multiIdxOrg = (ARM_MultiIndex *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(multiIdxOrg, ARM_MULTIINDEX) == 1)
			{
				if (multiIdxOrg)
				{
					delete multiIdxOrg;
					multiIdxOrg = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*) multiIdx, objId);

				return(ARM_OK);
			}
			else
			{
				if (multiIdx)
					delete multiIdx;
				multiIdx = NULL;
				
				result.setMsg ("ARM_ERR: previous object is not of a good type");

				return ARM_KO;
			}
		}

	}
		
	catch(Exception& x)
	{
		for(int i=0; i<20; i++)
		{
			if(Idx[i])
				delete Idx[i];
			Idx[i] = NULL;
		}

		if(vWeight = NULL)
			delete vWeight;
		vWeight = NULL;

		if (multiIdx)
			delete multiIdx;
		multiIdx = NULL;

		ARM_RESULT();
	}

}
