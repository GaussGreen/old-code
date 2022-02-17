#pragma warning(disable : 4541)


#ifndef ARMLOCAL_BARRIER_H
#define ARMLOCAL_BARRIER_H


long ARMLOCAL_BARRIER (long underlyingId,
			    	    long tAssetId,
				        long xStyleId,
				        long refValId1,
						long refValId2,
				        long upDownDouble,
				        long inOut,
				        long triggerVar,
				        double rebate,
						int isInArrear,
				        ARM_result& result,
				        long objId = -1);


long ARMLOCAL_CONSTBARRIER (long underlyingId,
					   long tAssetId,
					   double maturity,
					   double barrier1,
					   double barrier2,
					   long upDownDouble,
					   long inOut,
					   long triggerVar,
					   double rebate,
					   double firstX,
					   int isInArrear,
					   ARM_result& result,
					   long objId =-1);

#endif