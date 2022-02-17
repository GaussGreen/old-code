/*----------------------------------------------------------------------------*
 	bkshift.h
 
	Header for the ARM_BucketShift class, 
    
*----------------------------------------------------------------------------*/
#ifndef _BKSHIFT_H
#define _BKSHIFT_H




#include "armglob.h"



/*----------------------------------------------------------------------------*/



class ARM_BucketShift : public ARM_Object 
{
	private:

         double itsBPShift;           // shift of the yield curve in bp

         double itsBucketStartPeriod; // start period (in yearTerm) of 
                                      // the bucket

         double itsBucketEndPeriod;   // end period (in yearTerm) of the bucket



    
    public:

		ARM_BucketShift(void) 
        {
            Init();

            SetName(ARM_BUCKET_SHIFT);
        }

        //	Constructor 

		ARM_BucketShift(double bpShift, double bucketStartPeriod,
                       double bucketEndPeriod);



       ~ARM_BucketShift(void)
        {
        }

        void Init(void)
        {
            itsBPShift = 0;

            itsBucketStartPeriod = 0;

            itsBucketEndPeriod = 0;
        }

        void BitwiseCopy(const ARM_Object* srcShift)
        {
            ARM_BucketShift* BkShift = (ARM_BucketShift *) srcShift;


            itsBPShift = BkShift->itsBPShift;

            itsBucketStartPeriod = BkShift->itsBucketStartPeriod;

            itsBucketEndPeriod = BkShift->itsBucketEndPeriod;
        }
 
        void Copy(const ARM_Object* srcShift)
        {
            ARM_Object::Copy(srcShift);

            this->BitwiseCopy(srcShift);
        }

        ARM_Object* Clone(void)
        {
            ARM_BucketShift* theClone = new ARM_BucketShift();
 

            theClone->Copy(this);
  
            return(theClone);
        }

        ARM_CLASS_NAME GetRootName(void)
        {
            return(ARM_BUCKET_SHIFT);
        }


		void Print(void)
		{
			printf("\n\n ===> ARM_BucketShift");

            printf("\n BPShift : %lf", itsBPShift);

            printf("\n BucketStartPeriod : %lf", itsBucketStartPeriod);

            printf("\n BucketEndPeriod : %lf", itsBucketEndPeriod);

			printf("\n\n <=== ARM_BucketShift");
		}


        double	GetBPShift(void) const
        {
            return(itsBPShift); 
        }

        void SetBPShift(double bpShift)
		{
			itsBPShift = bpShift;
		}

        double  GetBucketStartPeriod(void) const
        {
            return(itsBucketStartPeriod);
        }
 
        void SetBucketStartPeriod(double BucketStartPeriod)
        {
            itsBucketStartPeriod = BucketStartPeriod;
        }

        double GetBucketEndPeriod(void) const
        {
            return(itsBucketEndPeriod);
        }
 
        void SetBucketEndPeriod(double BucketEndPeriod)
        {
            itsBucketEndPeriod = BucketEndPeriod;
        }
};


#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
