#ifndef _REFVALUE_H
#define _REFVALUE_H

#include "armglob.h"
#include "linalg.h"
#include "containr.h"
#include "dates.h"

class ARM_Security;
class ARM_ZeroCurve;
class ARM_Model;

/*----------------------------------------------------------------------------*/
class ARM_ReferenceValue : public ARM_Object 
{
    private:

        int itsValueType;    // 1 = price, 0 = yield
        /*  1 = if value must be converted from itsValueType 
		    when comparing value or when computing payoff, 0 otherwise  */
        
        int ConvertOrNot;    
        /* reference value calculation method
           K_CONSTANT    0 : constant value
           K_LININTERPOL 1 : linear interpolation meth from the discrete values
           K_DISCRETE    2 : values known only at discrete dates
           101 = average, 102 = MAX or MIN using model methods 
           1, 2,... other methods (codes to be implemented) */

        int itsCalculationMethod;   

        // Discrete Julian dates used for x
        ARM_Vector* itsDiscreteDates;    
        ARM_Vector* itsDiscreteValues2; // Discrete values used for the second dimension
        ARM_Vector* itsDiscreteValues; // Discrete values used for interpolation    
          
        // the following variables are used for MAX, MIN or Average functions
        int MaxOrMin;          // 1 = MAX, 0 = MIN

        // regularity of calculation
        // 1 = annual, 2 = semiannual, 4 = quarterly, 12 = monthly
        int itsCalculationFrequency;  
                                    
        // specify the day of the week or day of the month 
        // on which calculation is made
        int itsCalcDay;        

        // specify the compounding method used in zero price calculation
        int itsCompoundMeth;        

        double itsZeroRate;        

        int itsExtrapolMeth; // Default 0 : linear , 1 : Constant 

    public:

        ARM_ReferenceValue(void) 
        {
            Init();
        }

        // Constructor for constant reference value
        ARM_ReferenceValue(double value, int valueType = 1, 
                           int conversion = 0); 


        // Constructor for reference value obtained from lin interpol
        ARM_ReferenceValue(ARM_Vector* dates,
                           ARM_Vector* values,
                           int valueType = 1,
                           int conversion = 0,
						   int InterpolMeth = K_LININTERPOL_REF);


        ARM_ReferenceValue(ARM_Vector* dates, ARM_Vector* values, 
                           ARM_Vector* values2,
                           int valueType = 1, int conversion = 0); 

        ARM_ReferenceValue(ARM_Matrix* refValue);

        ARM_ReferenceValue(ARM_Date& endDate, double endValue,
                           double zeroRate, int compMeth);        

        // Constructor for reference value obtained from 
        // the payment dates of a security and a vector of values 
        ARM_ReferenceValue(ARM_Security* sec, ARM_Vector* values);


        ARM_ReferenceValue(const ARM_ReferenceValue& refValue);

        ARM_ReferenceValue& operator = (const ARM_ReferenceValue& refValue);

       ~ARM_ReferenceValue(void)
        {
            if (itsDiscreteValues)
               delete itsDiscreteValues;
            
            itsDiscreteValues = NULL;

            if (itsDiscreteValues2)
               delete itsDiscreteValues2;

            itsDiscreteValues2 = NULL;

            if (itsDiscreteDates)
               delete itsDiscreteDates;

            itsDiscreteDates = NULL;
        }

        void Init(void)
        {
            SetName(ARM_REFERENCE_VALUE);
            itsValueType			= K_PRICE;
            ConvertOrNot			= 0;
            itsCalculationMethod	= 0;
            MaxOrMin				= -10000;
            itsCalculationFrequency	= -10000;
            itsCalcDay				= -10000;
            itsDiscreteValues		= NULL;
            itsDiscreteValues2		= NULL;
            itsDiscreteDates		= NULL;
            itsCompoundMeth			= K_COMP_CONT;        
            itsZeroRate				= 0.0;        
            itsExtrapolMeth			= 0;
        }

        void BitwiseCopy(const ARM_Object* srcRef)
        {
            ARM_ReferenceValue* ref = (ARM_ReferenceValue *) srcRef;
			
			if (!ref) return;
            
			itsValueType			= ref->itsValueType;
            ConvertOrNot			= ref->ConvertOrNot;
            itsCalculationMethod	= ref->itsCalculationMethod;
            MaxOrMin				= ref->MaxOrMin;
            itsCalculationFrequency	= ref->itsCalculationFrequency;
            itsCalcDay				= ref->itsCalcDay;
            itsCompoundMeth			= ref->itsCompoundMeth;        
            itsZeroRate				= ref->itsZeroRate;        
            itsExtrapolMeth			= ref->itsExtrapolMeth;

            if (itsDiscreteValues)
               delete itsDiscreteValues;
			itsDiscreteValues		= ref->itsDiscreteValues ? (ARM_Vector *) ref->itsDiscreteValues->Clone() : NULL;

            if (itsDiscreteValues2)
               delete itsDiscreteValues2;
			itsDiscreteValues2		= ref->itsDiscreteValues2 ? (ARM_Vector *) ref->itsDiscreteValues2->Clone() : NULL;

            if (itsDiscreteDates)
               delete itsDiscreteDates;
			itsDiscreteDates		= ref->itsDiscreteDates ? (ARM_Vector *) ref->itsDiscreteDates->Clone() : NULL;
        }
 
        void Copy(const ARM_Object* srcRef)
        {
            ARM_Object::Copy(srcRef);
            this->BitwiseCopy(srcRef);
        }

        ARM_Object* Clone(void)
        {
            ARM_ReferenceValue* theClone = new ARM_ReferenceValue();
            theClone->Copy(this);
            return(theClone);
        }

        ARM_ReferenceValue& operator += (double val)
        {
            int sz = itsDiscreteValues->GetSize();

            for (int i = sz-1; i >= 0; i--)
            {
                itsDiscreteValues->Elt(i) += val;
            }

            return(*this);
        }

        inline ARM_CLASS_NAME GetRootName(void)
        {
            return(ARM_REFERENCE_VALUE);
        }

        int GetExtrapolMeth(void)
        {
            return(itsExtrapolMeth);
        }

        void SetExtrapolMeth(int meth)
        {
            itsExtrapolMeth = meth;
        }
        
		long GetSize(void) const
        {
            return (itsDiscreteDates == NULL)? 0 : itsDiscreteDates->GetSize();
        }
	
		inline size_t size() const { return (itsDiscreteDates == NULL)? 0 : itsDiscreteDates->GetSize();}

        long NbCurves(void) const
        {
            long nbcurves = 0;

            if ( itsDiscreteDates == NULL ) 
               return nbcurves;
            else
            {
               if ( itsDiscreteValues == NULL )
                  return(nbcurves);
               else
               {
                  nbcurves += 1;

                  if ( itsDiscreteValues2 != NULL )
                  {
                     nbcurves += 1;

                     return(nbcurves);
                  }
                  else
                     return(nbcurves);
               }
            }
        }

        // virtual double CptReferenceValue(const ARM_Date& date) const;
        virtual double CptReferenceValue(double JulianDate) const;
		double CptReferenceValue(const ARM_Date&date) const
		{
			return CptReferenceValue(date.GetJulian()); 
		}

		// virtual 
		// ARM_Vector* CptReferenceValues(ARM_Date& date);
        void InterpolateRefValue(ARM_Vector* date);

        ARM_ReferenceValue* CptReferenceValues(ARM_Vector* dates);

		double ComputePrice(ARM_Model* model, ARM_Date* dateValo = NULL);

        ARM_Vector* GetDiscreteValues(int index = 0) const
        {
            if ( index == 0 )
               return(itsDiscreteValues);
            else
               return(itsDiscreteValues2);
        }

        void SetDiscreteValues(ARM_Vector* discreteValues, int index = 0)
        {
            if (index == 0)
            {
                if (itsDiscreteValues)
                {
                    delete itsDiscreteValues;

                    itsDiscreteValues =NULL;
                }

                if (discreteValues)
                {
                   itsDiscreteValues = discreteValues;
                }
            }
            else
            {
                if (itsDiscreteValues2)
                {
                    delete itsDiscreteValues2;

                    itsDiscreteValues2 =NULL;
                }

                if (discreteValues)
                {
                   itsDiscreteValues2 = discreteValues;
                }
            }
        }

        void SetXYValue(int index, double X,double Y, int dim = 0);
        void SetXYValue(int index, double X,ARM_Vector* Y);

        inline void insert(double x,double y)
        {
            if (itsDiscreteDates)
            {
                itsDiscreteDates->insert(x);
                itsDiscreteValues->insert(y);
            }
            else
			{
				// for a const refValue !!
				// first Time
				SetCalcMethod(K_DISCRETE_REF);
				itsDiscreteDates = new ARM_Vector(1,x);
                
				if (itsDiscreteValues)
				{
					delete itsDiscreteValues;
				}
                itsDiscreteValues = new ARM_Vector(1,y);
			}
        }

        inline void update(double x,double y)
        {
            if (itsDiscreteDates == NULL)
            {
				insert(x, y);
            }
            else
			{
				int j = itsDiscreteDates->find(x);

				if (j != -1)
				{
					// on update la donnée
					itsDiscreteValues->Elt(j) = y;
				}
				else
				{
					// la donnée n'existe pas pour la date correspondante
					// on l'insert au bon endroit
					ARM_Vector newDiscreteDates(itsDiscreteDates->GetSize()+1);
					ARM_Vector* newDiscreteValues = new ARM_Vector(itsDiscreteValues->GetSize()+1);

					int i = 0;
					bool trouve = false;
					while (i < itsDiscreteDates->GetSize())
					{
						if ( x > itsDiscreteDates->Elt(i)) 
						{
							newDiscreteDates.Elt(i) = itsDiscreteDates->Elt(i);
							newDiscreteValues->Elt(i) = itsDiscreteValues->Elt(i);
							i++;
						}
						else
						{
							if (trouve == true)
							{
								newDiscreteDates.Elt(i+1) = itsDiscreteDates->Elt(i);
								newDiscreteValues->Elt(i+1) = itsDiscreteValues->Elt(i);
								i++;
							}
							else
							{
								newDiscreteDates.Elt(i) = x;
								newDiscreteValues->Elt(i) = y;
								trouve = true;
							}
						}
					}
					if (trouve == false)
					{
						newDiscreteDates.Elt(newDiscreteDates.GetSize()-1) = x;
						newDiscreteValues->Elt(newDiscreteValues->GetSize()-1) = y;
					}
					SetDiscreteDates(&newDiscreteDates);
					SetDiscreteValues(newDiscreteValues);
				}
			}
        }

        double Interpolate(ARM_Date& date, int dim = 0);
        double Interpolate(double JulianDate, int dim = 0);

        ARM_Vector* InterpolateMulti(double Juliandate);
        ARM_Vector* InterpolateMulti(ARM_Date& date);
        void View(char* id, FILE* ficOut);

		void DisplayScheduleDates(int datesType, int viewInitExch,
                                 char* id = NULL, FILE* fOut = NULL);

        void DisplayScheduleValues(int valuesType,
                                  char* id = NULL, FILE* fOut = NULL);

		int IsConstant(void) const;

        inline int GetCompoundMeth(void)
        {
            return(itsCompoundMeth);
        }

        /*inline int IsConstant(void) const
        {
            return( itsCalculationMethod == K_CONSTANT_REF );
        }*/

        inline void SetCompoundMeth(int compMeth)
        {
            itsCompoundMeth = compMeth;
        }

        inline double GetZeroRate(void)
        {
            return(itsZeroRate);
        }

        inline void SetZeroRate(double zRate)
        {
            itsZeroRate = zRate;
        }

        inline int GetCalcDay(void)
        {
            return(itsCalcDay);
        }

        inline void SetCalcDay(int day)
        {
            itsCalcDay = day;
        }

        inline int GetMaxOrMin(void)
        {
            return(MaxOrMin);
        }

        inline void SetMaxOrMin(int maxMin)
        {
            MaxOrMin = maxMin;
        }

        inline int GetValueType(void) const
        {
            return(itsValueType); 
        }

        inline void SetValueType(int vType)
        {
            itsValueType = vType;
        }

        inline int GetConvFlag(void) const
        {
            return(ConvertOrNot); 
        }

        inline void SetConvFlag(int conv)
        {
            ConvertOrNot = conv;
        }

        inline int GetCalcMethod(void) const
        {
            return(itsCalculationMethod); 
        }

        inline void SetCalcMethod(int calcMeth)
        {
            itsCalculationMethod = calcMeth;
        }

        inline int GetCalcFrequency(void) const
        {
            return(itsCalculationFrequency); 
        }

        inline void SetCalcFrequency(int calcFreq)
        {
            itsCalculationFrequency = calcFreq;
        }


        ARM_Vector* GetDiscreteDates(void) const 
        {    
            return(itsDiscreteDates); 
        }

        void SetDiscreteDates(ARM_Vector* discreteDates) 
        {
            if (itsDiscreteDates)
            {
               delete itsDiscreteDates;

               itsDiscreteDates =NULL;
            }

            if (discreteDates) 
            {
               itsDiscreteDates = (ARM_Vector *) discreteDates->Clone();
            }
        }

        ARM_ReferenceValue& operator *= (double scalar)
        {
            (*itsDiscreteValues) *= scalar;

            return(*this);
        }

		 ARM_ReferenceValue& operator /= (double scalar)
        {
            (*itsDiscreteValues) /= scalar;

            return(*this);
        }

        ARM_ReferenceValue operator + (ARM_ReferenceValue& refVal2);

		ARM_ReferenceValue operator * (ARM_ReferenceValue& refVal2);

		ARM_ReferenceValue operator / (ARM_ReferenceValue& refVal2);

        ARM_ReferenceValue* CreateStartResetAdjust(int spot, char* cur);
};




#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
