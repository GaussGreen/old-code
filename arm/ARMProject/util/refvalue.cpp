
#include "refvalue.h"
//#include "interpol.h"
//#include "security.h"
//#include "fromto.h"
//#include "merge.h"




/*----------------------------------------------------------------------------*
    Constructor for constant reference value
*----------------------------------------------------------------------------*/ 
ARM_ReferenceValue::ARM_ReferenceValue(double value, 
                                       int valueType, int conversion)
{
    Init();

    itsValueType = valueType;

    ConvertOrNot = conversion;

    itsDiscreteDates = NULL;

    itsDiscreteValues = new ARM_Vector(1, value);
    
    itsCalculationMethod = K_CONSTANT_REF;

    MaxOrMin = -10000;

    itsCalculationFrequency = -10000;

    itsCalcDay = -10000;
}



/*----------------------------------------------------------------------------*
    Constructor for reference value obtained by linear interpolation
*----------------------------------------------------------------------------*/ 

ARM_ReferenceValue::ARM_ReferenceValue(ARM_Vector* dates, 
                                       ARM_Vector* values,
                                       int valueType, 
                                       int conversion,
									   int InterpolMeth)
{
    Init();

    itsValueType = valueType;

    ConvertOrNot = conversion;

    itsDiscreteDates = dates;

    itsDiscreteValues = values;

    itsCalculationMethod = InterpolMeth;

    MaxOrMin = -10000;

    itsCalculationFrequency = -10000;

    itsCalcDay = -10000;
}


/*-----------------------------------------------------------------------------*
   Constructor for linear interpolation in 2 dimensions
*-----------------------------------------------------------------------------*/

ARM_ReferenceValue::ARM_ReferenceValue(ARM_Vector* dates,
                                       ARM_Vector* values,
                                       ARM_Vector* values2,
                                       int valueType, int conversion)
{
    Init();

    itsValueType = valueType;
    ConvertOrNot = conversion;
    itsDiscreteDates = dates;
    itsDiscreteValues = values;
    itsDiscreteValues2 = values2;
    itsCalculationMethod = K_LININTERPOL_REF;
    MaxOrMin = -10000;
    itsCalculationFrequency = -10000;
    itsCalcDay = -10000;
}


ARM_ReferenceValue::ARM_ReferenceValue(ARM_Matrix* refValue)
{
    Init();
    //if ( refValue == NULL)
    //   throw() /*Exception(__LINE__, __FILE__, ERR_PRICING_PB,
    //                   "The Matrix isn't allocated");*/
    //else
    //{
    //   itsDiscreteDates = refValue->GetColumn(0);

    //   if ( refValue->GetNumCols() >= 2 )
    //   {
    //      ARM_Vector* col1 = refValue->GetColumn(1);

    //      if ( col1 != NULL )
    //         itsDiscreteValues = col1;
    //      
    //   }
    //   else
    //   {
    //      throw() /*Exception(__LINE__, __FILE__, ERR_PRICING_PB,
    //                      "The Matrix doesn't have enough dimensions")*/;
    //   }

    //   if ( refValue->GetNumCols() >= 3 )
    //   {
    //      ARM_Vector* col2 = refValue->GetColumn(2);

    //      if ( col2 != NULL )
    //         itsDiscreteValues2 = col2;
    //   }

    //   itsCalculationMethod = K_LININTERPOL_REF;
    //}
}


/*----------------------------------------------------------------------------*
    Constructor for reference value computed as zero coupon price
    with given rate
ARM_ReferenceValue::ARM_ReferenceValue(ARM_Date& endDate,
                                       double endValue,
                                       ARM_Model* ycModel)
{
    Init();

    SetName(ARM_REFERENCE_VALUE);

    itsCalculationMethod = K_ZEROCOUPON_REF;

    itsModel = ycModel;

    itsDiscreteDates = new ARM_Vector(1, endDate.GetJulian());

    itsDiscreteValues = new ARM_Vector(1, endValue);
}
*----------------------------------------------------------------------------*/ 



/*----------------------------------------------------------------------------*
    Constructor for reference value from swap leg 
    payment dates and value vector
*----------------------------------------------------------------------------*/ 
ARM_ReferenceValue::ARM_ReferenceValue(ARM_Security* sec, 
                                       ARM_Vector* values)
{
    Init();

    //if ( values->GetSize() != sec->GetNumFlows() )
    //{
    //   throw()/* Exception(__LINE__, __FILE__, ERR_PRICING_PB,
    //   "The Values vector and the Security cash flows must have the same size")*/;
    //}
    //else
    //{
    //   itsValueType = K_PRICE;

    //   ConvertOrNot = 0;

    //   itsDiscreteDates = (ARM_Vector *) sec->GetPaymentDates()->Clone();

    //   itsDiscreteValues = values;
    //
    //   itsCalculationMethod = K_DISCRETE_REF;

    //   MaxOrMin = -10000;

    //   itsCalculationFrequency = -10000;

    //   itsCalcDay = -10000;
    //}
}



/*----------------------------------------------------------------------------*
    Assignment operator
/*----------------------------------------------------------------------------*/

ARM_ReferenceValue& ARM_ReferenceValue::operator = (const ARM_ReferenceValue& refValue)
{
    (*this).ARM_Object::operator = (refValue);

    this->BitwiseCopy(&refValue); 

    return(*this);
}



/*----------------------------------------------------------------------------*
    Constructor    (copy).
*----------------------------------------------------------------------------*/

ARM_ReferenceValue::ARM_ReferenceValue(const ARM_ReferenceValue& refValue)
                                      : ARM_Object(refValue)
{
    Init();

    this->BitwiseCopy(&refValue);
}

/*----------------------------------------------------------------------------*
    Test if all discrete values are equal
*----------------------------------------------------------------------------*/
int ARM_ReferenceValue::IsConstant() const
{
	int		i(0), maxSize;
	double	epsilon = 1.e-6;
	double	testValue;

	// First test on itsCalculationMethod

	if ( itsCalculationMethod == K_CONSTANT_REF )
	{
		return 1;
	}

	// Second test on discrete values

	maxSize	=	(itsDiscreteValues == NULL ? 0 : itsDiscreteValues->GetSize());
	
	if ( maxSize > 0 )
	{
		testValue	=	(*itsDiscreteValues)[0];
	}

	for (i=0;i<maxSize;i++)
	{
		if ( (fabs( (testValue - (*itsDiscreteValues )[i]) ) > epsilon)	)
		{
			return 0;
		}
	}

	maxSize	=	(itsDiscreteValues2 == NULL ? 0 : itsDiscreteValues2->GetSize());

	if ( maxSize > 0 )
	{
		testValue	=	(*itsDiscreteValues2)[0];
	}

	for (i=0;i<maxSize;i++)
	{
		if ( (fabs( (testValue - (*itsDiscreteValues2)[i]) ) > epsilon)	)
		{
			return 0;
		}
	}

	// To return 1, we must have itsDiscreteValues [0] = ... = itsDiscreteValues [maxSize-1]
	// and (if it exists)        itsDiscreteValues2[0] = ... = itsDiscreteValues2[maxSize-1]
	// The calculation method has to be equal to K_CONSTANT_REF

	return 1;
}

/*----------------------------------------------------------------------------*
    Compute the reference value at specified date

    When obtained for Maximum, Minimum or Averaging (path dependent),
    the reference value must to be calculated using the pricing model
    within the security payoff calculation method 

*----------------------------------------------------------------------------*/ 

/** 
double ARM_ReferenceValue::CptReferenceValue(const ARM_Date& date) const
{
    double value;


    // calculation method for linear interpolation

    switch(itsCalculationMethod)
    {
        case K_CONSTANT_REF:
        {
            value = itsDiscreteValues->Elt(0);

            return(value);
        };
        break;

        // calculation method for linear interpolation
        case K_LININTERPOL_REF :
        {
            if (itsDiscreteDates)
            {
               if (itsExtrapolMeth)
               {
                  value = linInterpol2(itsDiscreteDates, date.GetJulian(), 
                                       itsDiscreteValues);
               }
               else
               {
                  value = linInterpol(itsDiscreteDates, date.GetJulian(), 
                                      itsDiscreteValues);
               }

               return(value);
            }
            else if (itsDiscreteValues)
            {
               value = itsDiscreteValues->Elt(0);

               return(value);
            }
            else
            {
               throw Exception(__LINE__, __FILE__, ERR_PRICING_PB,
               "Wrong Computation Method for Reference Value Object");
            }
        };
        break;

        case K_PERFECT_DISCRETE_REF:
        {
            value = 0.;

            if (itsDiscreteDates)
            {
                int size = itsDiscreteDates->GetSize();

                for (int ind = 0; ind < size; ind++)
                {
                    if (itsDiscreteDates->Elt(ind) == date.GetJulian())
                    {
                        value = (*itsDiscreteValues)[ind];
                        break;
                    }
                }
            }
            else
                value = 0.;

            return value;
        };
        break;

        case K_DISCRETE_REF:
        {
            if (itsDiscreteDates)
            {
               int indx = locateIndex(itsDiscreteDates, date.GetJulian());

               if ( indx == -1 )
                  indx = 0;

               value = itsDiscreteValues->Elt(indx);

               return(value);
            }
            else if (itsDiscreteValues)
            {
               value = itsDiscreteValues->Elt(0);

               return(value);
            }
            else
            {
               throw Exception(__LINE__, __FILE__, ERR_PRICING_PB,
                    "Wrong Computation Method for Reference Value Object");
            }
        };
        break;

        case K_STEPUP_RIGHT :
        {
            if (itsDiscreteDates)
            {
               int indx = locateIndex(itsDiscreteDates, date.GetJulian());

               int newIndex;


               if ( indx == -1 )
               {
                  newIndex = 0;
               }
               else
               {
                  int sz = itsDiscreteValues->GetSize();

                  if ( indx >= sz-1 )
                  {
                     newIndex = sz-1;
                  }
                  else 
                  {
                     if (fabs(itsDiscreteDates->Elt(indx)
                         -date.GetJulian()) <= 1e-6 )
                     {
                        newIndex = indx;
                     }
                     else
                     {
                        newIndex = indx+1;
                     }
                  }
               }

               value = itsDiscreteValues->Elt(newIndex);

               return(value);
            }
            else if (itsDiscreteValues)
            {
               value = itsDiscreteValues->Elt(0);

               return(value);
            }
            else
            {
               throw Exception(__LINE__, __FILE__, ERR_PRICING_PB,
                   "Wrong Computation Method for Reference Value Object");
            }
        };
        break;
 

        case K_STEPUP_LEFT :
        {
            if (itsDiscreteDates)
            {
               int indx = locateIndex(itsDiscreteDates, date.GetJulian());

               if ( indx == -1 )
                  indx = 0;

               int sz = itsDiscreteValues->GetSize();

               int newIndex;  

               if ( indx <= 0 )
               {
                  newIndex = 0;
               }
               else
               {
                  newIndex = indx ; // TMP -1 
               }

               value = itsDiscreteValues->Elt(newIndex);

               return(value);
            }
            else if (itsDiscreteValues)
            {
               value = itsDiscreteValues->Elt(0);

               return(value);
            }
            else
            {
               throw Exception(__LINE__, __FILE__, ERR_PRICING_PB,
                   "Wrong Computation Method for Reference Value Object");
            }
        };
        break;

    }
    
    return(0.0);
}
**/ 

/* nouvelle methode plus rapide */
double ARM_ReferenceValue::CptReferenceValue(double Juliandate) const
{
    //double value;
    // calculation method for linear interpolation

  //  switch(itsCalculationMethod)
  //  {
  //      case K_CONSTANT_REF:
  //      {
  //          value = (*itsDiscreteValues)[0];

  //          return(value);
  //      };
  //      break;

  //      // calculation method for linear interpolation
  //      case K_LININTERPOL_REF:
  //      {
  //          if (itsDiscreteDates)
  //          {   
  //             if (itsExtrapolMeth)
  //             {
  //                value = linInterpol2(itsDiscreteDates, 
  //                                     Juliandate, 
  //                                     itsDiscreteValues);
  //             }
  //             else
  //             {
  //                value = linInterpol(itsDiscreteDates->GetElt(), 
  //                                    itsDiscreteDates->GetSize(),
  //                                    Juliandate, 
  //                                    itsDiscreteValues->GetElt());
  //             }

  //             return(value);
  //          }
  //          else if (itsDiscreteValues)
  //          {
  //             value = itsDiscreteValues->Elt(0);

  //             return(value);
  //          }
  //          else
  //          {
  //             throw() /*Exception(__LINE__, __FILE__, ERR_PRICING_PB,
  //                 "Wrong Computation Method for Reference Value Object")*/;
  //          }
  //      };
  //      break;

  //      case K_PERFECT_DISCRETE_REF:
  //      {
  //          value = 0.;

  //          if (itsDiscreteDates)
  //          {
  //              int size = itsDiscreteDates->GetSize();

  //              for (int ind = 0; ind < size; ind++)
  //              {
  //                  if (itsDiscreteDates->Elt(ind) == Juliandate)
  //                  {
  //                      value = (*itsDiscreteValues)[ind];
  //                      break;
  //                  }
  //              }
  //          }
  //          else
  //              value = 0.;

  //          return value;
  //      };    
  //      break;

  //      case K_DISCRETE_REF:
  //      {
  //          if (itsDiscreteDates)
  //          {
  //             int indx = locateIndex(itsDiscreteDates, Juliandate);

  //             if ( indx == -1 )
  //                indx = 0;

  //             value = (*itsDiscreteValues)[indx];

  //             return(value);
  //          }
  //          else if (itsDiscreteValues)
  //          {
  //             value = (*itsDiscreteValues)[0];

  //             return(value);
  //          }
  //          else
  //          {
  //             throw() /*Exception(__LINE__, __FILE__, ERR_PRICING_PB,
  //                  "Wrong Computation Method for Reference Value Object")*/;
  //          }
  //      };
  //      break;

  //      case K_STEPUP_RIGHT :
  //      {
  //          if (itsDiscreteDates)
  //          {
  //             int indx = locateIndex(itsDiscreteDates, Juliandate);

  //             int newIndex;

  //             if ( indx == -1 )
  //             {
  //                newIndex = 0;
  //             }
  //             else
  //             {
  //                int sz = itsDiscreteValues->GetSize();

  //                if ( indx >= sz-1 )
  //                {
  //                   newIndex = sz-1;
  //                }
  //                else
  //                {
  //                   if (fabs(itsDiscreteDates->Elt(indx)
  //                       -Juliandate) <= 1e-6 )
  //                   {
  //                      newIndex = indx;
  //                   }
  //                   else
  //                   {
  //                      newIndex = indx+1;
  //                   }
  //                }
  //             }

  //             value = itsDiscreteValues->Elt(newIndex);

  //             return(value);
  //          }
  //          else if (itsDiscreteValues)
  //          {
  //             value = itsDiscreteValues->Elt(0);

  //             return(value);
  //          }
  //          else
  //          {
  //             throw() /*Exception(__LINE__, __FILE__, ERR_PRICING_PB,
  //                 "Wrong Computation Method for Reference Value Object")*/;
  //          }
  //      };
  //      break;

  //      case K_STEPUP_LEFT :
  //      {
  //          if (itsDiscreteDates)
  //          {
  //             int indx = locateIndex(itsDiscreteDates, Juliandate);

  //             if ( indx == -1 )
  //                indx = 0;

  //             int sz = itsDiscreteValues->GetSize();

  //             int newIndex;

  //             if ( indx <= 0 )
  //             {
  //                newIndex = 0;
  //             }
  //             else
  //             {
  //                newIndex = indx /* TMP -1 */;
  //             }

  //             value = itsDiscreteValues->Elt(newIndex);

  //             return(value);
  //          }
  //          else if (itsDiscreteValues)
  //          {
  //             value = itsDiscreteValues->Elt(0);

  //             return(value);
  //          }
  //          else
  //          {
  //             throw() /*Exception(__LINE__, __FILE__, ERR_PRICING_PB,
  //                 "Wrong Computation Method for Reference Value Object")*/;
  //          }
  //      };
  //      break;
		//
		//case K_LIN_INTER_CST_EXTRA_REF:
		//{
		//	if (itsDiscreteDates)
  //          {
		//		value = linInterpol2(itsDiscreteDates,
		//			                 Juliandate,
		//							 GetDiscreteValues());

  //              return(value);
  //           }
  //           else
  //           {
  //              throw() /*Exception(__LINE__, __FILE__, ERR_PRICING_PB,
  //                   "You trying to interpolate without reference dates")*/;
  //           }
  //       };
		//break;
  //  }
    return(0.0);
}


/**
ARM_Vector* ARM_ReferenceValue::CptReferenceValues(ARM_Date& date)
{
    double value = CptReferenceValue(date);
    
    ARM_Vector* values = new ARM_Vector(1, value);
    
    return(values);
}
**/ 



ARM_ReferenceValue* ARM_ReferenceValue::CptReferenceValues(ARM_Vector* date)
{
    ARM_ReferenceValue* newRefValue = (ARM_ReferenceValue*) this->Clone();
    ARM_Vector* newValues = new ARM_Vector(date->GetSize());
    ARM_Vector* newValues2 = NULL;


    if (itsDiscreteValues2)
       newValues2 = new ARM_Vector(date->GetSize());

    for (int i = 0; i < date->GetSize(); i++)
    {
       ARM_Vector* value = InterpolateMulti(date->Elt(i));

        newValues->Elt(i) = value->Elt(0);

        if (itsDiscreteValues2)
           newValues2->Elt(i) = value->Elt(1);

        delete value;
    }

    newRefValue->SetDiscreteDates(date);
    newRefValue->SetDiscreteValues(newValues, 0);
    newRefValue->SetDiscreteValues(newValues2, 1);

    return(newRefValue);
}


void ARM_ReferenceValue::InterpolateRefValue(ARM_Vector* date)
{
    ARM_ReferenceValue* newRefValue = (ARM_ReferenceValue*) this->Clone();
    size_t size1 = date->GetSize();
    ARM_Vector* newValues = new ARM_Vector(size1);

    ARM_Vector* newValues2 = NULL;
    if (itsDiscreteValues2)
       newValues2 = new ARM_Vector(size1);

    for (int i = 0; i < size1; ++i)
    {
       ARM_Vector* value = InterpolateMulti(date->Elt(i));

        newValues->Elt(i) = value->Elt(0);

        if (itsDiscreteValues2)
           newValues2->Elt(i) = value->Elt(1);

        delete value;
    }


    delete itsDiscreteDates;
    itsDiscreteDates = (ARM_Vector*)date->Clone();

    delete itsDiscreteValues;
    itsDiscreteValues = newValues;    
    
    delete itsDiscreteValues2;
    itsDiscreteValues2 = newValues2;

}


double ARM_ReferenceValue::ComputePrice(ARM_Model* model, ARM_Date* dateValo)
{
	double price (0.0);

	if ( itsDiscreteDates == NULL ) 
       return price;

	double asOf;

	if (dateValo)
		asOf = (*dateValo).GetJulian();
	/*else
		asOf = model->GetZeroCurve()->GetAsOfDate().GetJulian();*/

	for (int i = 0; i < itsDiscreteDates->GetSize(); i++)
	{
		double yt = (itsDiscreteDates->Elt(i) - asOf) / 365.;
		price += 0;//itsDiscreteValues->Elt(i) * model->ZeroPrice(0.0,yt);
	}

	return price;
}


void ARM_ReferenceValue::View(char* id, FILE* ficOut)
{
    /*FILE* fOut;
    char fOutName[200];

 
 
    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
 
       (void) unlink(fOutName);

       fOut = fopen(fOutName, "w");
    }
    else
    {
       fOut = ficOut;
    }

    fprintf(fOut,"\n");

    switch(itsCalculationMethod)
    {
        case K_CONSTANT_REF:
        {
			fprintf(fOut, "        *****> Interpolation Mode: K_CONSTANT_REF \n\n"); 
            fprintf(fOut, "Constant  : %14.4lf\n",(*itsDiscreteValues)[0]);

            if (itsDiscreteValues2)
               fprintf(fOut,"%14.4lf",(*itsDiscreteValues2)[0]);

            fprintf(fOut,"\n\n");
        };
        break;

		case K_LIN_INTER_CST_EXTRA_REF :
		case K_PERFECT_DISCRETE_REF    :
		case K_STEPUP_RIGHT			   :
		case K_STEPUP_LEFT			   :
        case K_LININTERPOL_REF         :
        {
			switch(itsCalculationMethod)
			{
				case K_LIN_INTER_CST_EXTRA_REF   :
				{
        			fprintf(fOut, "        *****> Interpolation Mode: K_LIN_INTER_CST_EXTRA_REF \n\n"); 
				};
				break;
        
				case K_PERFECT_DISCRETE_REF    :
				{
					fprintf(fOut, "        *****> Interpolation Mode: K_PERFECT_DISCRETE_REF \n\n"); 
				};
				break;

				case K_STEPUP_RIGHT   :
				{
        			fprintf(fOut, "        *****> Interpolation Mode: K_STEPUP_RIGHT \n\n"); 
				};
				break;
        
				case K_STEPUP_LEFT    :
				{
					fprintf(fOut, "        *****> Interpolation Mode: K_STEPUP_LEFT \n\n"); 
				};
				break;
				
				case K_LININTERPOL_REF:
				{
					fprintf(fOut, "        *****> Interpolation Mode: K_LININTERPOL_REF \n\n"); 
				};
			}

            if (itsDiscreteDates)
            {
               fprintf(fOut, "%14s\t%14s\t", "Date", "Value1");

               if (itsDiscreteValues2)
                  fprintf(fOut, "Value2\n\n");
               else
                  fprintf(fOut, "\n\n");

               char date[20];

               for (int i = 0; i < itsDiscreteDates->GetSize(); i++)
               {
                   try
				   {
					   ((ARM_Date) (*itsDiscreteDates)[i]).JulianToStrDateDay(date);
					   fprintf(fOut," %14s\t%14.4lf\t", date,(*itsDiscreteValues)[i]);
				   }
				   catch(...)
				   {
					   fprintf(fOut," %14.4lf\t%14.4lf\t", (*itsDiscreteDates)[i], (*itsDiscreteValues)[i]);
				   }
                   if (itsDiscreteValues2)
                      fprintf(fOut,"%14.4lf",(*itsDiscreteValues2)[i]);
                   fprintf(fOut,"\n");
               }
            }
			else
			{
				fprintf(fOut," %14.4lf\t", (*itsDiscreteValues)[0]);
			}
        };
        break;

        case K_DISCRETE_REF:
        {
			fprintf(fOut, "        *****> Interpolation Mode: DISCRETE_REF \n"); 
        
            if (itsDiscreteDates)
            {
               fprintf(fOut, "Discrete Values     : \n");
               fprintf(fOut, "%14s\t%14s\t", "Date", "Value1");

               if (itsDiscreteValues2)
                  fprintf(fOut, "Value2\n\n");
               else
                  fprintf(fOut, "\n\n");

               char date[20];
 
               for (int i = 0; i < itsDiscreteDates->GetSize(); i++)
               {
                   try
				   {
					   ((ARM_Date) (*itsDiscreteDates)[i]).JulianToStrDateDay(date);
					   fprintf(fOut," %14s\t%14.4lf\t", date,(*itsDiscreteValues)[i]);
				   }
				   catch(...)
				   {
					   fprintf(fOut," %14.4lf\t%14.4lf\t", (*itsDiscreteDates)[i], (*itsDiscreteValues)[i]);
				   }
                   if(itsDiscreteValues2)
                       fprintf(fOut,"%14.4lf",(*itsDiscreteValues2)[i]);

                   fprintf(fOut,"\n");
               }
            }
            else if (itsDiscreteValues)
            {
               fprintf(fOut, "Constant  : %14.4lf\t",(*itsDiscreteValues)[0]);
               if(itsDiscreteValues2)
                       fprintf(fOut,"%14.4lf",(*itsDiscreteValues2)[0]);
                fprintf(fOut,"\n");
            }
        };
        break;
    }

    fprintf(fOut,"\n\n");

    if ( ficOut == NULL )
    {
       fclose(fOut);
    }*/
}



void ARM_ReferenceValue::DisplayScheduleDates(int datesType, int viewInitialExch,
                                       char* id, FILE* ficOut)
{
 //   FILE* fOut;
 //   char fOutName[200];

 //   ARM_Vector* julDates;
 //   int vectSize;

 //   if ( ficOut == NULL )
 //   {
 //      ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);

 //      (void) unlink(fOutName);

 //      fOut = fopen(fOutName, "w");
 //   }
 //   else
 //   {
 //      fOut = ficOut;
 //   }

 //   julDates = itsDiscreteDates;        

 //   if ( julDates == NULL )
 //   {
 //      vectSize = 0;
 //   }
 //   else
 //   {
 //      vectSize = julDates->GetSize();
 //   }

 //   // Now Print out Vect values in date format

 //   char d1[20];
	//fprintf(fOut,"%d\n", vectSize);

 //   for (int i = 0;  i< vectSize; i++)
 //   {
 //      ((ARM_Date) (julDates->Elt(i))).JulianToStrDate(d1);
 //      fprintf(fOut,"%s\n", d1);
 //   }

 //   if ( ficOut == NULL )
 //   {
 //      fclose(fOut);
 //   }
}



void ARM_ReferenceValue::DisplayScheduleValues(int valuesType,
                                        char* id, FILE* ficOut)
{
    //FILE* fOut;
    //char fOutName[200];

    //ARM_Vector* vectValues;
    //int vectSize;



    //if ( ficOut == NULL )
    //{
    //   ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);

    //   (void) unlink(fOutName);

    //   fOut = fopen(fOutName, "w");
    //}
    //else
    //{
    //   fOut = ficOut;
    //}


    //vectValues = itsDiscreteValues;
    //    
    //if ( vectValues == NULL )
    //{
    //   vectSize = 0;
    //}
    //else
    //{
    //   vectSize = vectValues->GetSize();
    //}

    //// Now Print out Vect values

    //fprintf(fOut,"%d\n", vectSize);

    //for (int i = 0; i < vectSize; i++)
    //{
    //       fprintf(fOut,"%lf\n", vectValues->Elt(i));
    //}

    //if ( ficOut == NULL )
    //{
    //   fclose(fOut);
    //}
}


/*******************************************************************/
/*                                                                 */
/*  ReferenceValue algebra                                         */
/*  Only possible for the following computation                    */
/*  CONST    - CONST                                               */
/*  CONST    - DISCRETE                                            */
/*  DISCRETE - CONST                                               */
/*  DISCRETE - DISCRETE                                            */
/*                                                                 */
/*******************************************************************/


//ARM_ReferenceValue ARM_ReferenceValue::operator + (ARM_ReferenceValue& refVal2)
//{
//    //if (( this->GetCalcMethod() == K_CONSTANT_REF ) 
//    //    || 
//    //    ( refVal2.GetCalcMethod() == K_CONSTANT_REF )
//    //   )
//    //{
//    //    ARM_ReferenceValue *constRef, *otherRef;
//
//    //    // works whatever the type of otherRef
//    //    if ( this->GetCalcMethod() == K_CONSTANT_REF )
//    //    {
//    //       constRef = this;
//
//    //       otherRef = &refVal2;
//    //    }
//    //    else
//    //    {
//    //       constRef = &refVal2;
//
//    //       otherRef = this;
//    //    }
//
//    //    double value = constRef->GetDiscreteValues()->Elt(0);
//
//    //    ARM_ReferenceValue newRef(*otherRef);
//
//    //    int size  = otherRef->GetDiscreteValues()->GetSize();
//
//    //    double* valuesElt = newRef.GetDiscreteValues()->GetElt();
//
//    //    for (int i = size-1; i>=0; i--)
//    //    {
//    //        valuesElt[i] += value;
//    //    }
//
//    //    return newRef;
//    //}
//    //else
//    //{
//    //   if (( this->GetCalcMethod() != K_DISCRETE_REF )
//    //       || 
//    //       ( refVal2.GetCalcMethod() != K_DISCRETE_REF )
//    //      )
//    //   {
//    //      throw() /*Exception(__LINE__, __FILE__, ERR_PRICING_PB,
//    //           "Cannot add refvalues other than stepup")*/;
//    //   }
//
//    //   ARM_Vector* mergedDates = NULL;
//
//    //   if (refVal2.GetDiscreteDates())
//    //   {
//    //      if (this->GetDiscreteDates())
//    //         MergeDates(&mergedDates,
//    //                 this->GetDiscreteDates(),
//    //                 refVal2.GetDiscreteDates());
//    //      else 
//    //         mergedDates = (ARM_Vector*) refVal2.GetDiscreteDates()->Clone();
//    //   }
//    //   else
//    //   {
//    //      if (this->GetDiscreteDates())
//    //         mergedDates = (ARM_Vector*) this->GetDiscreteDates()->Clone();
//    //   }
//
//    //   if (mergedDates)
//    //   {
//    //      int i, size = mergedDates->GetSize();
//
//    //      ARM_Vector* mergedValues = new ARM_Vector(size);
//
//    //      for (i = 0; i < size; i++)
//    //      {
//    //          mergedValues->Elt(i) = this->CptReferenceValue(mergedDates->Elt(i))
//    //                         +refVal2.CptReferenceValue(mergedDates->Elt(i));
//    //      }
//
//    //      ARM_ReferenceValue theResult(mergedDates, mergedValues);
//
//    //      theResult.SetCalcMethod(K_DISCRETE_REF);
//
//    //      return(theResult);            
//    //   }
//    //   else
//    //   {
//    //      double value =  this->CptReferenceValue(0)
//    //                              +refVal2.CptReferenceValue(0);
//
//    //      ARM_ReferenceValue theResult(value);
//
//    //      theResult.SetCalcMethod(K_DISCRETE_REF);
//
//    //      return(theResult);            
//    //   }
//    //}
//}
//        
//
//ARM_ReferenceValue ARM_ReferenceValue::operator * 
//                                      (ARM_ReferenceValue& refVal2)
//{
//    if (( this->GetCalcMethod() == K_CONSTANT_REF ) 
//        || 
//        ( refVal2.GetCalcMethod() == K_CONSTANT_REF )
//       ) 
//    {
//        ARM_ReferenceValue *constRef, *otherRef;
//        
//        // works whatever the type of otherRef
//        if ( this->GetCalcMethod() == K_CONSTANT_REF )
//        {
//           constRef = this;
//
//           otherRef = &refVal2;
//        }
//        else
//        {
//           constRef = &refVal2;
//
//           otherRef = this;
//        }
//
//        double value = constRef->GetDiscreteValues()->Elt(0);
//
//        ARM_ReferenceValue newRef(*otherRef);
//
//        int size  = otherRef->GetDiscreteValues()->GetSize();
//    
//        double* valuesElt = newRef.GetDiscreteValues()->GetElt();
//
//        for (int i = size-1; i>=0; i--)
//        {
//            valuesElt[i] *= value;
//        }
//
//        return newRef;
//    }
//    else
//    {
//       if (( this->GetCalcMethod() != K_DISCRETE_REF )
//           ||   
//           ( refVal2.GetCalcMethod() != K_DISCRETE_REF )
//          )
//       {
//          throw() /*Exception(__LINE__, __FILE__, ERR_PRICING_PB,
//             "Cannot add refvalues other than stepup")*/;
//        
//       }
//
//       ARM_Vector* mergedDates = NULL;
//
//       if (refVal2.GetDiscreteDates())
//       {
//          if (this->GetDiscreteDates())
//             MergeDates(&mergedDates, 
//                        this->GetDiscreteDates(), 
//                        refVal2.GetDiscreteDates());
//          else 
//             mergedDates = (ARM_Vector*) refVal2.GetDiscreteDates()->Clone();
//       }
//       else
//       {
//          if (this->GetDiscreteDates())
//             mergedDates = (ARM_Vector*) this->GetDiscreteDates()->Clone();
//       }
//
//       if (mergedDates)
//       {
//          int i, size = mergedDates->GetSize();
//
//          ARM_Vector* mergedValues = new ARM_Vector(size);
//
//          for (i = 0; i < size; i++)
//          {
//              mergedValues->Elt(i) = 
//                    this->CptReferenceValue(mergedDates->Elt(i))
//                        * refVal2.CptReferenceValue(mergedDates->Elt(i));
//          }
//
//          ARM_ReferenceValue theResult(mergedDates, mergedValues);
//  
//          theResult.SetCalcMethod(K_DISCRETE_REF);
//
//          return(theResult);
//       }
//       else
//       {
//          double value =  this->CptReferenceValue(0)
//                                   *refVal2.CptReferenceValue(0);
//
//          ARM_ReferenceValue theResult(value);
//
//          theResult.SetCalcMethod(K_DISCRETE_REF);
//
//          return(theResult);
//       }
//    }
//}
//
//ARM_ReferenceValue ARM_ReferenceValue::operator / 
//                                      (ARM_ReferenceValue& refVal2)
//{
//    if (( this->GetCalcMethod() == K_CONSTANT_REF ) 
//        || 
//        ( refVal2.GetCalcMethod() == K_CONSTANT_REF )
//       ) 
//    {
//        ARM_ReferenceValue *constRef, *otherRef;
//        
//        // works whatever the type of otherRef
//        if ( this->GetCalcMethod() == K_CONSTANT_REF )
//        {
//           constRef = this;
//
//           otherRef = &refVal2;
//        }
//        else
//        {
//           constRef = &refVal2;
//
//           otherRef = this;
//        }
//
//        double value = constRef->GetDiscreteValues()->Elt(0);
//
//        ARM_ReferenceValue newRef(*otherRef);
//
//        int size  = otherRef->GetDiscreteValues()->GetSize();
//    
//        double* valuesElt = newRef.GetDiscreteValues()->GetElt();
//
//        for (int i = size-1; i>=0; i--)
//        {
//            valuesElt[i] /= value;
//        }
//
//        return newRef;
//    }
//    else
//    {
//       if (( this->GetCalcMethod() != K_DISCRETE_REF )
//           ||   
//           ( refVal2.GetCalcMethod() != K_DISCRETE_REF )
//          )
//       {
//          throw() /*Exception(__LINE__, __FILE__, ERR_PRICING_PB,
//             "Cannot add refvalues other than stepup")*/;
//        
//       }
//
//       ARM_Vector* mergedDates = NULL;
//
//       if (refVal2.GetDiscreteDates())
//       {
//          if (this->GetDiscreteDates())
//             MergeDates(&mergedDates, 
//                        this->GetDiscreteDates(), 
//                        refVal2.GetDiscreteDates());
//          else 
//             mergedDates = (ARM_Vector*) refVal2.GetDiscreteDates()->Clone();
//       }
//       else
//       {
//          if (this->GetDiscreteDates())
//             mergedDates = (ARM_Vector*) this->GetDiscreteDates()->Clone();
//       }
//
//       if (mergedDates)
//       {
//          int i, size = mergedDates->GetSize();
//
//          ARM_Vector* mergedValues = new ARM_Vector(size);
//
//          for (i = 0; i < size; i++)
//          {
//              mergedValues->Elt(i) = 
//                    this->CptReferenceValue(mergedDates->Elt(i))
//                        / refVal2.CptReferenceValue(mergedDates->Elt(i));
//          }
//
//          ARM_ReferenceValue theResult(mergedDates, mergedValues);
//  
//          theResult.SetCalcMethod(K_DISCRETE_REF);
//
//          return(theResult);
//       }
//       else
//       {
//          double value =  this->CptReferenceValue(0)
//                                   *refVal2.CptReferenceValue(0);
//
//          ARM_ReferenceValue theResult(value);
//
//          theResult.SetCalcMethod(K_DISCRETE_REF);
//
//          return(theResult);
//       }
//    }
//}
//


/*------------------------------------------------------------------------------*
    Compute the reference value at specified date
    When obtained for Maximum, Minimum or Averaging (path dependent),
    the reference va lue must to be calculated
    using the pricing model within the security payoff calculation method
*-------------------------------------------------------------------------------*/
void ARM_ReferenceValue::SetXYValue(int index, double X,double Y, int dim)
{
    if ( itsDiscreteDates->GetSize() != GetDiscreteValues(dim)->GetSize() )
       return; //throw() /*Exception(__LINE__, __FILE__,0,"Incompatible!")*/;

    itsDiscreteDates->Elt(index)= X;

    GetDiscreteValues(dim)->Elt(index) = Y;
}



void ARM_ReferenceValue::SetXYValue(int index, double X,ARM_Vector* Y)
{
    if (( itsDiscreteDates->GetSize() != GetDiscreteValues(0)->GetSize()) 
        ||
        ( itsDiscreteDates->GetSize() != GetDiscreteValues(1)->GetSize()))
       return;// throw() /*Exception(__LINE__, __FILE__,0,"Incompatible!")*/;

    itsDiscreteDates->Elt(index)= X;
    GetDiscreteValues(0)->Elt(index) = Y->Elt(0);
    GetDiscreteValues(1)->Elt(index) = Y->Elt(1);
}




ARM_ReferenceValue* ARM_ReferenceValue::CreateStartResetAdjust(int spot,
                                                               char* ccy)
{
    ARM_ReferenceValue* newRef = (ARM_ReferenceValue*) this->Clone();

    ARM_Vector* DiscDates = newRef->GetDiscreteDates();


    if (!DiscDates)
       return newRef;

    double* newRefDates = DiscDates->GetElt(); 

    int i, size = itsDiscreteDates->GetSize();    

    for (i = 0; i < size; i++)
    {
        ARM_Date aDate(newRefDates[i]);

        newRefDates[i] = (aDate.PreviousBusinessDay(spot, ccy)).GetJulian();
    }

    return(newRef);
}



double ARM_ReferenceValue::Interpolate(ARM_Date& date, int dim)
{
    double value = Interpolate(date.GetJulian(), dim);

    return(value);
}



double ARM_ReferenceValue::Interpolate(double Juliandate, int dim)
{
    double value;

 //   // calculation method for linear interpolation
 //   switch(itsCalculationMethod)
 //   {
 //        case K_CONSTANT_REF:
 //        {
 //            value = (*GetDiscreteValues(dim))[0];

 //            return(value);
 //        };
 //        break;

 //        // calculation method for linear interpolation
 //        case K_LININTERPOL_REF:
 //        {
 //            if (itsDiscreteDates)
 //            {
 //               if (itsExtrapolMeth)
 //               {
 //                  value = linInterpol2(itsDiscreteDates,
 //                                       Juliandate, 
 //                                       GetDiscreteValues(dim));
 //               }
 //               else
 //               {
 //                  value = linInterpol(itsDiscreteDates->GetElt(),
 //                                      itsDiscreteDates->GetSize(),
 //                                      Juliandate, 
 //                                      GetDiscreteValues(dim)->GetElt());
 //               }

 //               return(value);
 //            }
 //            else
 //            {
 //               throw() /*Exception(__LINE__, __FILE__, ERR_PRICING_PB,
 //                    "Wrong Computat ion Method for Reference Value Object")*/;
 //            }
 //        };

 //        case K_LIN_INTER_CST_EXTRA_REF:
 //        {
 //            if (itsDiscreteDates)
 //            {
 //                  value = linInterpol2(itsDiscreteDates,
 //                                       Juliandate, 
 //                                       GetDiscreteValues(dim));
 //               return(value);
 //            }
 //            else
 //            {
 //               throw() /*Exception(__LINE__, __FILE__, ERR_PRICING_PB,
 //                    "You trying to interpolate without reference dates")*/;
 //            }
 //        };
 //        break;

 //        case K_DISCRETE_REF:
 //        {
 //            if (itsDiscreteDates)
 //            {
 //               int indx = locateIndex(itsDiscreteDates, Juliandate);

 //               if ( indx == -1 )
 //                  indx = 0;

 //               value = (*GetDiscreteValues(dim))[indx];

 //               return(value);
 //            }
 //            else if (GetDiscreteValues(dim))
 //            {
 //               value = (*GetDiscreteValues(dim))[0];

 //               return(value);
 //            }
 //            else
 //            {
 //               throw() /*Exception(__LINE__, __FILE__, ERR_PRICING_PB,
 //                     "Wrong Computat ion Method for Reference Value Object")*/;
 //            }
 //        };
 //        break;

 //        case K_STEPUP_RIGHT :
 //        {
 //            if (itsDiscreteDates)
 //            {
 //               int indx = locateIndex(itsDiscreteDates, Juliandate);

 //               int newIndex;

 //               if ( indx == -1 )
 //               {
 //                  newIndex = 0;
 //               }
 //               else
 //               {
 //                  int sz = GetDiscreteValues(dim)->GetSize();

 //                  if ( indx >= sz-1 )
 //                  {
 //                     newIndex = sz-1;
 //                  }
 //                  else
 //                  {
 //                     if (fabs(itsDiscreteDates->Elt(indx)
 //                         -Juliandate) <= 1e-6 )
 //                     {
 //                        newIndex = indx;
 //                     }
 //                     else
 //                     {
 //                        newIndex = indx+1;
 //                     }
 //                  }
 //              }

 //              value = GetDiscreteValues(dim)->Elt(newIndex);

 //              return(value);
 //            }
 //            else
 //               if (GetDiscreteValues(dim))
 //                  return  GetDiscreteValues(dim)->GetElt()[0];
 //               else
 //                  throw() /*Exception(__LINE__, __FILE__, ERR_PRICING_PB,
 //                 "Wrong Computation Method for Reference Value Object")*/;
 //        };
 //        break;


 //       case K_STEPUP_LEFT :
 //       {
 //           if (itsDiscreteDates)
 //           {
 //              int indx = locateIndex(itsDiscreteDates, Juliandate);

 //              if ( indx == -1 )
 //                 indx = 0;

 //              int sz = GetDiscreteValues(dim)->GetSize();

 //              int newIndex;

 //              if ( indx <= 0 )
 //              {
 //                 newIndex = 0;
 //              }
 //              else
 //              {
 //                 newIndex = indx /* TMP -1 */;
 //              }

 //              value = GetDiscreteValues(dim)->Elt(newIndex);

 //              return(value);
 //           }
	//		else
 //               if (GetDiscreteValues(dim))
 //                  return  GetDiscreteValues(dim)->GetElt()[0];
 //               else
 //                  throw() /*Exception(__LINE__, __FILE__, ERR_PRICING_PB,
 //                 "Wrong Computation Method for Reference Value Object")*/;
 //       };
 //       break;
	//}

    return(0.0);
}



ARM_Vector* ARM_ReferenceValue::InterpolateMulti(ARM_Date& date)
{
    ARM_Vector* X = InterpolateMulti(date.GetJulian());

    return (X);
}



ARM_Vector* ARM_ReferenceValue::InterpolateMulti(double Juliandate)
{
    // calculation method for linear interpolation
    //int size = 0;

    //if (itsDiscreteValues)
    //{
    //    size++;
    //}

    //if (itsDiscreteValues2)
    //{
    //    size++;
    //}

    //if ( size == 0 )
    //{
    //    throw() /*Exception(__LINE__, __FILE__, ERR_PRICING_PB,
    //                    "Reference Value Object expected")*/;
    //}

    //ARM_Vector* X = new ARM_Vector(size);

    //switch(itsCalculationMethod)
    //{
    //    case K_CONSTANT_REF:
    //    {
    //        if (itsDiscreteValues)
    //           X->Elt(0) = itsDiscreteValues->Elt(0);

    //        if (itsDiscreteValues2)
    //           X->Elt(1) = itsDiscreteValues2->Elt(0);

    //        return (X);
    //    };
    //    break;

    //    // calculation method for linear interpolation
    //    case K_LININTERPOL_REF:
    //    {
    //        if (itsDiscreteDates)
    //        {
    //           if (itsDiscreteValues)
    //           {
    //              if (itsExtrapolMeth)
    //              {
    //                 X->Elt(0) = linInterpol2(itsDiscreteDates,
    //                                          Juliandate,
    //                                          itsDiscreteValues);
    //                 if (itsDiscreteValues2)
    //                 {
    //                     X->Elt(1) = linInterpol2(itsDiscreteDates,
    //                                          Juliandate,
    //                                          itsDiscreteValues2);
    //                 }
    //              }
    //              else
    //              {
    //                 X->Elt(0) = linInterpol(itsDiscreteDates->GetElt(),
    //                                         itsDiscreteDates->GetSize(),
    //                                         Juliandate,
    //                                         itsDiscreteValues->GetElt());
    //                 if (itsDiscreteValues2)
    //                 {
    //                     X->Elt(1) = linInterpol(itsDiscreteDates->GetElt(),
    //                                         itsDiscreteDates->GetSize(),
    //                                         Juliandate,
    //                                         itsDiscreteValues2->GetElt());
    //                 }
    //              }
    //           }

    //           return (X);
    //        }
    //        else
    //        {
    //           throw() /*Exception(__LINE__, __FILE__, ERR_PRICING_PB,
    //               "Wrong Computat ion Method for Reference Value Object")*/;
    //            }
    //    };
    //    break;
    //    
    //    case K_DISCRETE_REF:
    //    {
    //        if (itsDiscreteDates)
    //        {
    //           if (itsDiscreteValues)
    //           {
    //               int indx = locateIndex(itsDiscreteDates, Juliandate);

    //               if ( indx == -1 )
    //                  indx = 0;

    //               X->Elt(0) = (*itsDiscreteValues)[indx];
    //               if (itsDiscreteValues2)
    //               {
    //                   X->Elt(0) = (*itsDiscreteValues2)[indx];
    //               }
    //               
    //            }

    //            return(X);
    //        }
    //        
    //        else
    //        {
    //           throw() /*Exception(__LINE__, __FILE__, ERR_PRICING_PB,
    //                "Wrong Computation Method for Reference Value Object")*/;
    //        }
    //    };
    //    break;

    //    case K_STEPUP_RIGHT:
    //    {
    //        if (itsDiscreteDates)
    //         {
    //            int indx = locateIndex(itsDiscreteDates, Juliandate);

    //            int newIndex;

    //            if ( indx == -1 )
    //            {
    //               newIndex = 0;
    //            }
    //            else
    //            {
    //               int sz = itsDiscreteDates->GetSize();

    //               if ( indx >= sz-1 )
    //               {
    //                  newIndex = sz-1;
    //               }
    //               else
    //               {
    //                  if (fabs(itsDiscreteDates->Elt(indx)
    //                      -Juliandate) <= 1e-6 )
    //                  {
    //                     newIndex = indx;
    //                  }
    //                  else
    //                  {
    //                     newIndex = indx+1;
    //                  }
    //               }
    //            }
    //            X->Elt(0) = itsDiscreteValues->Elt(newIndex);
    //            if (itsDiscreteValues2)
    //            {
    //                X->Elt(1) = itsDiscreteValues2->Elt(newIndex);
    //            }
    //            return(X);
    //        }
    //   
    //  };
    //  break;
    //    case K_STEPUP_LEFT:
    //    {
    //        if (itsDiscreteDates)
    //        {
    //            if (itsDiscreteValues)
    //            {
    //               int indx = locateIndex(itsDiscreteDates, Juliandate);

    //               if ( indx == -1 )
    //                  indx = 0;

    //               int sz = itsDiscreteValues->GetSize();

    //               int newIndex;

    //               if ( indx <= 0 )
    //               {
    //                  newIndex = 0;
    //               }
    //               else
    //               {
    //                  newIndex = indx;
    //               }

    //               X->Elt(0) = itsDiscreteValues->Elt(newIndex);
    //               if (itsDiscreteValues2)
    //               {
    //                   X->Elt(1) = itsDiscreteValues2->Elt(newIndex);
    //               }

    //               return(X);
    //            }
    //        }            
    //        else
    //        {
    //           throw() /*Exception(__LINE__, __FILE__, ERR_PRICING_PB,
    //               "Wrong Computation Method for Reference Value Object")*/;
    //        }
    //    };
    //    break;
    //}

    return(NULL);
}





/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
