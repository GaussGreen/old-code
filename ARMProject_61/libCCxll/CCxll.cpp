/*-------------------------------------------------------------------------*
 
   File: %M% 
   Path: %P% 
   Description: Interface de lecture/ecriture des donnees Excel (type LPXLOPER)
   Created: 00/01/03
   Author: Charles-Emmanuel MUSY
   Modified: %E% %U% 
   Last maintained by: Charles-Emmanuel MUSY
   Revision: %I% 
 
*--------------------------------------------------------------------------*/



#define CCxll_cc
#include "CCxll.h"
#include <limits.h>



int IsVBA()
{
	XLOPER xRes;

	Excel4(xlfCaller, &xRes, 0);

	if ( (xRes.xltype ==  xltypeErr) || (xRes.xltype ==  xltypeStr) )
		return 1;
	else
		return 0;
}

int XL_getNumCell(LPXLOPER xlTab, double* cell, double* defaultValue)
{
  	LPXLOPER px;
	XLOPER xMulti;
  	int error = XL_NO_ERROR;


	
   try
   {
     /// this is very fundamental as there is no guarantee
	 /// that the pointer is not null!
     // But this should never happen

     if ( xlTab != NULL )
	 {
		switch(xlTab->xltype) 
		{
  			case xltypeNil:
			case xltypeMissing:
            {
				if (defaultValue)
				{
				   *cell = *defaultValue;
				}
				else
				{
					error = xlerrValue;
				}
            };
			break;
			
            case xltypeBool:
				*cell = xlTab->val.boolean;
			break;
			
            case xltypeNum:
				*cell = xlTab->val.num;
			break;
			
            case xltypeRef:
			case xltypeSRef:
			case xltypeMulti:
            {
				if ( xlretUncalced == Excel (xlCoerce, &xMulti, 2, (LPXLOPER)xlTab, TempNum (xltypeMulti))) 
				{
					error = xlerrValue;
					return error;
				}
				
				if (xMulti.val.array.rows * xMulti.val.array.columns > 1)
				{
					error = xlerrValue;
					return error;
				}

                px = xMulti.val.array.lparray;
				
                switch (px->xltype) 
				{
					
                    case xltypeNum:
						*cell = px->val.num;
					break;
					
                    case xltypeErr:
						error = px->val.err;
					break;
					
                    case xltypeNil:
					case xltypeMissing:
                    {
						if (defaultValue)
						{
						   *cell = *defaultValue;
						}
						else
						{
						   error = xlerrValue;
						}
                    };
					break;
					
                    case xltypeBool:
						*cell = px->val.boolean;
					break;
					
                    default:
						error = xlerrValue;
					break;
				}
			
                Excel (xlFree, 0, 1, (LPXLOPER)&xMulti);
            };
            break;
			
            case xltypeErr:
				error = xlTab->val.err;
			break;
		
            default:
				error = xlerrValue;
			break;
  		}
     }
   }
   
   catch(...)
   {
	   error = xlerrValue;

	   return(error);
   }

   return(error);
}



int XL_getNumArray(LPXLOPER xlTab, double** array, long* size) 
{
	LPXLOPER px;
  	XLOPER xMulti;
  	WORD i;
  	int error = XL_NO_ERROR;
	double* temp = NULL;

	*size = 0;

   try
   {
	 /// this is very fundamental as there is no guarantee
	 /// that the pointer is not null!
	 if ( xlTab != NULL )
     {
		switch(xlTab->xltype) 
		{
  			case xltypeMissing:
				error = xlerrValue;
			break;
	
            case xltypeNum:
				error = xlerrValue;
			break;
			
            case xltypeRef:
			case xltypeSRef:
			case xltypeMulti:
            {
				if(xlretUncalced == Excel (xlCoerce, &xMulti, 2, (LPXLOPER)xlTab, TempNum (xltypeMulti))) 
				{
					error = xlerrValue;
					return error;
				}
			
                temp = (double*)malloc (xMulti.val.array.rows * xMulti.val.array.columns * sizeof (double)); 
			
                for(i = 0; i < (xMulti.val.array.rows * xMulti.val.array.columns); i++) 
				{
					px = xMulti.val.array.lparray + i;
					switch(px->xltype) 
					{
						case xltypeNum:
							temp[*size] = px->val.num;
							(*size)++;
						break;
						
                        case xltypeErr:
							error = px->val.err;
						break;
						
                        case xltypeNil:
							error = xlerrValue;
						break;
						
                        default:
							error = xlerrValue;
						break;
					}
				}

				Excel (xlFree, 0, 1, (LPXLOPER)&xMulti);
            };
            break;
			
            case xltypeErr:
				error = xlTab->val.err;
			break;
			
            default:
				error = xlerrValue;
			break;
		}
	
        if ( error != XL_NO_ERROR)
		{
			if(temp != NULL)
			{
				free (temp);
			}

			*array = NULL;
		}
		else
		{
			*array = temp;
		}
     }
   }
	
   catch(...)
   {
	   error = xlerrValue;

	   return(error);
   }
	
   return(error);
}



int XL_getNumVector(LPXLOPER xlTab, VECTOR<double>& vec, VECTOR<double>* vecDefault) 
{
  	LPXLOPER px;
  	XLOPER xMulti;
  	WORD i;
  	int error = XL_NO_ERROR;

	vec.clear ();
	

	/// this is very fundamental as there is no guarantee
	/// that the pointor is not null!

   try
   {
     if ( xlTab != NULL )
     {
		switch(xlTab->xltype) 
		{
  			case xltypeNil:
  			case xltypeMissing:
            {
				if (vecDefault)
				{
					for (int i=0;i<vecDefault->size();i++)
						vec.push_back((*vecDefault)[i]);
				}
				else
				{
					error = xlerrValue;
				}
            };
			break;

			case xltypeNum:
				error = xlerrValue;
			break;
			
            case xltypeRef:
			case xltypeSRef:
			case xltypeMulti:
            {
				if ( xlretUncalced == Excel (xlCoerce, &xMulti, 2, (LPXLOPER)xlTab, TempNum (xltypeMulti))) 
				{
					error = xlerrValue;
					return error;
				}

				for (i = 0; i < (xMulti.val.array.rows * xMulti.val.array.columns); i++) 
				{
					px = xMulti.val.array.lparray + i;
					switch(px->xltype) 
					{
						case xltypeNum:
							vec.push_back (px->val.num);
						break;
						
                        case xltypeErr:
							error = px->val.err;
						break;
						
                        case xltypeNil:
							if (vecDefault)
							{
								vec.resize(0);
								for (int i=0;i<vecDefault->size();i++)
									vec.push_back((*vecDefault)[i]);
							}
							else
							{
								error = xlerrValue;
							}
						break;
						
                        default:
							error = xlerrValue;
							break;
					}
				}
			
                Excel (xlFree, 0, 1, (LPXLOPER)&xMulti);
            };
            break;
			
            case xltypeErr:
				error = xlTab->val.err;
			break;

            default:
				error = xlerrValue;
			break;
		}
     }
   }

   catch(...)
   {
	   error = xlerrValue;

	   return(error);
   }

   return(error);
}

// Fonction permetttant de recuperer un vecteur comportant des trous. Remplit les cases vides avec la valeur CREDIT_DEFAULT_VALUE
int XL_getNumVectorWithHole(LPXLOPER xlTab, VECTOR<double>& vec, VECTOR<double>* vecDefault) 
{
  	LPXLOPER px;
  	XLOPER xMulti;
  	WORD i;
  	int error = XL_NO_ERROR;

	vec.clear ();
	

	/// this is very fundamental as there is no guarantee
	/// that the pointor is not null!

   try
   {
     if ( xlTab != NULL )
     {
		switch(xlTab->xltype) 
		{
  			case xltypeNil:
  			case xltypeMissing:
            {
				if (vecDefault)
				{
					for (int i=0;i<vecDefault->size();i++)
						vec.push_back((*vecDefault)[i]);
				}
				else
				{
					error = xlerrValue;
				}
            };
			break;

			case xltypeNum:
				error = xlerrValue;
			break;
			
            case xltypeRef:
			case xltypeSRef:
			case xltypeMulti:
            {
				if ( xlretUncalced == Excel (xlCoerce, &xMulti, 2, (LPXLOPER)xlTab, TempNum (xltypeMulti))) 
				{
					error = xlerrValue;
					return error;
				}

				for (i = 0; i < (xMulti.val.array.rows * xMulti.val.array.columns); i++) 
				{
					px = xMulti.val.array.lparray + i;
					switch(px->xltype) 
					{
						case xltypeNum:
							vec.push_back (px->val.num);
						break;
						
                        case xltypeErr:
							error = px->val.err;
						break;
						
                        case xltypeNil:
							vec.push_back (XL_TYPE_DEFAULT_VALUE);
						break;
						
                        default:
							error = xlerrValue;
							break;
					}
				}
			
                Excel (xlFree, 0, 1, (LPXLOPER)&xMulti);
            };
            break;
			
            case xltypeErr:
				error = xlTab->val.err;
			break;

            default:
				error = xlerrValue;
			break;
		}
     }
   }

   catch(...)
   {
	   error = xlerrValue;

	   return(error);
   }

   return(error);
}

int XL_getNumVectorAndSize(LPXLOPER xlTab, long& nbRows, long& nbColumns, VECTOR<double>& vec) 
{
  	LPXLOPER px;
  	XLOPER xMulti;
  	WORD i;
  	int error = XL_NO_ERROR;

	/// to make sure we return some result, we first initialise this to
	/// 0 and clear the vector of result!
	nbRows = nbColumns = 0;
	vec.clear();
	
	/// this is very fundamental as there is no guarantee
	/// that the pointor is not null!

   try
   {
     if ( xlTab != NULL )
     {
		switch(xlTab->xltype) 
		{
  			case xltypeNil:
  			case xltypeMissing:
				error = xlerrValue;
			break;
			
            case xltypeNum:
				nbRows = 1;
				nbColumns = 1;
				vec.push_back(xlTab->val.num);
			break;

			case xltypeRef:
			case xltypeSRef:
			case xltypeMulti:
            {
				if(xlretUncalced == Excel (xlCoerce, &xMulti, 2, (LPXLOPER)xlTab, TempNum (xltypeMulti))) 
				{
					error = xlerrValue;
					return error;
				}
			
                for(i = 0; i < (xMulti.val.array.rows * xMulti.val.array.columns); i++) 
				{
					px = xMulti.val.array.lparray + i;
				
                    switch(px->xltype) 
					{
						case xltypeNum:
							vec.push_back (px->val.num);
					    break;
						
                        case xltypeErr:
							error = px->val.err;
						break;
						
                        case xltypeNil:
							error = xlerrValue;
						break;
						
                        default:
							error = xlerrValue;
						break;
					}
				}
			
                nbRows = xMulti.val.array.rows;
				nbColumns = xMulti.val.array.columns;

				Excel (xlFree, 0, 1, (LPXLOPER)&xMulti);
            };
            break;
			
            case xltypeErr:
				error = xlTab->val.err;
			break;

            default:
				error = xlerrValue;
			break;
		}
     }
   }

   catch(...)
   {
	   error = xlerrValue;

	   return(error);
   }

   return(error); 
}



int XL_getNumVectorVector (LPXLOPER xlTab, VECTOR< VECTOR<double> >& vecvec) 
{
  	LPXLOPER px;
  	XLOPER xMulti;
  	WORD i,j,k,size;
  	int error = XL_NO_ERROR;

	/// to make sure we return some result, we first initialise this to
	/// 0 and clear the vector of result!
	vecvec.clear();
	

   try
   {
     /// this is very fundamental as there is no guaranty
	 /// that the pointer is not null!

     if ( xlTab != NULL )
     {
		switch(xlTab->xltype) 
		{
  			case xltypeNil:
  			case xltypeMissing:
				//error = xlerrValue;
			break;
			
            case xltypeNum:
				error = xlerrValue;
			break;
			
            case xltypeRef:
			case xltypeSRef:
			case xltypeMulti:
            {
				if (xlretUncalced == Excel (xlCoerce, &xMulti, 2, (LPXLOPER)xlTab, TempNum (xltypeMulti))) 
				{
					error = xlerrValue;
					return error;
				}
			
                k = 0;
				size = xMulti.val.array.rows;
				vecvec.resize(size);
				for (i = 0; i < xMulti.val.array.rows; i++) 
				{                
					for(j = 0; j <  xMulti.val.array.columns; j++) 
					{
						px = xMulti.val.array.lparray + k;
						++k;
						switch(px->xltype) 
						{
							case xltypeNum:
								vecvec[i].push_back (px->val.num);
							break;
							
                            case xltypeErr:
								error = px->val.err;
							break;
							
                            case xltypeNil:
								error = xlerrValue;
							break;
			
                            default:
								error = xlerrValue;
							break;
						}
					}
				}			

				Excel (xlFree, 0, 1, (LPXLOPER)&xMulti);
            };	
            break;
			
            case xltypeErr:
				error = xlTab->val.err;
			break;
			
            default:
				error = xlerrValue;
			break;
		}
     }
   }

   catch(...)
   {
	   error = xlerrValue;

	   return(error);
   }

   return(error);
}



int XL_getStrCell (LPXLOPER str, int* error, char** reason, 
				   char* r_missing, char* r_inconsistency, 
				   char* r_is_xl_err, CCString& result, 
                   const char* defaultValue, int long_or_float)
{
	XLOPER xMulti;
	LPXLOPER px;
	char* temp = NULL;
	
	// CC_ASSERT (error, MSG_F_MIS_INCONSISTENT_VALUE);
	// CC_ASSERT (reason, MSG_F_MIS_INCONSISTENT_VALUE);

	*error = XL_NO_ERROR;
   
   try
   {
	 /// this is very fundamental as there is no guarantee
	 /// that the pointor is not null!

     if ( str != NULL )
     {
		switch(str->xltype)
		{
			case xltypeRef:
			case xltypeSRef:
			case xltypeMulti:
            {
				if (xlretUncalced == Excel (xlCoerce, &xMulti, 2, (LPXLOPER)str, TempInt (xltypeMulti)))
				{
					*error = xlerrValue;
					*reason = r_inconsistency;
				}
				else
				{
					px = xMulti.val.array.lparray;
				}
            };
			break;

			default:
				px = str;
			break;
		}

		if (*error == XL_NO_ERROR)
		{
			switch(px->xltype)
			{
				case xltypeNum:
                {
					if(long_or_float == LONG_TYPE )
					{
						if( fabs(px->val.num ) > ULONG_MAX )
						{
							Excel (xlFree, 0, 1, (LPXLOPER)&xMulti);
							*error = xlerrValue;
							return (-1);
						}

						temp = new char[(CHARS_TO_PRINT_BYTES (sizeof (unsigned long)) * sizeof (char) + 1)];
						
                        if(!temp)
						{
							Excel (xlFree, 0, 1, (LPXLOPER)&xMulti);
							*error = xlerrValue;
							return (-1);
						}
						
						sprintf (temp, "%ld", (unsigned long)px->val.num );
					}
					else
					{
						const int MAX_SIZE		= 255;
						const double BIG_DOUBLE = 1E+239;	/// 255 - 15 (decimal)- 1(.) = 239 to allow 15 digits for the decimal part
						
                        temp = new char[MAX_SIZE+1];		/// +1 for the null character at the end!
				
                        if(!temp)
						{
							Excel (xlFree, 0, 1, (LPXLOPER)&xMulti);
							*error = xlerrValue;
							return (-1);
						}

						if( fabs( px->val.num ) > BIG_DOUBLE )
						{
							Excel (xlFree, 0, 1, (LPXLOPER)&xMulti);
							*error = xlerrValue;
							return (-1);
						}

						sprintf (temp, "%.15lf", (double)px->val.num);
					}
                };	
                break;

				case xltypeInt:
                {
					/// because we use delete
					/// we should never ever
					/// use malloc
					/// NEVER EVER
					temp = new char[(CHARS_TO_PRINT_BYTES (sizeof (unsigned long)) * sizeof (char) + 1)];
					if(!temp)
					{
						Excel (xlFree, 0, 1, (LPXLOPER)&xMulti);
						*error = xlerrValue;
						return (-1);
					}
				
                    sprintf (temp, "%ld", (unsigned long)px->val.w);
                };
                break;

				case xltypeStr:
                {
					if(!(temp = (char*)XL_StrPascal2StrC (px->val.str)))
					{
						*error = xlerrValue;
						*reason = r_inconsistency;
						break;
					}
                };
                break;

				case xltypeMissing:
				case xltypeNil:
                {
					if(defaultValue)
					{
						/// because we use delete
						/// we should never ever
						/// use malloc
						/// NEVER EVER

						temp = new char[strlen(defaultValue)+1];
						strcpy(temp,defaultValue);
						if(!(temp))
						{
							*error = xlerrValue;
							*reason = r_inconsistency;
						}
					}
					else
					{
						*error = xlerrValue;
						*reason = r_inconsistency;
					}
                };
                break;

				case xltypeErr:
					*error = px->val.err;
					*reason = r_is_xl_err;
				break;

				default:
					*error = xlerrValue;
					*reason = r_missing;
				break;
			}
		}

		switch(str->xltype)
		{
			case xltypeRef:
			case xltypeSRef:
			case xltypeMulti:
				Excel (xlFree, 0, 1, (LPXLOPER)&xMulti);
			break;
			
            default:
			break;
		}

		result.Set (temp);
		
        if (temp)
		{
			delete [] temp;
			temp = NULL;
		}
     }
   }

   catch(...)
   {
	   *error = xlerrValue;

	   return(-1);
   }

   return(0);
}



char** XL_getStrArray (LPXLOPER tabstr, long* size, int* error, char** reason, 
                       char* r_missing, char* r_inconsistency, 
                       char* r_is_xl_err, int long_or_float)
{
	LPXLOPER px;
	XLOPER xMulti;
	int i;
	char** result;

	*error = XL_NO_ERROR;
	*size = 0;

	i = 0;

   try
   {
	 if ( tabstr != NULL )
     {
		switch(tabstr->xltype)
		{
			case xltypeMissing:
			case xltypeNum:
			case xltypeStr:
				*error = xlerrValue;
				*reason = r_missing;
			break;
			
            case xltypeRef:
			case xltypeSRef:
			case xltypeMulti:
            {
				if (xlretUncalced == Excel (xlCoerce, &xMulti, 2, (LPXLOPER)tabstr, TempInt (xltypeMulti)))
				{
					*error = xlerrValue;
					*reason = r_inconsistency;
				}
				else
				{
					*size = xMulti.val.array.rows * xMulti.val.array.columns;
					result = (char**)malloc (*size * sizeof (char*));
					if(!result)
					{
						Excel (xlFree, 0, 1, (LPXLOPER)&xMulti);
						return (NULL);
					}
			
                    for(i = 0; i < *size; i++)
					{
						px = xMulti.val.array.lparray + i;
						CCString str;
						XL_getStrCell(px, error, reason, r_missing, r_inconsistency, r_is_xl_err, str, NULL, long_or_float);
						result[i] = (char*)str;
				
                        if(*error != XL_NO_ERROR)
						{
							Excel (xlFree, 0, 1, (LPXLOPER)&xMulti);
							return (NULL);
						}
					}
				}
	
                Excel (xlFree, 0, 1, (LPXLOPER)&xMulti);
            };	
            break;

			case xltypeErr:
				*error = tabstr->val.err;
				*reason = r_is_xl_err;
		    break;
			
            default:
				*error = xlerrValue;
				*reason = r_missing;
		    break;
		}
     }
   }

   catch(...)
   {
	   *error = xlerrValue;

	   return(NULL);
   }

   return(result);
}



int XL_getStrVectorWD (LPXLOPER tabstr, char** reason, char* r_missing, 
                       char* r_inconsistency, char* r_is_xl_err, 
                       VECTOR<CCString>& vec, VECTOR<CCString>* vecDefault,
                       int long_or_float, const char* singleDefault)
{
	XLOPER xMulti;
	LPXLOPER px;
	WORD i;
	int size = 0;
	int error = XL_NO_ERROR;

	vec.clear ();
	i = 0;

   try
   {
	 if ( tabstr != NULL )
     {
		switch(tabstr->xltype)
		{
  			case xltypeNil:
			case xltypeMissing:
            {
				if (vecDefault)
				{
					for (int i=0;i<vecDefault->size();i++)
						vec.push_back((*vecDefault)[i]);
				}
				else
				{
					error = xlerrValue;
				}
            };
            break;

			case xltypeStr:
			{
				CCString str;
				XL_getStrCell (tabstr, &error, reason, r_missing, r_inconsistency, r_is_xl_err, str, singleDefault, long_or_float);
				vec.push_back (str);
			};
			break;
			case xltypeNum:			
				error = xlerrValue;
				*reason = r_missing;
			break;
			
            case xltypeRef:
			case xltypeSRef:
			case xltypeMulti:
            {
				if (xlretUncalced == Excel (xlCoerce, &xMulti, 2, (LPXLOPER)tabstr, TempInt (xltypeMulti)))
				{
					error = xlerrValue;
					*reason = r_inconsistency;
				}
				else
				{
					size = xMulti.val.array.rows * xMulti.val.array.columns;
					for(i = 0; i < size; i++)
					{
						px = xMulti.val.array.lparray + i;
						CCString str;
						XL_getStrCell (px, &error, reason, r_missing, r_inconsistency, r_is_xl_err, str, singleDefault, long_or_float);
						vec.push_back (str);
						if(error != XL_NO_ERROR)
						{
							Excel (xlFree, 0, 1, (LPXLOPER)&xMulti);
							return error;
						}
					}
				}
			
                Excel (xlFree, 0, 1, (LPXLOPER)&xMulti);
            };
            break;
			
            case xltypeErr:
				error = tabstr->val.err;
				*reason = r_is_xl_err;
			break;
			
            default:
				error = xlerrValue;
				*reason = r_missing;
			break;
		}
     }
   }

   catch(...)
   {
	   error = xlerrValue;

	   return(error);
   }

   return(error);
}



int XL_getStrVector(LPXLOPER tabstr, char** reason, char* r_missing, 
                    char* r_inconsistency, char* r_is_xl_err, VECTOR<CCString>& vec, VECTOR<CCString>* vecDefault, int long_or_float)
{
	return XL_getStrVectorWD( tabstr, reason, r_missing, r_inconsistency, r_is_xl_err, vec, vecDefault, long_or_float, NULL );
}



/// in addition to the function XL_getStrVector
/// this returns the number of rows and columns
int XL_getStrVectorAndSizeWD(LPXLOPER tabstr, char** reason, char* r_missing,
                             char* r_inconsistency, char* r_is_xl_err, 
						     long& nbRows, long& nbColumns,
						     VECTOR<CCString>& vec, 
                             VECTOR<CCString>* vecDefault, int long_or_float, const char* singleDefault )
{
	/// get first the number of columns and rows  	
	XLOPER xMulti;
	int error = XL_NO_ERROR;

	/// to make sure we return some result, we first initialise this to 0
	nbRows = nbColumns = 0;

   try
   {
	 if ( tabstr != NULL )
     {
		switch(tabstr->xltype)
		{
  			case xltypeNil:
			case xltypeMissing:
            {
				if (vecDefault)
				{
					nbRows		= vecDefault->size();
					nbColumns	= 1;
				}
				else
				{
					*reason = r_inconsistency;
					return xlerrValue;
				}
            };
            break;

			/// if it is a reference or a Sref or a Variant...ok
			case xltypeRef:
			case xltypeSRef:
			case xltypeMulti:
            {
				if (xlretUncalced == Excel (xlCoerce, &xMulti, 2, (LPXLOPER)tabstr, TempInt (xltypeMulti)))
				{
					*reason = r_inconsistency;
					return xlerrValue;
				}
				else
				{
					nbRows		= xMulti.val.array.rows;
					nbColumns	= xMulti.val.array.columns;
				}	
			
                Excel (xlFree, 0, 1, (LPXLOPER)&xMulti);
            };
            break;

            //// otherwise, return error
			default:
				*reason = r_missing;
			return xlerrValue;
		}
     }
   }

   catch(...)
   {
	   error = xlerrValue;

	   return(error);
   }

   /// use the standard function
   return XL_getStrVectorWD(tabstr, reason, r_missing, r_inconsistency, r_is_xl_err, vec, vecDefault, long_or_float, singleDefault );
}				



int XL_getStrVectorAndSize(LPXLOPER tabstr, char** reason, char* r_missing,
                           char* r_inconsistency, char* r_is_xl_err, 
						   long& nbRows, long& nbColumns,
						   VECTOR<CCString>& vec, 
                           VECTOR<CCString>* vecDefault, int long_or_float)
{
	return XL_getStrVectorAndSizeWD(tabstr, reason, r_missing, r_inconsistency, r_is_xl_err, 
			nbRows, nbColumns, vec, vecDefault, long_or_float, NULL );
}



CCString XL_StrPascal2StrC(LPSTR str)
{
	int size, i;
	char *temp;

	if (!str) 
	{
       return(NULL);
	}
	
  	size = (int)str[0];

	/// this is to allow long string from excel
	/// since something longer than 128 has a size = 256-size
	if( size <0 )
		size += 256;

	temp = (char *) malloc(sizeof (char) * (size + 1));
	
    if (!temp)
	{
		return NULL;
	}
	
    for(i = 0; i < size; i++) 
	{
		temp[i] = str[i + 1];
	}

    temp[size] = '\0';

	CCString result (temp);

	if(temp)
	{
		free (temp);
	}
	
	return result;
}



LPSTR XL_StrC2StrPascal(const CCString& str)
{
	int str_long, i;
	char* str_result = NULL;
	
	if(!str)
	{
		return NULL;
	}

	str_long = strlen ((const char*)str);

	str_result = (char*)malloc(sizeof (char) * (str_long + 2));
	str_result[0] = (BYTE)str_long;

    for(i = 1; i < str_long + 1; i++)
	{
		str_result[i] = str[i - 1];
	}
	
    str_result[str_long + 1] = '\0';
		
	return ((LPSTR)str_result);
}



int XL_Coordonnate2Rank (int x, int y, int nbcol)
{
	return (nbcol * x + y);
}


BOOL CALLBACK EnumProc (HWND hwnd, EnumStruct* pEnum)
{
  	char rgsz[CLASS_NAME_BUFFER];

  	GetClassName (hwnd, rgsz, CLASS_NAME_BUFFER);

  	if(!lstrcmpi(rgsz, "XLMAIN")) 
    {
    	if(LOWORD((DWORD)hwnd) == pEnum->wLoword) 
		{
	  		pEnum->hwnd = hwnd;
	  		return FALSE;
		}
    }

  	return TRUE;
}


BOOL GetHwnd (HWND* pHwnd)
{
  	XLOPER x;
  
	if (Excel4 (xlGetHwnd, &x, 0) == xlretSuccess) 
    {
    	EnumStruct enm;
    	enm.hwnd = NULL;
    	enm.wLoword = x.val.w;

    	EnumWindows ((WNDENUMPROC)EnumProc, (LPARAM)&enm);

		Excel (xlFree, 0, 1, (LPXLOPER)&x);

      	if(enm.hwnd != NULL) 
		{
	  		*pHwnd = enm.hwnd;
	  		return TRUE;
		}
    }

  	return FALSE;
}



BOOL GetInst(HWND* pInst)
{
  	XLOPER x;
  
	if(Excel4 (xlGetInst, &x, 0) == xlretSuccess) 
    {
    	EnumStruct enm;
    	enm.hwnd = NULL;
    	enm.wLoword = x.val.w;

    	EnumWindows ((WNDENUMPROC)EnumProc, (LPARAM)&enm);

		Excel (xlFree, 0, 1, (LPXLOPER)&x);

      	if(enm.hwnd != NULL) 
		{
	  		*pInst = enm.hwnd;
	  		return TRUE;
		}
    }

  	return FALSE;
}


CCString XL_getCaller ()
{
	XLOPER xSheetName, xRef;
	HWND hWnd;
	char buffer[XL_CALLER_SIZE];

	Excel4 (xlfCaller, (LPXLOPER)&xRef, 0);

	if ( ( (xRef.xltype == xltypeSRef) || (xRef.xltype == xltypeRef) )&& (GetHwnd (&hWnd) == TRUE))
	{
		int lig, col;
		if (xRef.xltype == xltypeSRef)
		{
			lig = xRef.val.sref.ref.rwFirst;
			col = xRef.val.sref.ref.colFirst;
		}
		else
		{
			lig = xRef.val.mref.lpmref->reftbl[0].rwFirst;
			col = xRef.val.mref.lpmref->reftbl[0].colFirst;
		}
		
		Excel4 (xlSheetNm, &xSheetName, 1, (LPXLOPER)&xRef);
		CCString SheetName = XL_StrPascal2StrC (xSheetName.val.str);
				
		sprintf (buffer, "%d$%s_%d_%d", hWnd, (const char*)SheetName, lig, col);

		Excel (xlFree, 0, 1, (LPXLOPER)&xSheetName);
		Excel (xlFree, 0, 1, (LPXLOPER)&xRef);

		return CCString (buffer);
	}

	Excel (xlFree, 0, 1, (LPXLOPER)&xRef);

	strcpy(buffer, "caller unknown");

	return CCString (buffer);
}



CCString XL_getActiveCell ()
{
	XLOPER xSheetName, xRef;
	HWND hWnd;
	char buffer[XL_CALLER_SIZE];

	Excel4 (xlfActiveCell, (LPXLOPER)&xRef, 0);

	if ((xRef.xltype == xltypeSRef) && (GetHwnd (&hWnd) == TRUE))
	{
		int lig = xRef.val.sref.ref.rwFirst;
		int col = xRef.val.sref.ref.colFirst;
		
		Excel4 (xlSheetNm, &xSheetName, 1, (LPXLOPER)&xRef);
		CCString SheetName = XL_StrPascal2StrC (xSheetName.val.str);
				
		sprintf (buffer, "%d$%s_%d_%d", hWnd, (const char*)SheetName, lig, col);

		Excel (xlFree, 0, 1, (LPXLOPER)&xSheetName);
		Excel (xlFree, 0, 1, (LPXLOPER)&xRef);

		return CCString (buffer);
	}

	Excel (xlFree, 0, 1, (LPXLOPER)&xRef);

	strcpy (buffer, "caller unknown");

	return CCString (buffer);
}



void XL_getActiveCellContent (LPXLOPER xRef)
{
	Excel4 (xlfActiveCell, xRef, 0);
}



BOOL CALLBACK EnumProcWiz (HWND hwnd, LPEnumStructWiz pEnum)
{
	char rgsz[CLASS_NAME_BUFFER];

	GetClassName (hwnd, (LPSTR)rgsz, CLASS_NAME_BUFFER);

    if(2 == CompareString (MAKELCID(MAKELANGID(LANG_ENGLISH,SUBLANG_ENGLISH_US),SORT_DEFAULT),
		NORM_IGNORECASE, (LPSTR)rgsz, (lstrlen ((LPSTR)rgsz) > lstrlen ("bosa_sdm_XL"))
		? lstrlen ("bosa_sdm_XL") : -1, "bosa_sdm_XL", -1))
	{
		if(LOWORD((DWORD)GetParent (hwnd)) == pEnum->hwndXLMain)
		{
			pEnum->bFuncWiz = TRUE;
			return FALSE;
		}
	}

	return TRUE;
}



BOOL XL_IsCalledByFuncWiz (void)
{
	XLOPER xHwndMain;
	EnumStructWiz enm;

	if(Excel4 (xlGetHwnd, &xHwndMain, 0) == xlretSuccess)
	{
		enm.bFuncWiz = FALSE;
		enm.hwndXLMain = xHwndMain.val.w;
		EnumWindows ((WNDENUMPROC)EnumProcWiz, (LPARAM)((LPEnumStructWiz)&enm));
		return enm.bFuncWiz;
	}

	return FALSE;
}



int XL_setNumCell( XLOPER& XL_result, double value)
{
	XL_result.xltype	= xltypeNum;
	XL_result.val.num	= value;

	return XL_NO_ERROR;
}


/// function to free vector of string (pascal type)
/// hence allocated with malloc hence the free
void XLOPER_Holder::FreeStringVec( XLOPER& XL_result )
{
	if ( XL_result.xltype == xltypeMulti )
	{
		WORD i,j,
			rows	= XL_result.val.array.columns,
			columns	= XL_result.val.array.rows;
		LPXLOPER pxArray = XL_result.val.array.lparray;

		for( i=0; i<rows; ++i)
			for( j=0; j<columns; ++j )
			{
				/// test the xltypeStr bit!
				if( xltypeStr == ( xltypeStr & pxArray[XL_Coordonnate2Rank (i, j, columns)].xltype ) )
				{
					free( pxArray[XL_Coordonnate2Rank (i, j, columns)].val.str );
				}
			}
	}
}

/// function to free an xloper of type xltypeMulti
int XLOPER_Holder::FreeXLOPER_Multi( XLOPER& XL_result )
{
	/// tries to free memory
	if ( XL_result.xltype == xltypeMulti )
	{
		XLOPER_Holder::FreeStringVec( XL_result );
		HGLOBAL ReturnPtr = GlobalFree( (HGLOBAL) XL_result.val.array.lparray );
	
        if(ReturnPtr)
		{
			/// problem in memory free
			/// need to put a break point here 
			/// and look at the error code!
			int errorcode = GetLastError();
			return XL_ERROR;
		}
	}

	/// other cases are ok!
	return XL_NO_ERROR;
}


/// when constructing the object set the type to nil
XLOPER_Holder::XLOPER_Holder()
{
	itsResult.xltype = xltypeNil;
}



/// Destructor!
XLOPER_Holder::~XLOPER_Holder()
{
	XLOPER_Holder::FreeXLOPER_Multi(itsResult);
}



/// function to set a numeric vector
int XL_setNumVector( XLOPER& XL_result, const VECTOR<double>& vec, 
	int additionalLinesNb, bool fillWithBlank, bool filterSpecificValue, double specificValue )
{
	/// avoid negative additional blank lines!
	if( additionalLinesNb < 0 )
		additionalLinesNb = 0;

	/// function to free an xloper
	if( XLOPER_Holder::FreeXLOPER_Multi(XL_result) == XL_ERROR )
		return XL_ERROR;

	int nbVecRows	= vec.size();
	int nbRows		= nbVecRows + additionalLinesNb;
	int nbColumns	= 1;

	// gestion du cas une seule ligne, pour rajouter des #NA
/*	if (nbRows == 1)
		nbRows++;
*/
	/// allocate the memory
	LPXLOPER pxArray;
	if( !( XL_result.val.array.lparray = pxArray = (LPXLOPER) GlobalAlloc (GMEM_ZEROINIT, nbRows* nbColumns * sizeof (XLOPER) ) ) )
	{
		/// problem with memory allocation!
		return XL_ERROR; 
	}

	XL_result.xltype			= xltypeMulti;
	XL_result.val.array.columns	= nbColumns;
	XL_result.val.array.rows	= nbRows; 

	int i;

	for(i = 0; i<nbVecRows; ++i)
	{
		if(filterSpecificValue && fillWithBlank && vec[i]== specificValue )
		{
			pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].val.str= XL_StrC2StrPascal("");
			pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].xltype |= xlbitDLLFree;
		}
		else
		{
			pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].xltype	= xltypeNum;
			pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].val.num	= vec[i]; 
		}
	}

	/// fill the rest with nothing (typeString with no string)
	for( ; i<nbRows; ++i )
	{
		if( fillWithBlank )
		{
			pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].val.str= XL_StrC2StrPascal("");
			pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].xltype |= xlbitDLLFree;
		}
		else
		{
			pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].xltype = xltypeErr;
			pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].val.err = xlerrNA;
		}
	}

    return XL_NO_ERROR;
}



/// set a vector of string
int XL_setStrVector( XLOPER& XL_result, const VECTOR<CCString>& vecStr, 
	int additionalLinesNb, bool fillWithBlank )
{
	/// avoid negative additional blank lines!
	if( additionalLinesNb < 0 )
		additionalLinesNb = 0;

	/// function to free an xloper
	if( XLOPER_Holder::FreeXLOPER_Multi(XL_result) == XL_ERROR )
		return XL_ERROR;

	int nbVecRows	= vecStr.size();
	int nbRows		= nbVecRows + additionalLinesNb;
	int nbColumns	= 1;

	LPXLOPER pxArray;
	if( !(XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbRows * nbColumns * sizeof (XLOPER)) ) )
	{
		/// problem with memory allocation!
		return XL_ERROR; 
	}

	XL_result.xltype			= xltypeMulti;
	XL_result.val.array.columns = nbColumns;
	XL_result.val.array.rows	= nbRows; 

	int i;
	for(i=0; i <nbVecRows; ++i)
	{
		pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].xltype = xltypeStr;
		pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].val.str = XL_StrC2StrPascal (vecStr[i]);
		pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].xltype |= xlbitDLLFree;
	}

	/// fill the rest with nothing (xltypeNil)
	for( ; i<nbRows; ++i )
	{
		if(fillWithBlank)
		{
			pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].val.str = XL_StrC2StrPascal ("");
			pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].xltype |= xlbitDLLFree;
		}
		else
		{
			pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].xltype = xltypeNil;
		}
	}

	return XL_NO_ERROR;
}

/// line vector version of the XL_setStrVector function

int XL_setStrLineVector( XLOPER& XL_result, const VECTOR<CCString>& vecStr, 
	int additionalColsNb, bool fillWithBlank )
{
	/// avoid negative additional blank lines!
	if( additionalColsNb < 0 )
		additionalColsNb = 0;

	/// function to free an xloper
	if( XLOPER_Holder::FreeXLOPER_Multi(XL_result) == XL_ERROR )
		return XL_ERROR;

	int nbVecCols	= vecStr.size();
	int nbCols		= nbVecCols + additionalColsNb;
	int nbLines	= 1;

	LPXLOPER pxArray;
	if( !(XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbLines * nbCols * sizeof (XLOPER)) ) )
	{
		/// problem with memory allocation!
		return XL_ERROR; 
	}

	XL_result.xltype			= xltypeMulti;
	XL_result.val.array.columns = nbCols;
	XL_result.val.array.rows	= nbLines; 

	int i;
	for(i=0; i <nbVecCols; ++i)
	{
		pxArray[XL_Coordonnate2Rank (0, i, nbLines)].xltype = xltypeStr;
		pxArray[XL_Coordonnate2Rank (0, i, nbLines)].val.str = XL_StrC2StrPascal (vecStr[i]);
		pxArray[XL_Coordonnate2Rank (0, i, nbLines)].xltype |= xlbitDLLFree;
	}

	/// fill the rest with nothing (xltypeNil)
	for( ; i<nbCols; ++i )
	{
		if(fillWithBlank)
		{
			pxArray[XL_Coordonnate2Rank (0, i, nbLines)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (0, i, nbLines)].val.str = XL_StrC2StrPascal ("");
			pxArray[XL_Coordonnate2Rank (0, i, nbLines)].xltype |= xlbitDLLFree;
		}
		else
		{
			pxArray[XL_Coordonnate2Rank (0, i, nbLines)].xltype = xltypeNil;
		}
	}

	return XL_NO_ERROR;
}

/// set two columns first vector of string and vector of numeric
int XL_setStrAndNumVector(XLOPER& XL_result, const VECTOR<CCString>& vecStr, 
                          const VECTOR<double>& vec, 
	                      int additionalLinesNb, bool fillWithBlank )
{
	/// avoid negative additional blank lines!
	if( additionalLinesNb < 0 )
		additionalLinesNb = 0;

	/// function to free an xloper
	if( XLOPER_Holder::FreeXLOPER_Multi(XL_result) == XL_ERROR )
		return XL_ERROR;

	if( vecStr.size() != vec.size() )
		return XL_ERROR;

	int nbVecRows	= vec.size();
	int nbRows		= nbVecRows + additionalLinesNb;
	int nbColumns	= 2;

	LPXLOPER pxArray;
	if (!(XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbRows * nbColumns * sizeof (XLOPER)) ) )
	{
		/// problem with memory allocation!
		return XL_ERROR; 
	}
	XL_result.xltype			= xltypeMulti;
	XL_result.val.array.columns = nbColumns;
	XL_result.val.array.rows	= nbRows; 


	int i;
	for(i = 0; i < nbVecRows; ++i)
	{
		pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].xltype  = xltypeStr;
		pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].val.str = XL_StrC2StrPascal (vecStr[i]);
		pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].xltype |= xlbitDLLFree;
		pxArray[XL_Coordonnate2Rank (i, 1, nbColumns)].xltype  = xltypeNum;
		pxArray[XL_Coordonnate2Rank (i, 1, nbColumns)].val.num = vec[i]; 
	}

	/// fill the rest with nothing (typeString with no string)
	for( ; i<nbRows; ++i )
	{
		if(fillWithBlank)
		{
			pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].val.str = XL_StrC2StrPascal ("");
			pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].xltype |= xlbitDLLFree;
			
			pxArray[XL_Coordonnate2Rank (i, 1, nbColumns)].xltype = xltypeStr;
			pxArray[XL_Coordonnate2Rank (i, 1, nbColumns)].val.str = XL_StrC2StrPascal ("");
			pxArray[XL_Coordonnate2Rank (i, 1, nbColumns)].xltype |= xlbitDLLFree;
		}
		else
		{
			pxArray[XL_Coordonnate2Rank (i, 0, nbColumns)].xltype = xltypeNil;
			pxArray[XL_Coordonnate2Rank (i, 1, nbColumns)].xltype = xltypeNil;
		}
	}

    return XL_NO_ERROR;
}



/// in addition to the function XL_getStrVector
/// this returns the number of rows and columns
int XL_getStrVectorSizeAndType(LPXLOPER tabstr, char** reason, char* r_missing,
                               char* r_inconsistency, char* r_is_xl_err,
                               long& nbRows, long& nbColumns,  
                               VECTOR<CCString>& vec, 
                               VECTOR<CCString>* vecDefault, 
                               VECTOR<long>& type, 
                               VECTOR<long>* typeDefault,
                               int long_or_float)
{
	/// check first that it is really a matrix 
	/// and get the corresponding number of columns and rows  	
	XLOPER xMulti;
	int error = XL_NO_ERROR;

   try
   {
	 switch(tabstr->xltype)
     {
		/// if it is a reference or a Sref or a Variant...ok
		case xltypeMissing:
        {
			if (vecDefault && typeDefault)
			{
				int i = 0;
				for (i=0;i<vecDefault->size();i++)
					vec.push_back((*vecDefault)[i]);
				for (i=0;i<typeDefault->size();i++)
					type.push_back((*typeDefault)[i]);
			}
			else
			{
				error = xlerrValue;
			}
        };
		break;

        case xltypeRef:
		case xltypeSRef:
		case xltypeMulti:
        {
			if (xlretUncalced == Excel (xlCoerce, &xMulti, 2, (LPXLOPER)tabstr, TempInt (xltypeMulti)))
			{
				*reason = r_inconsistency;
				return xlerrValue;
			}
			else
			{
				nbRows		= xMulti.val.array.rows;
				nbColumns	= xMulti.val.array.columns;

				int i;
				int size = nbRows * nbColumns;

				/// clear and allocate memory for vec and type!
				vec.clear();
				vec.reserve( size);
				type.clear();
				type.reserve( size);

				/// single cell object
				LPXLOPER px;
				
				/// fills the object
				for(i = 0; i < size; i++)
				{
					px = xMulti.val.array.lparray + i;
					CCString str;

					/// read the single cell
					XL_getStrCell (px, &error, reason, r_missing, r_inconsistency, r_is_xl_err, str, "", long_or_float);
				
                    if (!str.GetLen()&& px->xltype == xlbitXLFree)
					{
						Excel (xlFree, 0, 1, (LPXLOPER)&xMulti);
						*reason = new char[255];
						sprintf( *reason, ", char too long in cell row=%ld col=%ld max 256", i/nbColumns, i%nbColumns);
						return XL_ERROR;
					}

					vec.push_back (str);

					/// get the single cell type!
					type.push_back(px->xltype);

					if ( error != XL_NO_ERROR )
					{
						Excel (xlFree, 0, 1, (LPXLOPER)&xMulti);
						return error;
					}
				}
			}	

			/// free memory
			Excel (xlFree, 0, 1, (LPXLOPER)&xMulti);
        };
        break;
		
		//// otherwise, return error
		default:
			*reason = r_missing;
			return xlerrValue;
     }
   }

   catch(...)
   {
       error = xlerrValue;

       return(error);
   }

   return(error);
}


/// Set a matrix of string
/// A vector, saving string rows by rows, is passed and converted to an LPXLOPER
int XL_setStrMatrixSizeAndType( XLOPER& XL_result,
                               const VECTOR<CCString>& values,
                               const VECTOR<long>& types,
                               long nbRows, long nbCols )
{
    /// Check vector size consistency
	if ( nbRows*nbCols != values.size()  || values.size() != types.size())
	{
		return XL_ERROR; 
	}

	/// function to free an xloper
	if( XLOPER_Holder::FreeXLOPER_Multi(XL_result) == XL_ERROR )
		return XL_ERROR;

	LPXLOPER pxArray;
	
    if ( !(XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbRows * nbCols * sizeof (XLOPER)) ) )
	{
		/// problem with memory allocation!
		return XL_ERROR; 
	}

	XL_result.xltype			= xltypeMulti;
	XL_result.val.array.columns = nbCols;
	XL_result.val.array.rows	= nbRows; 

	int rowIdx,colIdx;
    int idx=0;
    char* resStr=NULL;
    long theType;

	for(rowIdx=0;rowIdx<nbRows;++rowIdx)
	{
	    for(colIdx=0;colIdx<nbCols;++colIdx)
	    {
            theType = types[idx];
            switch(theType)
            {
              case xltypeStr :
		        pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, nbCols)].val.str = XL_StrC2StrPascal (values[idx]);
		        pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, nbCols)].xltype |= xlbitDLLFree;
              break;

              case xltypeNum :
                resStr=values[idx].c_str();
		        pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, nbCols)].val.num = atof(resStr);
              break;

              case xltypeInt :
                resStr=values[idx].c_str();
		        pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, nbCols)].val.w = (short int)(atof(resStr));
              break;

              case xltypeBool :
                resStr=values[idx].c_str();
                pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, nbCols)].val.boolean = ((int)(atof(resStr)) != 0 ? 1 : 0);
              break;

              case xltypeMissing :
		        pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, nbCols)].val.str = XL_StrC2StrPascal ("");
		        pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, nbCols)].xltype |= xlbitDLLFree;
                theType = xltypeStr;
              break;
            }

            ++idx;

		    pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, nbCols)].xltype = theType;

            delete resStr;
            resStr=NULL;
        }
	}

	return XL_NO_ERROR;
}


/// Set a matrix of string
/// A vector, saving string rows by rows, is passed and converted to an LPXLOPER
/// and complete with blank
int XL_setStrMatrixSizeAndTypeWithOptions( XLOPER& XL_result,
                               const VECTOR<CCString>& values,
                               const VECTOR<long>& types,
                               long nbRows, long nbCols,
							   int additionalLinesNb,
							   bool fillWithBlank )
{
    /// Check vector size consistency
	if ( nbRows*nbCols != values.size()  || values.size() != types.size())
	{
		return XL_ERROR; 
	}

	/// function to free an xloper
	if( XLOPER_Holder::FreeXLOPER_Multi(XL_result) == XL_ERROR )
		return XL_ERROR;

	LPXLOPER pxArray;
	int NbRowsWithBlank = nbRows + additionalLinesNb;
	
    if ( !(XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, NbRowsWithBlank * nbCols * sizeof (XLOPER)) ) )
	{
		/// problem with memory allocation!
		return XL_ERROR; 
	}

	XL_result.xltype			= xltypeMulti;
	XL_result.val.array.columns = nbCols;
	XL_result.val.array.rows	= NbRowsWithBlank; 

	int rowIdx,colIdx;
    int idx=0;
    char* resStr=NULL;
    long theType;

	for(rowIdx=0;rowIdx<nbRows;++rowIdx)
	{
	    for(colIdx=0;colIdx<nbCols;++colIdx)
	    {
            theType = types[idx];
            switch(theType)
            {
              case xltypeStr :
		        pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, nbCols)].val.str = XL_StrC2StrPascal (values[idx]);
		        pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, nbCols)].xltype |= xlbitDLLFree;
              break;

              case xltypeNum :
                resStr=values[idx].c_str();
		        pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, nbCols)].val.num = atof(resStr);
              break;

              case xltypeInt :
                resStr=values[idx].c_str();
		        pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, nbCols)].val.w = (short int)(atof(resStr));
              break;

              case xltypeBool :
                resStr=values[idx].c_str();
                pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, nbCols)].val.boolean = ((int)(atof(resStr)) != 0 ? 1 : 0);
              break;

              case xltypeMissing :
		        pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, nbCols)].val.str = XL_StrC2StrPascal ("");
		        pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, nbCols)].xltype |= xlbitDLLFree;
                theType = xltypeStr;
              break;
            }

            ++idx;

		    pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, nbCols)].xltype = theType;

            delete resStr;
            resStr=NULL;
        }
	}


	/// fill the rest with nothing (typeString with no string)
	for( ; rowIdx<NbRowsWithBlank; ++rowIdx )
	{
		for(colIdx=0;colIdx<nbCols;++colIdx)
	    {
			if( fillWithBlank )
			{
				pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, nbCols)].xltype = xltypeStr;
				pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, nbCols)].val.str= XL_StrC2StrPascal("");
				pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, nbCols)].xltype |= xlbitDLLFree;
			}
			else
			{
				pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, nbCols)].xltype = xltypeNil;
			}
		}
	}
	
	return XL_NO_ERROR;
}


/// Set a matrix of string
/// A vector, saving string rows by rows, is passed and converted to an LPXLOPER
int XL_setNumMatrixSize( XLOPER& XL_result,
                               const VECTOR<double>& values,
                               long nbRows, long nbCols )
{
    /// Check vector size consistency
	if ( nbRows*nbCols != values.size())
	{
		return XL_ERROR; 
	}

	/// function to free an xloper
	if( XLOPER_Holder::FreeXLOPER_Multi(XL_result) == XL_ERROR )
		return XL_ERROR;

	LPXLOPER pxArray;
	
    if ( !(XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, nbRows * nbCols * sizeof (XLOPER)) ) )
	{
		/// problem with memory allocation!
		return XL_ERROR; 
	}

	XL_result.xltype			= xltypeMulti;
	XL_result.val.array.columns = nbCols;
	XL_result.val.array.rows	= nbRows; 

	int rowIdx,colIdx;
    int idx=0;

	for(rowIdx=0;rowIdx<nbRows;++rowIdx)
	{
	    for(colIdx=0;colIdx<nbCols;++colIdx)
	    {
		    pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, nbCols)].val.num = values[idx];
		    pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, nbCols)].xltype = xltypeNum;
			
			++idx;
        }
	}

	return XL_NO_ERROR;
}

/// Set a matrix of string
/// A vector, saving string rows by rows, is passed and converted to an LPXLOPER
int XL_setNumMatrixSizeWithOptions( XLOPER& XL_result,
                               const VECTOR<double>& values,
                               long nbRows, long nbCols,
                               int additionalLinesNb,
                               bool fillWithBlank)
{
    /// Check vector size consistency
	if ( nbRows*nbCols != values.size())
	{
		return XL_ERROR; 
	}

	/// function to free an xloper
	if( XLOPER_Holder::FreeXLOPER_Multi(XL_result) == XL_ERROR )
		return XL_ERROR;

	XLOPER xRes;
	Excel4(xlfCaller, &xRes, 0);

	/// limit ourselves to the maximum numbers of rows of excel 65536
	int callerRowSize =  1; ///65536;
	int callerColSize =  1; /// 256;

	if( xRes.xltype == xltypeSRef && fillWithBlank)
	{
		callerRowSize = xRes.val.sref.ref.rwLast - xRes.val.sref.ref.rwFirst +1;
		callerColSize = xRes.val.sref.ref.colLast- xRes.val.sref.ref.colFirst+1;
	}
	else
	{
		callerRowSize = nbRows;
		callerColSize = nbCols;
	}

	int filledRowSize = callerRowSize<nbRows? callerRowSize : nbRows;
	int filledColSize = callerColSize<nbCols? callerColSize : nbCols;


	LPXLOPER pxArray;
    if ( !(XL_result.val.array.lparray = pxArray = (LPXLOPER)GlobalAlloc (GMEM_ZEROINIT, callerRowSize * callerColSize * sizeof (XLOPER)) ) )
	{
		/// problem with memory allocation!
		return XL_ERROR; 
	}

	XL_result.xltype			= xltypeMulti;
	XL_result.val.array.rows	= callerRowSize; 
	XL_result.val.array.columns = callerColSize;

	int rowIdx,colIdx;
    int idx=0;

	for(rowIdx=0;rowIdx<filledRowSize;++rowIdx)
	{
		idx = rowIdx*nbCols;

	    for(colIdx=0;colIdx<filledColSize;++colIdx)
	    {
		    pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, callerColSize)].val.num = values[idx];
		    pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, callerColSize)].xltype = xltypeNum;
			++idx;
        }
        
		for(;colIdx<callerColSize;++colIdx)
	    {
			if( fillWithBlank )
			{
				pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, callerColSize)].xltype = xltypeStr;
				pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, callerColSize)].val.str= XL_StrC2StrPascal("");
				pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, callerColSize)].xltype |= xlbitDLLFree;
			}
			else
			{
				pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, callerColSize)].xltype = xltypeNil;
			}
		}
	}

    for( ; rowIdx<callerRowSize; ++rowIdx )
	{
		for(colIdx=0;colIdx<callerColSize;++colIdx)
	    {
			if( fillWithBlank )
			{
				pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, callerColSize)].xltype = xltypeStr;
				pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, callerColSize)].val.str= XL_StrC2StrPascal("");
				pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, callerColSize)].xltype |= xlbitDLLFree;
			}
			else
			{
				pxArray[XL_Coordonnate2Rank (rowIdx, colIdx, callerColSize)].xltype = xltypeNil;
			}
		}
	}

	return XL_NO_ERROR;
}


// EOF %M%
