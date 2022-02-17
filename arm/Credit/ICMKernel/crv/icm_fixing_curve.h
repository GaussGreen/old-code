#ifndef _ICM_Fixing_Curve_CURVE_H
#define _ICM_Fixing_Curve_CURVE_H

#include <string>
#include <map>

#include "ICMKernel\util\icm_macro.h"
// #include "irindex.h"
class ARM_IRIndex; 

/*********************************************************************************/
/*! \class  ICM_Fixing_Curve ICM_Fixing_Curve.h "ICM_Fixing_Curve.h"
 *  \author cChalaye
 *	\version 1.0
 *	\date   22 Nov 2006
 *	\file   ICM_Fixing_Curve.h
 *	\brief  Default Curve */
/***********************************************************************************/

class ICM_Fixing_Curve : public ARM_Object
{
    private:
		std::string itsIndexName;	
		ARM_IRIndex* itsIndex;			//	0-1 , aggregated
		ARM_Date itsAsOfDate;
		std::map<ARM_Date, double>  itspMFixing;

    public:
		ICM_Fixing_Curve(const std::string& Index,
						const ARM_Date& AsOfDate,
						const std::map<ARM_Date, double>& aMap,const  ARM_IRIndex*index);
		ICM_Fixing_Curve(const std::string& Index,
						const ARM_Date& iAsOfDate,
						const ARM_Date& aDate, double aValue,const  ARM_IRIndex*index );
		ICM_Fixing_Curve(void) {Init();}

		void Set (const std::string& Index,
						const ARM_Date& AsOfDate,
						const std::map<ARM_Date, double>& aMap,const  ARM_IRIndex*index );
		void Set (const std::string& Index,
						const ARM_Date& AsOfDate,
						const ARM_Date& aDate, double aValue,const  ARM_IRIndex*index );


		void Init();

		void View(char* id = NULL, FILE* ficOut = NULL );

		double getValue(const ARM_Date&valueDate) const; 
		bool hasValue(const ARM_Date&valueDate) const; 
		/**const std::map<ARM_Date, double>& getFixingMap() const {
			return *itspMFixing;
		}
		
		std::map<ARM_Date, double>* getPFixingMap() const {
			return itspMFixing;
		}**/ 


		const std::string& GetIndexName() const { return itsIndexName;}
		const ARM_Date& GetAsOfDate() const { return itsAsOfDate; }

		bool hasIndex() const { return itsIndex!=NULL; }
		const ARM_IRIndex& getIndex() const 
		{ 
			if (!itsIndex) ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Fixing_Curve: index is not defined"); 
			return *itsIndex; 
		}
	private:
		void SetIndexName(const std::string& aIndexName) { itsIndexName = aIndexName;}  
		void SetAsOfDate(const ARM_Date& aDate) { itsAsOfDate = aDate;} 
	public:
		void SetMap(const std::map<ARM_Date, double>&apMap ) 
		{
			itspMFixing = apMap; 
		}

		void AddDateValue(const ARM_Date& aDate, double value) ;
		
		void BitwiseCopy(const ARM_Object* src);
        void Copy(const ARM_Object* srcZc)
		{
			ARM_Object::Copy(srcZc);
			BitwiseCopy(srcZc);
		}

        ARM_Object* Clone(void)
        {
            ICM_Fixing_Curve* theClone = new ICM_Fixing_Curve();
            theClone->Copy(this);
            return(theClone);
        }
		
		~ICM_Fixing_Curve() ; 
private:
	ICM_Fixing_Curve(const ICM_Fixing_Curve&ref);			//NA 
	ICM_Fixing_Curve& operator=(const ICM_Fixing_Curve&) ; //NA 
};


#endif /*---- End of file ----*/
