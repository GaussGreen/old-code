#pragma warning (disable : 4786 )

#ifndef _ICM_MATRIX_H_
#define _ICM_MATRIX_H_

#include "ARMKernel\glob\linalg.h"
#include "ARMKernel\glob\dates.h"
#include "ARMKernel\util\merge.h"
#include "ICMKernel/glob/icm_enums.h"
#include "ICMKernel\util\icm_qmatrix.h"
#include "ICMKernel/util/icm_utils.h"

#include <string.h>
#include <stdio.h>
#define unlink _unlink


/*********************************************************************************/
/*! \class  ICM_Matrix icm_matrix.h "icm_matrix.h"
 *  \author 
 *	\version 1.0
 *	\date   June 2003
 *	\file   icm_matrix.h
 *	\brief  defines a matrix with ARM_Vector columns */
/***********************************************************************************/


template <class T> class ICM_Matrix : public ARM_Object 
{
    private:
	vector<T*>	itsMatrix;	
	vector<string> itsColNames;

    public:

		ICM_Matrix():ARM_Object()
		{
			Init();
			SetName(ICM_PWCSHORTCURVE);
		}
		bool empty() const { return itsMatrix.empty() ;}
		void Init(void)
		{
			itsMatrix.clear();
			itsColNames.clear();
		}
		

       virtual ~ICM_Matrix(void)
        {
		   itsMatrix.clear();
		   itsColNames.clear();	
        }


        void BitwiseCopy(const ARM_Object* srcmatx)
        {
            ICM_Matrix<T>* matx = (ICM_Matrix<T>*) srcmatx;
			Empty();

			int size = matx->itsMatrix.size();
			if (size)
			{itsMatrix.resize(size);
			itsColNames.resize(size);}

			for (int i=0;i<size;i++)
			{
				if (matx->itsMatrix[i]) 
				{	
					itsMatrix[i] = matx->itsMatrix[i];
					itsColNames[i] = matx->itsColNames[i];
				}
			}
       }

		inline void	Push(T* value, char* name = NULL)
        {
			int size = itsColNames.size();
			itsMatrix.resize(size+1);
			itsColNames.resize(size+1);
			itsMatrix[size] = value;
			itsColNames[size] = (string)(const char*)name;
        }
		
		// JLA Added
		inline void	Push(const T& value, char* name = NULL)
        {
			int size = itsColNames.size();
			itsMatrix.resize(size+1);
			itsColNames.resize(size+1);
			itsMatrix[size] = new T(value) ;
			itsColNames[size] = (string)(const char*)name;
        }
		inline void	Push(T* value, const std::string&name)
        {
			int size = itsColNames.size();
			itsMatrix.resize(size+1);
			itsColNames.resize(size+1);
			itsMatrix[size] = value ;
			itsColNames[size] = name; 
        }
		
		inline void Empty(void)
        {
		   itsMatrix.clear();
		   itsColNames.clear();	
		}

		inline vector<T*>& GetValue(void) {return itsMatrix;}

		inline T* GetCol(int n) const 
        {
			int size = itsColNames.size(); 
			if ((n<size) && (n>=0))
				return itsMatrix[n];
			else
				return NULL;
		}

		inline T* GetColVect(const char* nomcol) const 
        {
			T* Result;
			int nbcol = GetCol(nomcol);
			Result = GetCol(nbcol);

			return (Result);
		}
		inline T* GetColVect(const std::string& nomcol) const 
        { 
			return GetColVect(nomcol.c_str()); 
		}

		inline int GetCol(const char* nomcol,int mode = 0) const 
        {
			int size = itsColNames.size();
			for (int i = 0; i<size; i++)
			{
			if (!_stricmp(itsColNames[i].c_str(),nomcol))
				return (i);
			}
			return (-1);
		}
		
		inline const char* GetColName(int n) const 
        {
			if (n<itsColNames.size())
				return itsColNames[n].c_str();
			else
				return NULL;
		}

		inline long GetNbCol(void) {return itsColNames.size();}

        void Copy(const ARM_Object* src)
        {
            ARM_Object::Copy(src);

            BitwiseCopy(src);
        }


        virtual ARM_Object* Clone(void)
        {
            ICM_Matrix<T>* theClone = new ICM_Matrix<T>();

            theClone->Copy(this);

            return(theClone);
        }

        //	JLA Added 
		ARM_Object* Clone(void) const 
        {
            ICM_Matrix<T>* theClone = new ICM_Matrix<T>();

            theClone->Copy(this);

            return(theClone);
        }

		void View(char* id, FILE* ficOut)
		{
			FILE* fOut;
			char  fOutName[200];
			int i = 0;
			char tmp[ARM_NB_MAX_CHAR_TERMS];

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

			fprintf(fOut, "\n ======> Matrix :\n\n");
		
			for (i = 0; i<itsColNames.size(); i++)
			{
				fprintf(fOut, "\t %s ", itsColNames[i].c_str()); 
			}

			fprintf(fOut, "\n");
			
			if (itsMatrix.size()>0)
			{
			int size = itsMatrix[0]->GetSize();
			int k = itsMatrix[0]->GetSize();

			for (i = 0; i<size; i++)
			{
				for (k = 0; k<itsMatrix.size(); k++)
				{
				strcpy(tmp,itsColNames[k].c_str());
				tmp[2] = '\0';

				if (itsMatrix[0]->GetName()==ARM_VECTOR)
					fprintf(fOut, "\t %f ", (itsMatrix[k])->Elt(i)); 
				}

				fprintf(fOut, "\n");
			}

			}

			if ( ficOut == NULL )
			{
				fclose(fOut);
			}
		}

};

class ICM_Parameters : public ARM_Object
{
private:
	ICM_Matrix<ARM_Vector>					itsDblParams;
	ICM_Matrix< ICM_QMatrix<string> >		itsStrParams;

public:
	ICM_Parameters() {Init();}
	void Init()	{SetName(ICM_PARAMETERS);}

	inline void	Push(ARM_Object* value,char* name)
    {
		if (value->GetName()==ARM_VECTOR)
			itsDblParams.Push((ARM_Vector*)value,name);
		else itsStrParams.Push((ICM_QMatrix<string>*)value,name);
    }
	inline void	Push(ARM_Object* value,const std::string& name)
    {
		if (value->GetName()==ARM_VECTOR)
			itsDblParams.Push((ARM_Vector*)value,name);
		else itsStrParams.Push((ICM_QMatrix<string>*)value,name);
    }
	bool empty() const { return itsDblParams.empty() && itsStrParams.empty(); }
	inline ARM_Vector* GetColVect(char* nomcol) { return  itsDblParams.GetColVect(nomcol);}
	inline ARM_Vector* GetColVect(const std::string&nomcol) const { return  itsDblParams.GetColVect(nomcol);}
	inline ICM_QMatrix<string>* GetColVectStr(char* nomcol) { return  itsStrParams.GetColVect(nomcol);}
	inline ICM_QMatrix<string>* GetColVectStr(const std::string& nomcol) const { return  itsStrParams.GetColVect(nomcol);}

    void Copy(const ARM_Object* src)
    {	ARM_Object::Copy(src);
		BitwiseCopy(src);}

    ARM_Object* Clone(void)
    { ICM_Parameters* theClone = new ICM_Parameters();
      theClone->Copy(this);
      return(theClone);}

    void BitwiseCopy(const ARM_Object* srcmatx)
    {ICM_Parameters* matx = (ICM_Parameters*) srcmatx;
	itsDblParams = 	matx->itsDblParams;
	itsStrParams = 	matx->itsStrParams;}

	void View(char* id, FILE* ficOut)
	{
		FILE* fOut;
		char  fOutName[200];
		int i = 0;
		if ( ficOut == NULL )
		{
			ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
			(void) unlink(fOutName);
			fOut = fopen(fOutName, "w"); 
		}
		else
		{	fOut = ficOut;} 

		fprintf(fOut, "\n-----------------------------------------------------------\n");
		fprintf(fOut, "\n Parameters :\n");
		fprintf(fOut, "\n-----------------------------------------------------------\n\n");
		fprintf(fOut, "\n Parameters double :\n\n");
		itsDblParams.View(id,fOut);

		fprintf(fOut, "\n Parameters string :\n\n");
		itsStrParams.View(id,fOut);

		if ( ficOut == NULL ){fclose(fOut);}
	}

	//	JLA Added 
	ARM_Object* Clone(void) const 
	{
		ICM_Parameters* theClone = new ICM_Parameters();
		theClone->Copy(this);
		return(theClone);
	}

	ICM_Matrix<ARM_Vector>*	GetDblParams() { return &itsDblParams;}
	ICM_Matrix< ICM_QMatrix<string> >* GetStrParams() { return &itsStrParams;}

	//JLA. Simple Helpers to get/set params
	//		returns true/false or throw an error if not found, 
	//		according to throwOnError value.
	bool getParam(const std::string&paramName,double&ret,bool throwOnError=true) const ;
	bool getParam(const std::string&paramName,long&ret,bool throwOnError=true) const; 
	bool getParam(const std::string&paramName,ARM_Vector&ret,bool throwOnError=true) const; 
	//	those versions will throw 
	long	getLong(const std::string&paramName) const
	{
		long ret; getParam(paramName,ret,true); return ret; 
	}
	double  getDouble(const std::string&paramName) const
	{
		double ret; getParam(paramName,ret,true); return ret; 
	}
	//	
	template <class T> void setParam(const std::string&paramName,const T& v)
	{
		ARM_Vector* item = new ARM_Vector(1); 
		item->Elt(0)=v ;
		Push(item,paramName); 
	}
};


class ICM_PropertyList : public ARM_Object
{
	class abstract_t 
	{ 
	public: virtual ~abstract_t() {}
	public: virtual abstract_t* clone() const =0; 
	private: virtual void toStream(std::ostream&o) const=0; 
	friend std::ostream& operator<<(std::ostream&o,const abstract_t&item) { item.toStream(o) ; return o; }
	}; 
	template<class T> class concrete_t : public abstract_t
	{
	private:
		T t_; 
	public:
		concrete_t(const T&item) :  t_(item) {}
		virtual ~concrete_t() {}
		concrete_t(const concrete_t&ref) : abstract_t(ref), t_(ref.t_) {}
		virtual abstract_t* clone() const { return new concrete_t(*this); }
		virtual void toStream(std::ostream&o) const { o<<t_; }
		T& item() const { return t_; }
	private:
		concrete_t& operator=(const concrete_t&ref); 
	}; 	
	typedef std::map<std::string,abstract_t*> cont_t; 
private:
	cont_t itsValues; 
public:
	ICM_PropertyList() ;
	virtual ~ICM_PropertyList(); 
	ICM_PropertyList(const ICM_PropertyList &ref); 
private:ICM_PropertyList& operator=(const ICM_PropertyList &ref); 
public:
	virtual ARM_Object* Clone() ; 
public:
	template<class T> void set(const std::string&name,T&ref)
	{
		abstract_t* item = new concrete_t<T>(ref); 
		cont_t::iterator it = itsValues.find(name); 
		if (it==itsValues.end()) itsValues[name]=item;  
		else { delete  it->second ; it->second = item; }
	}
	template<class T> void get(const std::string&name,T&ref) const 
	{
		cont_t::const_iterator it = itsValues::find(name); 
		if (it==itsValues.end()) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_PropertyList: "<<name<<" attribute does not exist"); 
		concrete_t<T>* item = dynamic_cast<concrete_t<T>*>(*it); 
		it (item==0) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_PropertyList: "<<name<<" is not of proper type"); 
		ref=item.item(); 
	}
	bool exists(const std::string&name) const; 
	virtual void View(char* id = NULL, FILE* ficOut = NULL) ;
}; 

//	--------------------------------------------------------------------
inline 	
ICM_PropertyList::ICM_PropertyList() : ARM_Object()
{}
//	--------------------------------------------------------------------
inline 
bool 
ICM_PropertyList::exists(const std::string& name) const 
{
	cont_t::const_iterator it = itsValues.find(name); 
	return (it!=itsValues.end()) ;
}
 
#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
