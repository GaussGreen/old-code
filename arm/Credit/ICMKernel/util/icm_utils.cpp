#pragma warning (disable:4786)

#include "ICMKernel\inst\icm_credit_index.h"
#include "ICMKernel\util\icm_utils.h"
#include "ICMKernel\util\icm_macro.h"
#include "ICMKernel\glob\icm_enums.h"
#include "ARMKernel\inst\swapleg.h"
#include "ARMKernel\crv\volint.h"
#include "ARMKernel\util\merge.h"
#include "ICMKernel\glob\icm_maths.h"
#include "ICMKernel\crv\icm_defaultcurve.h"
#include "ICMKernel\glob\icm_correlation.h"

#include <set>


//Used to delete char**
/**
void FreePointerTabChar(char**& Tab,const int& size)
{
	if (Tab == NULL) return;

	for (int i=0; i<size; i++)
	{
		if (Tab[i])	delete[] Tab[i];
		Tab[i] = NULL;
	}

	if (Tab) delete[] Tab;

	Tab = NULL;
}

//Used to copy char**
char** CopyTabString(char** Tab,const int& size)
{
	if (Tab == NULL) return NULL;

	char** NewTab = new char*[size];

	for (int i=0; i<size; i++)
	{
		NewTab[i] = new char[_size_zclabel_];
		strcpy(NewTab[i],Tab[i]);
	}

	return (NewTab);
}

int FindRowInTab(char* chaine, char** Tab,const int& size)
{
	if (Tab == NULL) return NULL;

	for (int i=0; i<size; i++)
	{
		if (strcmp(Tab[i],chaine)==NULL)
			return (i);
	}

	return -999;
}

**/ 
void MergeVector(const vector<double>& v1,
				 const vector<double>& v2,
				 vector<double>& vout)
{
	int i=0;
	set<ORDER_> mOrder;
	vout.clear();

	for (i=0; i<v1.size(); i++)
	{
		ORDER_ order(v1[i]);
		if (mOrder.find(order) == mOrder.end()) mOrder.insert(order);
	}

	for (i=0; i<v2.size(); i++)
	{
		ORDER_ order(v2[i]);
		if (mOrder.find(order) == mOrder.end()) mOrder.insert(order);
	}

	vout.resize(mOrder.size());
	i=0;
	for(set<ORDER_>::iterator iterateur = mOrder.begin(); iterateur != mOrder.end(); ++ iterateur)
	{
		vout[i] = iterateur->id;
		i++;
	}	
}	


/** 
void ComputeFwdDates(ARM_Vector* flowStartDates, 
					 ARM_Vector* flowEndDates,
                     ARM_Vector* paymentDates, 
					 ARM_Vector* resetDates,
					 ARM_IRIndex* irIndex,
					 ARM_Vector*& FwdStartDates,
					 ARM_Vector*& FwdEndDates)
{
	int size = flowStartDates->GetSize();
	ARM_Vector* intDays = new ARM_Vector(size,30.);
	ARM_Vector* fwdRates = new ARM_Vector(size,1.);

	ARM_ReferenceValue* Notional = new ARM_ReferenceValue(100 , 1 // price 
	, 0 // 
	);


    ARM_SwapLeg  SwapLeg((ARM_Vector*)flowStartDates->Clone(), 
						(ARM_Vector*)flowEndDates->Clone(),
						(ARM_Vector*)paymentDates->Clone(), 
						(ARM_Vector*)resetDates->Clone(), 
						intDays, 
						fwdRates,
						Notional, 
						irIndex);

	FwdStartDates = (ARM_Vector*) SwapLeg.GetFwdRateStartDates()->Clone();
	FwdEndDates = (ARM_Vector*) SwapLeg.GetFwdRateEndDates()->Clone();
}
**/ 

 
// JLA 
//	---------------------------------------------------------------------------------------
std::ostream& operator<<(std::ostream&o,const ARM_Vector&v)
{
	o<<"ARM_Vector "<<v.GetNumLines()<< "x" << v.GetNumCols() << std::endl; 
	for(unsigned int i=0;i<v.GetNumCols();i++) o<<i<<":"<<v.Elt(i)<<std::endl ;
	return o ;
}
//	---------------------------------------------------------------------------------------
/*
std::ostream& operator<<(std::ostream&o,const ARM_Date&d)
{
	o<<d.toString('/'); 
	return o ;
}
*/

 // Onin
double Mean(const vector<double>& V)
{
	double M = 0.;
	int size = V.size();
	if (size != 0)
	{
		for (int i = 0; i<size;i++)
			M += V[i];

		M *= 1./(double)size;
	}

	return M;
}

/** 
double vector_max(const vector<double>& V)
{
	if (V.size()==0) return 0; 
	double M = V[0];

	//if (V.size() != 0)
	//{
	for (int i = 0; i<V.size();i++) {if (V[i]>M) M=V[i];}
//}

	return M;
}

double vector_min(const vector<double>& V)
{
	if (V.size()==0) return 0; 
	double M = V[0];

	//if (V.size() != 0)
	//{
		for (int i = 0; i<V.size();i++) {if (V[i]<M) M=V[i];}
	//}

	return M;
}

double vector_sum(const vector<double>& V)
{
	double out=0.;
	for (int i = 0; i<V.size();i++) {out+=V[i];}
	return (out);
}
**/ 

double Sigma(const vector<double>& V, double Moy)
{
	double Sigma = 0.;
	double M = Moy;
	
	if (fabs(M + 999.)<1e-8)
		M = Mean(V);

	int size = V.size();
	if (size > 1)
	{
		for (int i =0; i<size;i++)
			Sigma += (V[i] - M) * (V[i] - M);

		Sigma *= 1.0/(size - 1.0);
		Sigma = sqrt(Sigma);
	}
	return Sigma;
}

//CC Sort of ARM_Vector for the vector of string.
vector<string> SortVectorString(const vector<string>& vToSort, 
									 vector<int>& vPermut,
									 const  int  incOrdec)
{
    
	int i, j, inc, k=0, d=0;
    vector<string> array = vToSort; // copy
    int n = array.size();
	vPermut.reserve(n);
	for (i=0;i< n; i++) 
		vPermut[i]=(int) i;
    string   v;
    inc = 1;
    do
    {
        inc *= 3; 
        inc++;
    }
    while (inc <= n);
    do
    {
        inc /= 3; 
        for (i=inc+1; i<=n; i++) 
        {
            v = array[i-1];
            d = vPermut[i-1];
            j=i;       
            while (array[j-inc-1]  > v )
            {
                array[j-1] = array[j-inc-1];
				vPermut[j-1] = vPermut[j-inc-1];
				j -= inc;
                if (j <= inc) break;
            }
            array[j-1] = v;
			vPermut[j-1] = d;
			
        }
    }
    while (inc>1);

	if (incOrdec == K_DECREASING)
	{
		for (i = 0; i < n / 2; i++)
		{
			string  tmp = array[i];
			int ti = vPermut[i];
			array[i] = array[n-i-1];
			vPermut[i] = vPermut[n-i-1];
			array[n-i-1] = tmp;
			vPermut[n-i-1] = ti;
		}
	}	
    return(array);
}
// CC from ARM_Vector::Sort of ARM
ARM_Vector SortWithConstraints(const ARM_Vector& vToSort, const vector<string>& vConstraintes, 
							   ARM_Vector& vPermut, const int incOrdec)
{
	int i, j, inc, k=0, d=0;
	ARM_Vector array(vToSort);
	vector<string> dynConst(vConstraintes);
	int n = array.GetSize();
	vPermut.Resize(n);
	for (i=0;i< n; i++) 
		vPermut[i]=(int) i;
	double    v;
	string ll;
	inc = 1;
	do
	{
        inc *= 3; 
        inc++;
    }
    while (inc <= n);

    do
    {
        inc /= 3; 
        for (i=inc+1; i<=n; i++) 
        {
            v = array[i-1];
            d = vPermut[i-1];
			ll = dynConst[i-1];
            j=i;
            
            while (Inf2Contraintes(array[j-inc-1] , v, dynConst[j-inc-1], ll) )
            {
                array[j-1] = array[j-inc-1];
				dynConst[j-1] = dynConst[j-inc-1];
				vPermut[j-1] = vPermut[j-inc-1];
				j -= inc;
                if (j <= inc) break;
            }
            array[j-1] = v;
			vPermut[j-1] = d;
			dynConst[j-1] = ll;
			
        }
	}
	while (inc>1);

	if (incOrdec == K_DECREASING)
	{
		for (i = 0; i < n / 2; i++)
		{
			double tmp = array[i];
			int ti = vPermut[i];
			array[i] = array[n-i-1];
			vPermut[i] = vPermut[n-i-1];
			array[n-i-1] = tmp;
			vPermut[n-i-1] = ti;
		}
	}
	
    return(array);

}
bool Inf2Contraintes( const double& dref, const double& dtest, 
					 const string& sref, const string& stest)
{
	if ( dtest  - dref> DB_TOL) {
		return true;
	}else if (fabs (dtest  - dref)< DB_TOL) {
		if (sref < stest) // test on the issuers names
			return true;
		return false;
	}
	return false;
}
 
// vecteur colonne de vecteur ligne // pas testé !!
vector<string> concatStringVector(const vector<vector<string> >& VVString){

	if ( VVString.size() == 0)
		ICMTHROW(ERR_INVALID_ARGUMENT,"concatStringVector : vector size is nulle ");
	vector<string> resultVector(VVString.size());
	for (int r= 0; r< VVString.size(); r++){
		string stringConct("");
		for ( int c = 0; c< VVString[r].size() -1 ; c++) {
			stringConct += VVString[r][c] + "_";
		}
		resultVector[r] = stringConct + VVString[r][c];
	}
	return resultVector;
}

void SearchCoordForVolatility(ARM_VolLInterpol* volcurve,const std::string&  info,int& nthLine,int& nthCol)
{
	nthLine = nthCol = -1;

	//On a une structure ISSUER_STRIKE_MATU

	int  separate = '|';
	char TmpStrike[255];memset(TmpStrike,'\0',sizeof(TmpStrike));
	char TmpMatu[255];memset(TmpMatu,'\0',sizeof(TmpMatu));
	strcpy(TmpStrike,info.c_str());
	char* pdest = strchr(TmpStrike, separate );
	if (pdest == NULL) return;
	*pdest = '\0';
	strcpy(TmpMatu,pdest+1);

	ARM_Vector* Strikes = volcurve->GetStrikes();
	ARM_CRV_TERMS* ExpiryTerms = &(volcurve->itsYearTermsX);
	ARM_Matrix* Volatilities = volcurve->GetVolatilities();
	int nLines = Volatilities->GetNumLines();
    int nCols  = Volatilities->GetNumCols();	

	double Strike = atof(TmpStrike);

	for (int i=0; i<nCols; i++)
	{
		if (fabs(Strike-Strikes->Elt(i))<DB_TOL)
		{
			nthCol = i;
			for (int j=0; j<nLines; j++)
			{
			if (strcmp(TmpMatu,(*ExpiryTerms)[j])==NULL)
			{
				nthLine = j;
				return;
			}
			}
		}
	}
}

//
//		[K:1]					for Tenor (...)
//		[date:d:20060612]		for Expiry Maturity
//		[yf:0.25687]			for Expiry Year fraction	
//		[term:1Y]				for Expiry Tenor
//
//		[I#:1]					for row (expiry)
//		[J#:1]					for col (expiry)
//	output
//		line,col>0
//		line=0	: entire line
//		col=0	: entire col
//
//		otherwise Exception. 
//
//	also supported: 
//		[K:-]				returns 0 for Col
//		[M:d:-]				returns 0 for Row
//	

static inline std::string extract(const std::string& src,const std::string& key)
{
	std::string::size_type KBegin = src.find(key); 
	if (KBegin==std::string::npos) return ""; 
	std::string::size_type KEnd = src.find("]",KBegin); 
	if (KBegin==std::string::npos) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Can't read "<<key<<" for "<<src); 
	return src.substr(KBegin+key.size(),KEnd-KBegin-key.size()); 
}
void 
SearchCoordForVolatility(const ARM_VolCurve&vol_,const std::string& info,int& foundLine,int& foundCol)
{
	// (arm design) 
	const ARM_VolLInterpol&vol=dynamic_cast<const ARM_VolLInterpol&>(vol_);  

	// handle the old case. 
	if (info.find("|")!=std::string::npos) 
		SearchCoordForVolatility(&unconst(vol),info,foundLine,foundCol); 

	//
	//	looking for strikes. 
	//
	foundCol=-2 ; // not found yet 
	std::string item ; 
	item = extract(info,"[K:"); 
	if (!item.empty())
	{
		if (item=="-") foundCol=-1; 
		else 
		{
			std::stringstream sstr; sstr<<item; 
			double strikeValue; sstr>>strikeValue; 
			int pos =-1 ; 
			if (vol.GetStrikes()) pos = vol.GetStrikes()->find(strikeValue); 
			if (pos==-1) 
				ICMTHROW(ERR_INVALID_ARGUMENT,"Can't find strike value: "<<item); 
			foundCol=pos; 
		}
	}

	item = extract(info,"[J#:"); 
	if (!item.empty() && foundCol==-2) 
	{
		if (item=="-") foundCol=-1; 
		else 
		{
			std::stringstream sstr; sstr<<item; 
			int strikePos; sstr>>strikePos; 
			foundCol=strikePos -1 ; 
				// bounds will be check at heding time ?
		}
	}
	// finally:
	//
	if (foundCol==-2) 
		ICMTHROW(ERR_INVALID_ARGUMENT,"Can't find strike: "<<info); 
	//
	//	now looking for matu
	//
	foundLine=-2; 
	item = extract(info,"[term:"); 
	if (!item.empty()) 
	{
		if (item=="-") foundLine=-1;  
		else 
		{
			double yf= CountYears(KACTUAL_365,vol.GetAsOfDate(),
				AddPeriod(
					vol.GetAsOfDate(),
					item,
					vol.GetCurrency()->GetCcyName(),
					false,
					dynamic_cast<const ICM_Credit_Index&>(vol.GetIndex()).GetAdjForTenor())
					);
			int pos = -1; 
			if (vol.GetExpiryTerms()) pos=vol.GetExpiryTerms()->find(yf); 
			if (pos==-1) 
				ICMTHROW(ERR_INVALID_ARGUMENT,"Can't find term value: "<<item<<" yf="<<yf); 
			foundLine=pos; 
		}
		return; 
	}
	item = extract(info,"[date:"); 
	if (!item.empty() && foundLine==-2) 
	{
		if(item=="-") foundLine=-1; 
		else 
		{
			//	We shall transpose from YYYYMMDD to YearFractions
			ARM_Date arg(item.c_str(),"YYYYMMDD"); 
			double yf= CountYears(KACTUAL_365,vol.GetAsOfDate(),arg); 
			int pos = -1; 
			if (vol.GetExpiryTerms()) pos=vol.GetExpiryTerms()->find(yf); 
			if (pos==-1) 
				ICMTHROW(ERR_INVALID_ARGUMENT,"Can't find term value: "<<item<<" yf="<<yf); 
			foundLine=pos; 
		}
		return; 
	}
	item = extract(info,"[yf:"); 
	if (!item.empty() && foundLine==-2) 
	{
		if (item=="-") foundLine=-1; 
		{
			//	This is direct yf input. 
			std::stringstream sstr ; sstr<<item; 
			double yf ; sstr>>yf; 
			int pos = -1; 
			if (vol.GetExpiryTerms()) pos=vol.GetExpiryTerms()->find(yf); 
			if (pos==-1) 
				ICMTHROW(ERR_INVALID_ARGUMENT,"Can't find term value: "<<item<<" yf="<<yf); 
			foundLine=pos; 
		}
		return; 

	}
	item = extract(info,"[I#:"); 
	if (!item.empty() && foundLine==-2) 
	{
		if (item=="-") foundLine=-1; 
		else 
		{
			//	This is direct # input. 
			std::stringstream sstr ; sstr<<item; 
			int  pos ; sstr>>pos; 
			foundLine=pos-1; 
		}
		return; 
	}
	ICMTHROW(ERR_INVALID_ARGUMENT,"Can't find maturities : "<<info); 	
}



double VectorInterpol(const vector<double>& X,const vector<double>& Y, const double& x, int Method)
{
	switch (Method)
	{
	case K_LINEAR:
	default:
		return LinearVectorInterpol(X,Y,x);
	}
}

double LinearVectorInterpol(const vector<double>& X,const vector<double>& Y, const double& x)
{
	if (X.size()!=Y.size()) ICMTHROW(ERR_INVALID_ARGUMENT,"LinearVectorInterpol:vector have different sizes");
	if (X.size()==0) ICMTHROW(ERR_INVALID_ARGUMENT,"LinearVectorInterpol:empty vector");
	if (X.size()==1) return Y[0];
	int il=0,inf=CREDIT_DEFAULT_VALUE,sup=CREDIT_DEFAULT_VALUE; 

	//recherche de l'intervalle ]X[inf],X[sup][ contenant x
	for (il=0; il<X.size(); il++)
	{ if CHECK_EQUAL(X[il],x) return Y[il];
	  if (X[il]<x) inf=il; else break;}

	int inf2 = 0;
	if (inf>=0) inf2= inf; 

	for (il=inf2; il<X.size(); il++)
	{if (X[il]> x) {sup=il;break;}}

	if (inf==CREDIT_DEFAULT_VALUE)
	{inf = 0;
	 sup = inf+1;}		
	else if (sup==CREDIT_DEFAULT_VALUE)
	{ sup = X.size()-1;
	  inf = sup-1;}	

	if  (!(X[sup]-X[inf])) return Y[inf];
	return (Y[sup]-Y[inf])/(X[sup]-X[inf])*(x-X[inf])+Y[inf];
}


double FlatVectorInterpol(const vector<double>& X,const vector<double>& Y, const double& x, bool valueinf)
{
	if (X.size()!=Y.size()) ICMTHROW(ERR_INVALID_ARGUMENT,"LinearVectorInterpol:vector have different sizes");
	if (X.size()==0) ICMTHROW(ERR_INVALID_ARGUMENT,"LinearVectorInterpol:empty vector");
	if (X.size()==1) return Y[0];
	int il=0,inf=CREDIT_DEFAULT_VALUE,sup=CREDIT_DEFAULT_VALUE; 

	//recherche de l'intervalle ]X[inf],X[sup][ contenant x
	for (il=0; il<X.size(); il++)
	{ if CHECK_EQUAL(X[il],x) return Y[il];
	  if (X[il]<x) inf=il; else break;}

	int inf2 = 0;
	if (inf>=0) inf2= inf; 

	for (il=inf2; il<X.size(); il++)
	{if (X[il]> x) {sup=il;break;}}

	if (inf==CREDIT_DEFAULT_VALUE)
	{inf = 0;
	 sup = inf+1;}		
	else if (sup==CREDIT_DEFAULT_VALUE)
	{ sup = X.size()-1;
	  inf = sup-1;}	

	if (valueinf) 
	{ return Y[inf];}
	else
	{ return Y[sup];}
}

void Bornes(const vector<double>& X,const double& x, int & _inf, int & _sup,bool &equal)
{
	if (X.size()==0) ICMTHROW(ERR_INVALID_ARGUMENT,"LinearVectorInterpol:empty vector");
	if (X.size()==1)
	{
		_inf = _sup = 0 ;
		equal = true;
		return;
	}
	int il=0,inf=CREDIT_DEFAULT_VALUE,sup=CREDIT_DEFAULT_VALUE; 

	//recherche de l'intervalle ]X[inf],X[sup][ contenant x
	for (il=0; il<X.size(); il++)
	{ if CHECK_EQUAL(X[il],x) {_inf = _sup = il;
								equal = true;
								return;};
	  if (X[il]<x) inf=il; else break;}

	int inf2 = 0;
	if (inf>=0) inf2= inf; 

	for (il=inf2; il<X.size(); il++)
	{if (X[il]> x) {sup=il;break;}}

	_inf = inf;	
	_sup = sup;	
}

// ----------------------------------------------------------------
// Generation du schedule d'integration
// ----------------------------------------------------------------
ARM_Vector* GenerateIntSch(const double& start,const double& end, const int& step)
{
	int Nb_Step = (int) ((end - start + 1.) / (double)step);
	double schedule = start;
	ARM_Vector* Integration_Schedule = NULL;

	if ((Nb_Step*step+(int)start)>((int)end)) Nb_Step--;

	if (Nb_Step<2)
	{	
		Integration_Schedule = new ARM_Vector(2,0.);
		Integration_Schedule->Elt(0)=start;
		Integration_Schedule->Elt(1)=end;
	}
	else
	{
		Integration_Schedule = new ARM_Vector(Nb_Step,0.);

		for (int i=0; i<Nb_Step; i++)
		{
		Integration_Schedule->Elt(i)=schedule;
		schedule += step;
		}

		Integration_Schedule->Elt(Nb_Step-1)=end;
	}

	return (Integration_Schedule);
}

ARM_Vector* GenerateIntSchInludingOtherSchedule(const double& start,
												const double& end, 
												ARM_Vector* sch,
												const int& step)
{
	ARM_Vector* Schedule2 = NULL;
	ARM_Vector* gensch = GenerateIntSch(start,end,step);
	MergeDates(&Schedule2,gensch,sch);
	int size = 0,i=0,k=0,index = 0;

	ARM_Vector ID(Schedule2->GetSize(),1.);
	for (index=0; index<sch->GetSize();index++)
	for (i=0; i<Schedule2->GetSize();i++)
	{if (fabs(Schedule2->Elt(i)-sch->Elt(index))<DB_TOL) {ID.Elt(i)=0.;}}

	for (index=0; index<sch->GetSize();index++)
	for (i=0; i<Schedule2->GetSize();i++)
	{	if ((fabs(Schedule2->Elt(i)-sch->Elt(index))<step) && (ID.Elt(i)))
		{size++;
		 Schedule2->Elt(i) = CREDIT_DEFAULT_VALUE;}}

	ARM_Vector* Schedule3 = new ARM_Vector(Schedule2->GetSize()-size,0.);
	for (i=0; i<Schedule2->GetSize();i++)
	{
		if (Schedule2->Elt(i)!=CREDIT_DEFAULT_VALUE)
		{Schedule3->Elt(k)=Schedule2->Elt(i);
		k++;}
	}

	if (Schedule2) delete Schedule2;
	if (gensch) delete gensch;

	return (Schedule3);
}
//	-----------------------------------------------------------------
//		returns the roll date for the specified convention
ARM_Date DateRoll(const ARM_Date& asof, qCDS_ADJ adj) 
{
	static ARM_Date refDate_CDSINDZ5(20,12,2005); 
	static ARM_Date refDate_CDSINDM5(20,06,2005); 
	static ARM_Date refDate_CDSINDH5(20,03,2005); 
	static ARM_Date refDate_CDSINDM6(20,06,2006); 
	static ARM_Date refDate_CDSINDZ6(20,12,2006); 
	static ARM_Date refDate_CDSINDM7(20,06,2007); 
	static ARM_Date refDate_CDSINDZ7(20,12,2007); 
	static ARM_Date refDate_CDSINDM8(20,06,2008); 
	static ARM_Date refDate_CDSINDZ8(20,12,2008); 
	static ARM_Date refDate_CDSINDM4(20,06,2004); 
	static ARM_Date refDate_CDSINDU4(20,9,2004); 
	static ARM_Date refDate_CDSINDZ4(20,12,2004); 

	switch(adj)
	{
		case qCredit_CDSDTRX:	
		case qCredit_CDSINDH5:	return refDate_CDSINDH5; 
		case qCredit_CDSDIND:	
		case qCredit_CDSINDM5:	return refDate_CDSINDM5; 
		case qCredit_CDSINDZ:
		case qCredit_CDSINDZ5:	return refDate_CDSINDZ5; 
		case qCredit_CDSINDM6:	return refDate_CDSINDM6; 
		case qCredit_CDSINDZ6:	return refDate_CDSINDZ6 ;
		case qCredit_CDSINDM7:	return refDate_CDSINDM7; 
		case qCredit_CDSINDZ7:	return refDate_CDSINDZ7 ;
		case qCredit_CDSINDM8:	return refDate_CDSINDM8; 
		case qCredit_CDSINDZ8:	return refDate_CDSINDZ8 ;
		case qCredit_CDSINDU4:	return refDate_CDSINDU4 ;
		case qCredit_CDSINDM4:	return refDate_CDSINDM4 ;
		case qCredit_CDSINDZ4:	return refDate_CDSINDZ4 ;
		case qCredit_Adjust20:
		/*case qCredit_Special_None_Date:
		case qCredit_Special_None_YF:*/
			{
				if		(asof.GetMonth()< 3 || (asof.GetMonth()== 3 && asof.GetDay()<20)) return ARM_Date(20,3,asof.GetYear()); 
				else if (asof.GetMonth()< 6 || (asof.GetMonth()== 6 && asof.GetDay()<20)) return ARM_Date(20,6,asof.GetYear()); 
				else if (asof.GetMonth()< 9 || (asof.GetMonth()== 9 && asof.GetDay()<20)) return ARM_Date(20,9,asof.GetYear()); 
				else if (asof.GetMonth()<12 || (asof.GetMonth()==12 && asof.GetDay()<20)) return ARM_Date(20,12,asof.GetYear()); 
				else return  ARM_Date(20,3,asof.GetYear()+1); 
			} ;
		case qCredit_Adjust20N: 
			{
				if		(asof.GetMonth()< 3 || (asof.GetMonth()== 3 &&  asof.GetDay()<20)) return ARM_Date(20,6,asof.GetYear()); 
				else if (asof.GetMonth()< 6 || (asof.GetMonth()== 6 &&  asof.GetDay()<20)) return ARM_Date(20,9,asof.GetYear()); 
				else if (asof.GetMonth()< 9 || (asof.GetMonth()== 9 &&  asof.GetDay()<20)) return ARM_Date(20,12,asof.GetYear()); 
				else if (asof.GetMonth()<12 || (asof.GetMonth()==12 &&  asof.GetDay()<20)) return ARM_Date(20,3,asof.GetYear()+1); 
				else return  ARM_Date(20,6,asof.GetYear()+1); 
			} ;
		case qCredit_Adjust20SA:
			{
				if		(asof.GetMonth()<3 || (asof.GetMonth()==3 && asof.GetDay()<20)) return ARM_Date(20,3,asof.GetYear()); 
				else if (asof.GetMonth()<9 || (asof.GetMonth()==9 && asof.GetDay()<20)) return ARM_Date(20,9,asof.GetYear()); 
				else return  ARM_Date(20,3,asof.GetYear()+1); 
			}; 
		case qCredit_STDINDEX:
			{
				if		(asof.GetMonth()<3 || (asof.GetMonth()==3 &&  asof.GetDay()<20)) return ARM_Date(20,12,asof.GetYear()-1); 
				else if (asof.GetMonth()<9 || (asof.GetMonth()==9 &&  asof.GetDay()<20)) return ARM_Date(20,6,asof.GetYear()); 
				else return  ARM_Date(20,12,asof.GetYear()); 
			}; 

		case qCredit_Default :
			return asof ;
	} 
	ICMTHROW(ERR_INVALID_ARGUMENT, "Can't compute DateRoll for convention "<<adj); 
}


ARM_ZeroCurve*	GenerateIRCurveMovedInTime(ARM_ZeroCurve* IRCurve, ARM_Date TheDate)
{
	if (IRCurve == NULL)
		ICMTHROW(ERR_INVALID_ARGUMENT, "NULL IR CURVE: Can't move it!");

	ARM_ZeroCurve*	New_IRCurve;

	New_IRCurve	=	(ARM_ZeroCurve*)	IRCurve->Clone();

	// if one day, I want to use a Business Day Shift
//	ARM_Date NextBusDay = IRCurve->GetAsOfDate().NextBusinessDay(1) ;
			
	New_IRCurve->SetAsOfDate(TheDate) ;

	return	New_IRCurve;
}


// --------------------------------------------------------------------
// Quick CDO pricing with start date & end date (for single correlation)
// --------------------------------------------------------------------
void GenerateScheduleYF(ARM_Date& AsOf, ARM_Date& startdate, ARM_Date& enddate,long frequency,
						ARM_Vector& YFstartdates,ARM_Vector& YFenddates)
{
	YFstartdates.Resize(0);
	YFenddates.Resize(0);

	ARM_Vector* startdates=NULL;
	ARM_Vector* enddates=NULL;

	int i=0;

	startdates = CptStartDates(startdate,enddate,frequency,
                                            K_MOD_FOLLOWING, K_LONGSTART, K_MATUNADJUSTED,
                                            "EUR",0);

	enddates = new ARM_Vector(startdates);

	YFstartdates.Resize(startdates->size());
	YFenddates.Resize(startdates->size());

	for (i=0;i<enddates->size()-1;i++)
	{enddates->Elt(i)=enddates->Elt(i+1);}

	enddates->Elt(enddates->size()-1)=enddate.GetJulian();

	for (i=0;i<enddates->size();i++)
	{
		YFstartdates.Elt(i)=(startdates->Elt(i)-AsOf.GetJulian())/365.;
		if (YFstartdates.Elt(i)<0.) {YFstartdates.Elt(i)=0.;}
		YFenddates.Elt(i)=(enddates->Elt(i)-AsOf.GetJulian())/365.;
		if (YFenddates.Elt(i)<0.) {YFenddates.Elt(i)=0.;}
	}

	if (startdates) delete startdates;
	if (enddates) delete enddates;
}

double ComputeEL_CDO_LHP(const double& tranche_down, const double& tranche_up, const double& beta_down,
							   const double& beta_up,  const double& Pdef,  double recovery)
{
	double its_FH_Barrier=0.;

	if (fabs(Pdef) < DB_TOL) its_FH_Barrier = -10.;
	else if (fabs(Pdef-1.0) < DB_TOL) its_FH_Barrier = 10.;
	else its_FH_Barrier = NAG_deviates_normal_dist(Pdef);

	double ExpectedLoss = 0.,PrbLoss_down=0.,PrbLoss_up=0.;
	double K1=0.,K2=0.;
	double DevNorm1=0.,DevNorm2=0.;

	K1 = tranche_down/(1.-recovery);
	if (CHECK_EQUAL(K1, 0.0)) {DevNorm1 = _MINUS_INFINITY_;}
	else if (CHECK_EQUAL(K1, 1.0)) {DevNorm1 = _PLUS_INFINITY_;}
	else DevNorm1 = NAG_deviates_normal(K1);

	K2 = tranche_up/(1.-recovery);
	if (CHECK_EQUAL(K2, 0.0)) {DevNorm2 = _MINUS_INFINITY_;}
	else if (CHECK_EQUAL(K2, 1.0)) {DevNorm2 = _PLUS_INFINITY_;}
	else DevNorm2 = NAG_deviates_normal(K2);

	PrbLoss_up = NAG_bivariate_normal_dist(-DevNorm2,
							its_FH_Barrier,-sqrt(1.-beta_up*beta_up));

	PrbLoss_down = NAG_bivariate_normal_dist(-DevNorm1,
							its_FH_Barrier,-sqrt(1.-beta_down*beta_down));

	ExpectedLoss=(PrbLoss_down-PrbLoss_up)/(K2-K1);

	return (ExpectedLoss);
}

double ComputeEL_LHP(const double& x, const double& beta,  const double& Pdef,  double recovery)
{
	double its_FH_Barrier=0.,X=0.;

	if (fabs(Pdef) < DB_TOL) its_FH_Barrier = -10.;
	else if (fabs(Pdef-1.0) < DB_TOL) its_FH_Barrier = 10.;
	else its_FH_Barrier = NAG_deviates_normal_dist(Pdef);

	if (fabs(x) < DB_TOL) X = -10.;
	else if (fabs(x-1.0) < DB_TOL) X = 10.;
	else X = NAG_deviates_normal_dist(x);

	double ExpectedLoss = NAG_cumul_normal((1./beta)*(sqrt(1.-beta*beta)*X-its_FH_Barrier));

	return (ExpectedLoss);
}

void ExtractVectorBeetwen2values(ARM_Vector* v,double start,double end, vector<double>& vout,int paytiming)
{
	vout.clear();

	if (paytiming==K_ARREARS){
	for (int i=0; i<v->GetSize();i++)
	{
		if ((v->Elt(i)<=end) && (v->Elt(i)>start))
		{vout.push_back(v->Elt(i));}
		else if (CHECK_EQUAL(v->Elt(i),end))
		{vout.push_back(v->Elt(i));}
	}} else
	{
	for (int i=0; i<v->GetSize();i++)
	{
		if ((v->Elt(i)<end) && (v->Elt(i)>=start))
		{vout.push_back(v->Elt(i));}
		else if (CHECK_EQUAL(v->Elt(i),start))
		{vout.push_back(v->Elt(i));}
	}}
}

int getDefLegInterestRule(int feelegIntRule)
{
	switch (feelegIntRule)
	{
	case K_ADJUSTED: return K_MATUNADJUSTED ;
	case K_UNADJUSTED: return K_UNADJUSTED ;
	case K_MATUNADJUSTED: return K_MATUNADJUSTED ;
	}
	ICMTHROW(ERR_INVALID_ARGUMENT,"getDefLegInterestRule: can't handle rule "<<feelegIntRule);  
}
void CptPeriodDates(const double& start,const double& end,
					const int& resetgap,const int& resettiming,const int& paytiming,ARM_Currency* ccy,
					double yfmaturity,qCDS_ADJ fwdadj,int paygap,double ForcePayDate,
					int resetweekday, int resetoccur,
					double& paydate,double& reset,double& fwdstart,double& fwdend)
{
	ARM_Date D_reset;

	if (resettiming==K_ADVANCE)
		D_reset = (ARM_Date) start;
	else
		D_reset = (ARM_Date) end;

	//determinationde la reset date
	if (resetoccur!=0)
	{
		int resetDay = D_reset.GetDayOfWeek() ;
		if (resetoccur>=0) 
			reset = D_reset.GetJulian() + MODULO((resetweekday-resetDay-1),7)+1+(resetoccur-1)*7 ; 
		else 
			reset = D_reset.GetJulian() + MODULO((resetweekday-resetDay),7)+resetoccur*7 ; 
		reset = ARM_Date(reset).GapBusinessDay(0,ccy->GetCcyName()).GetJulian(); 
	}else
	{
		reset = D_reset.GetJulian();
	}

	char tmp[10];
	itoa((int)yfmaturity,tmp,10);

	string IndexMty = tmp;
	IndexMty = IndexMty + "Y";

	int fwdGap=0 ; 
	if( strcmp( ccy->GetCcyName(), "AUD") == 0 )fwdGap = 0;
	else if( strcmp( ccy->GetCcyName(), "CAD") == 0 )fwdGap= 2;
	else fwdGap= ccy->GetSpotDays();

	fwdstart = ARM_Date(reset).NextBusinessDay(fwdGap,ccy->GetCcyName()).GetJulian(); 			
	fwdend=AddPeriod(ARM_Date(reset),IndexMty,ccy->GetCcyName(),false,fwdadj).GetJulian() ;

	if (ForcePayDate!=0)
	{
		ARM_Date D_Paydate=(ARM_Date)start; //K_ADVANCE

		if (paytiming==K_ARREARS)
			{D_Paydate=(ARM_Date) end;}

		D_Paydate.NextBusinessDay(paygap,ccy->GetCcyName());
		paydate =D_Paydate.GetJulian(); 
	}
	else
	{
		paydate=ForcePayDate;
	}
}

void CptResetDateWeekly(const double& start,const double& end,int resetweekday,double yfmaturity, ARM_Currency* ccy, qCDS_ADJ fwdadj,
							  ARM_Vector& VresetDates, ARM_Vector& Vfwdstart, ARM_Vector& Vfwdend)
{

	ARM_Date D_reset;
	D_reset = (ARM_Date) start;
	int resetDay = D_reset.GetDayOfWeek() ;
	bool ended = false;
	// first reset
	int between = resetweekday - resetDay;
	if (between <  0) {
		// case + 1 w
		between = 7-resetDay + resetweekday ;
		D_reset.AddDays(between);
	}else {
		D_reset.AddDays(between);
	}
	VresetDates.push_back(D_reset.GetJulian());
	// all reset
	while (!ended){
		D_reset.AddDays(7);	
		if ( D_reset.GetJulian()  >end) {
			ended = true;
			continue;
		}
		VresetDates.push_back(D_reset.GetJulian());
	}
	// test Non business day
	for (int i=0; i< VresetDates.size(); i++) {
		VresetDates[i] = ARM_Date(VresetDates[i]).GapBusinessDay(0,ccy->GetCcyName()).GetJulian(); 
	}
	Vfwdstart.Resize(VresetDates.size());
	Vfwdend.Resize(VresetDates.size());
	// fwd start && end
		char tmp[10];
		double fwdstart =0;
		double fwdend = 0;
		itoa((int)yfmaturity,tmp,10);
		string IndexMty = tmp;
		IndexMty = IndexMty + "Y";
		int fwdGap=0 ; 
		if( strcmp( ccy->GetCcyName(), "AUD") == 0 )fwdGap = 0;
		else if( strcmp( ccy->GetCcyName(), "CAD") == 0 )fwdGap= 2;
		else fwdGap= ccy->GetSpotDays();

		for (int j=0; j< VresetDates.size(); j++) {
			ARM_Date dateRes = ARM_Date(VresetDates[j]);
			Vfwdstart[j] = (dateRes.NextBusinessDay(fwdGap,ccy->GetCcyName())).GetJulian(); 			
			Vfwdend[j]=AddPeriod(dateRes,IndexMty,ccy->GetCcyName(),false,fwdadj).GetJulian() ;
		}
		
}


double FlatVectorInterpol(const ARM_Vector& X,const ARM_Vector& Y, const double& x, bool valueinf)
{
	if (X.size()!=Y.size()) ICMTHROW(ERR_INVALID_ARGUMENT,"LinearVectorInterpol:vector have different sizes");
	if (X.size()==0) ICMTHROW(ERR_INVALID_ARGUMENT,"LinearVectorInterpol:empty vector");
	if (X.size()==1) return Y[0];
	int il=0,inf=CREDIT_DEFAULT_VALUE,sup=CREDIT_DEFAULT_VALUE; 

	//recherche de l'intervalle ]X[inf],X[sup][ contenant x
	for (il=0; il<X.size(); il++)
	{ if CHECK_EQUAL(X[il],x) return Y[il];
	  if (X[il]<x) inf=il; else break;}

	int inf2 = 0;
	if (inf>=0) inf2= inf; 

	for (il=inf2; il<X.size(); il++)
	{if (X[il]> x) {sup=il;break;}}

	if (inf==CREDIT_DEFAULT_VALUE)
	{inf = 0;
	 sup = inf+1;}		
	else if (sup==CREDIT_DEFAULT_VALUE)
	{ sup = X.size()-1;
	  inf = sup-1;}	

	if (valueinf) 
	{ return Y[inf];}
	else
	{ return Y[sup];}
}


void DeduceYearTermsForTSR(ICM_Correlation* correlation, 
						   double yf,
						   double& yt1_corr,
						   double& yt2_corr,
						   double& yt1_pdef,
						   double& yt2_pdef)
{
	yt1_corr=0.;yt2_corr=0.;yt1_pdef=0.;yt2_pdef=0.;

	ARM_Vector BaseCorrelSchedule; 
	correlation->GetCorrelationTerms(BaseCorrelSchedule);
	if (yf<=BaseCorrelSchedule.Elt(0))
	{yt1_corr=BaseCorrelSchedule.Elt(0);}
	else{
	yt1_corr = FlatVectorInterpol(BaseCorrelSchedule,BaseCorrelSchedule,yf,true);
	yt2_corr = FlatVectorInterpol(BaseCorrelSchedule,BaseCorrelSchedule,yf,false);}

	if (yt1_corr==yt2_corr)
	{yt1_corr = FlatVectorInterpol(BaseCorrelSchedule,BaseCorrelSchedule,yf-0.01,true);}

	if (yf<=BaseCorrelSchedule.Elt(0))
	{yt1_pdef = yf;}
	else{
	yt1_pdef = FlatVectorInterpol(BaseCorrelSchedule,BaseCorrelSchedule,yf,true);
	yt2_pdef = yf;}

	if (yt1_pdef==yt2_pdef && (yf))
	{yt1_pdef = FlatVectorInterpol(BaseCorrelSchedule,BaseCorrelSchedule,yf-0.01,true);}
}


std::string GetIndexName(const ARM_IRIndex&index)
{ 
	const ICM_Credit_Index* creditIndex =
		dynamic_cast<const ICM_Credit_Index*>(&index); 
	
	if (creditIndex) return creditIndex->GetLabel(0); 
	return ARM_ParamView::GetMappingName(S_INDEX_TYPES, index.GetIndexType()) ;
}

void vc2stdv(int size,double* V,vector<double>& out)
{
	out.clear();
	for (int i=0;i<size;i++) {out.push_back(V[i]);}
}


// --------------------------------------------------------------------
// Quick cds pricing 
// --------------------------------------------------------------------
double QuickCDSpv(ARM_Vector& schedule, const int& periodmatu,  double rate,double recovery,
					  const vector<double>& el,double notfeeleg,double notdefleg, 
					  ARM_ZeroCurve* zc, double& feepv,double& defpv)
{
	int nbperiods = MIN(schedule.GetSize()-1,periodmatu);

	feepv = defpv = 0.;

	for (int i=0;i<nbperiods;i++)
	{
		double begin = schedule.Elt(i);
		double end = schedule.Elt(i+1);

		if (end<=0.)
		{	continue; }
		else if ((end>0.) && (begin<0.))
		{	begin=0.; }

		double disc = zc->DiscountPrice(end);

		double deltaT = (365./360.)*(end-begin);
		feepv += deltaT*rate*notfeeleg*(1.-el[i+1])*disc;
		defpv += (1.-recovery)*notdefleg*(-el[i]+el[i+1])*disc;
	}

	return (feepv-defpv);

}

double QuickFwdSpread(ARM_Vector& curveschedule,
					  ARM_Vector& defaultproba,
					  double yfstart,
					  double yfend,
					  ARM_ZeroCurve* zc,
					  double recovery,
					  double step)
{
	int i=0;
	int indexstart=0.;

	vector<double> X;X.resize(curveschedule.size());
	vector<double> Y;Y.resize(curveschedule.size());

	ARM_Vector cdssched;
	vector<double> Pdef;
	
	for (i=0; i<curveschedule.size(); i++)
	{
		X[i] = curveschedule.Elt(i);
		Y[i] = defaultproba.Elt(i);
	}

	if (!CHECK_NULL(Y[0])) 
		{return CREDIT_DEFAULT_VALUE;}

	indexstart=0;
	double yf = 0;
	double pdef = 0.;
	for (i=0,yf=0.; yf<yfstart; i++)
	{
		pdef = LinearVectorInterpol(X,Y,yf);
		cdssched.push_back(yf);
		Pdef.push_back(pdef);
		yf += step;
	}

	indexstart = i;

	for (i=0,yf=yfstart; i<curveschedule.size() && (yf<=yfend); i++)
	{
		pdef = LinearVectorInterpol(X,Y,yf);
		cdssched.push_back(yf);
		Pdef.push_back(pdef);
		yf += step;
	}

	if (!CHECK_EQUAL(cdssched[cdssched.size()-1],yfend))
	{
		pdef = LinearVectorInterpol(X,Y,yfend);
		cdssched.push_back(yfend);
		Pdef.push_back(pdef);
	}

	double be1=0,duration1=0,feeleg1=0,defleg1=0,pv1=0.;
	be1 = duration1 = feeleg1 = defleg1 = 0.;

	if (!CHECK_NULL(yfstart))
	{
	pv1 = QuickCDSpv(cdssched,indexstart,1.,recovery,Pdef,1.,1., 
					  zc,feeleg1,defleg1);
	be1=defleg1/feeleg1;
	}

	double be2,duration2,feeleg2,defleg2;
	be2 = duration2 = feeleg2 = defleg2 = 0.;
	double pv2 = QuickCDSpv(cdssched,Pdef.size(),1.,recovery,Pdef,1.,1., 
					  zc,feeleg2,defleg2);
	be2=defleg2/feeleg2;

	double fwd = (be2*feeleg2 - be1*feeleg1)/(feeleg2 - feeleg1);

	if (feeleg1==0.)
	{fwd = be2;};

	return (fwd);

}

int FindIndexInVector(ARM_Vector& V,double value)
{
	vector<double> X;X.resize(V.size());
	vector<double> Y;Y.resize(V.size());

	for (int i=0;i<V.size();i++)
	{	
		X[i]=V[i];
		Y[i]=(double)i;
	}
	
	double val = FlatVectorInterpol(X,Y,value,true);

	return (int)val;
}

double QuickCdsDuration(ARM_Vector& curveschedule,
					  ARM_Vector& defaultproba,
					  double yfstart,
					  double yfend,
					  ARM_ZeroCurve* zc,
					  double recovery,
					  double step)
{
	int i=0;
	int indexstart=0.;

	vector<double> X;X.resize(curveschedule.size());
	vector<double> Y;Y.resize(curveschedule.size());

	ARM_Vector cdssched;
	vector<double> Pdef;
	
	for (i=0; i<curveschedule.size(); i++)
	{
		X[i] = curveschedule.Elt(i);
		Y[i] = defaultproba.Elt(i);
	}

	if (!CHECK_NULL(Y[0])) 
		{return CREDIT_DEFAULT_VALUE;}

	indexstart=0;
	double yf = 0;
	double pdef = 0.;
	for (i=0,yf=0.; yf<yfstart; i++)
	{
		pdef = LinearVectorInterpol(X,Y,yf);
		cdssched.push_back(yf);
		Pdef.push_back(pdef);
		yf += step;
	}

	indexstart = i;

	for (i=0,yf=yfstart; i<curveschedule.size() && (yf<=yfend); i++)
	{
		pdef = LinearVectorInterpol(X,Y,yf);
		cdssched.push_back(yf);
		Pdef.push_back(pdef);
		yf += step;
	}

	if (!CHECK_EQUAL(cdssched[cdssched.size()-1],yfend))
	{
		pdef = LinearVectorInterpol(X,Y,yfend);
		cdssched.push_back(yfend);
		Pdef.push_back(pdef);
	}

	double be1=0,duration1=0,feeleg1=0,defleg1=0,pv1=0.;
	be1 = duration1 = feeleg1 = defleg1 = 0.;

	if (!CHECK_NULL(yfstart))
	{
	pv1 = QuickCDSpv(cdssched,indexstart,1.,recovery,Pdef,1.,1., 
					  zc,feeleg1,defleg1);
	be1=defleg1/feeleg1;
	}

	double be2,duration2,feeleg2,defleg2;
	be2 = duration2 = feeleg2 = defleg2 = 0.;
	double pv2 = QuickCDSpv(cdssched,Pdef.size(),1.,recovery,Pdef,1.,1., 
					  zc,feeleg2,defleg2);
	be2=defleg2/feeleg2;

	double duration = feeleg2;

	return (feeleg2);

}

