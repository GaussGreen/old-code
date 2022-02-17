#if !defined(AFX_COMMON_H__6A7521CB_E62B_11D2_9FE5_080009C16B4A__INCLUDED_)
#define AFX_COMMON_H__6A7521CB_E62B_11D2_9FE5_080009C16B4A__INCLUDED_

// Structure and Type Definitions

const int CHF_TRANCHE_ID_LEN = 30;
const int CHF_MAX_FORWARD_RATES = 480;

#define PP_MAXTERM				13
#define PP_MAXINCENTIVE			37
#define PP_MAX_ENUM_LOANTYPE	64
#define PP_MAX_CPR				99.9

//////////////////////////////////////////////////////////////////////////////
/*This is a list of valid products types (handeled by the prepayment model) and their names
Example: CF30Yr (Conventional Fixed 30Yr), P771 (7/1 Arm), P751 (5/1 Arm), etc

Each Tranche or bucket of loans will be assigned to one of these products.
*/

typedef struct
{
	CString product_type_str;
	int		product_type;
} product_type_struct;


///////////////////////////////////////////////////////////////////
class ProductTypeArray : public CPtrArray
{
public:
	ProductTypeArray() : CPtrArray() {};
	virtual	~ProductTypeArray()
	{
		ClearAll();
	};
	product_type_struct * GetAt(int i) {return (product_type_struct *)CPtrArray::GetAt(i);};
	void ClearAll()
	{
		for (int i = 0; i < GetSize(); i++)
			if(GetAt(i))
				delete GetAt(i);
		CPtrArray::RemoveAll();
	};

	void add(int id, PSTR szProduct)
	{
		product_type_struct *pProductTypeStruct;
		pProductTypeStruct = new product_type_struct;
		pProductTypeStruct->product_type = id;
		pProductTypeStruct->product_type_str = szProduct;
		CPtrArray::Add((void *)pProductTypeStruct);
	}

};

////////////////////////////////////////////////////////////////////////

enum e_CHF_error_code
{
	success, 
	general_error, 
	param_value_error, 
	database_error, 
	memory_error, 
	math_error
};

enum e_CHF_valuation_type
{
	servicing,
	whole_loan
};

/* This is the RATES Structure. This is the entire treasury curve, the swap curve and the mortgage rate
at any given simulation month.

At t = 0, this will be the sport interest rate cuve that is loaded into the valuation engine. 

An array of CHF_RATE_STRUCTURE with 480 elements is created in the DLL during initialization. A pointer to
the array is sent to the valuation engine on demand, when CHFGetForwardRatesPointer() is called.
*/
typedef struct tag_RATES
{
	double dSwap_1m;
	double dSwap_3m;
	double dSwap_6m;
	double dSwap_9m;
	double dSwap_12m;
	double dSwap_2y;
	double dSwap_3y;
	double dSwap_4y;
	double dSwap_5y;
	double dSwap_7y;
	double dSwap_10y;
	double dSwap_15y;
	double dSwap_20y;
	double dSwap_30y;
	double dTreasury_3m;
	double dTreasury_6m;
	double dTreasury_1y;
	double dTreasury_2y;
	double dTreasury_5y;
	double dTreasury_10y;
	double dTreasury_30y;
	double dMortgage_rate;
} CHF_RATE_STRUCTURE;

/*This structure defines some of the parameters that would change for each simulation month. Hence
this is passed to the CHFCalcPrepay() function.
*/

typedef struct tag_CHFPREPAYMENTMODELSTRUCT
{
	int			simulation_month;
	int			simulation_year;
	
	char			szTranche_id[CHF_TRANCHE_ID_LEN+1];
	int			loan_type;
	enum e_CHF_valuation_type	valuation_type;
	
	double		gross_wac;
	double		original_term;
	double		remain_term;
	double		wala;
	double		state_factor;
	double		remain_balance;
	long		loan_count;
	double		delinquency30;
	double		delinquency60;
	double		delinquency90;
	double		foreclosure;
	
	enum e_CHF_error_code	error_code;

} CHF_PREPAYMENTMODELSTRUCTURE;


// Exported Functions

/*****************************************************
PURPOSE:
Sets date and spot rates and signals to re-load prepay parameters.

INPUT PARAMETERS:
date;
pRateStruct - pointer to CHF_RATE_STRUCTURE.
path: Path to Portfolio Database (Path to PrepayDB)
reload: True to mean reload the prepay db
When to call: This is called whenever the portfolio is opened or loaded. 

Called for every tranche. During load and also during pricing.
RETURNS : True/False. 
*****************************************************/
BOOL WINAPI CHFSetInitParams(CHF_RATE_STRUCTURE *pRateStruct, COleDateTime date, LPSTR path, BOOL reload);

/*****************************************************
PURPOSE:
Gets pointer to forecast rates structure. Valuation Engine will populate this with rates. This is an array of CHF_RATE_STRUCTURE,
where the size of the array is given by CHF_MAX_FORWARD_RATES (480 for now).

When to call: This can be once whenever the portfolio is reloaded or opened.
However, the rates need to be populated once for each OAS path. 
For simplicity, it may be better to call this function once for each OAS path.

RETURNS : A valid pointer to array of CHF_RATE_STRUCTURE. 
*****************************************************/
CHF_RATE_STRUCTURE * WINAPI CHFGetForwardRatesPointer();

/*****************************************************
PURPOSE:
Called to setup prepayment parameters and constants. 
Pops up custom setup screen. 

NOTE: This modifies parameters for the entire portfolio.

INPUT PARAMETERS: None
*****************************************************/
void WINAPI CHFSetDefaults();

/*****************************************************
PURPOSE:
Fills array with available product types. 
This function can be used by the valuation engine to associate a product type to 
a tranche or bucket of loans. Valuation engine would store the product id (an integer)
as part of the tranche information.

INPUT PARAMETERS:
Empty ProductTypeArray.
 
RETURNS : TRUE - success. 
/*****************************************************/
ProductTypeArray* WINAPI CHFGetProductList(int &iNumProducts);


/*****************************************************
PURPOSE:
Calculates Prepayment rate CPR. To be called for each OAS path, for each tranche 
and for each simulation month. 

INPUT PARAMETERS:
pPrepayStruct - pointer to CHF_PREPAYMENTMODELSTRUCTURE.

RETURNS : CPR. 
*****************************************************/
double WINAPI CHFCalcPrepay(CHF_PREPAYMENTMODELSTRUCTURE *pPrepayStruct);

#endif // !defined(AFX_COMMON_H__6A7521CB_E62B_11D2_9FE5_080009C16B4A__INCLUDED_)
