#ifndef STRUCTS_H
#define STRUCTS_H

//----------------------------------------
#include <stdio.h>
#include <math.h>
#include<sstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <string>
#include "assert.h"
#include <conio.h>
#include <time.h>
#include "vector"
#include <algorithm>

using namespace std;  //introduces namespace std


//#define bitSetNumOfMeta 3//16 //16//65
//#define bitSetNumOfReac 6//24 //24//95

#define NumberOfMETABOLITEs		12//31	//1655	//S6::4  S5::4  S4::3  S3::3  S2::6  S1::6	 S::3			sample :: 16	cho::12		revTest_1000i :: 1000
#define NumberOfREACTIONSs		18//27	//2381	//S6::5  S5::5  S4::4  S3::5  S2::6  S1::7	 S::6			sample :: 24	cho::18		revTest_1000i :: 1001
#define NumberOfREACTIONSsPlus	19//28	//2382	//S6::6  S5::6  S4::5  S3::6  S2::7  S1::8	 S::7			sample :: 25	cho::19		revTest_1000i :: 1002

#define NumberOfCandidates	48//48

#define MaxNumOfRecInOut	4	//4
#define	MaxNumOfMetInOut	4	//4

#define MaxDepth			4	//4

#define	STACKSIZE			1000
#define LOOPLimit			10

int DONE = 0;

//----------------------------------------- Structs
struct FLUX {
	int		reationName;
	float	reactionFlux;
};
//...............................
struct INPUT {
	vector<int> metabolitNamesIn;
	vector<int>	metabolitNamesOut;			
	int		reacNum;
	float	Coef;
	bool	direction;

	int		type;  // added for passOne --> which type of in/out is in a path 

	void setParams(vector<int> m, vector<int> s, int r, float c, bool d, int t/*, bitset<bitSetNumOfMeta> sp, bitset<bitSetNumOfMeta*3> sp3*/)
	{	
		metabolitNamesIn	= m;
		metabolitNamesOut	= s; 
		reacNum = r; 
		Coef = c;
		direction = d; 
		type = t;
		/*SetPoints = sp;
		SetPoints3 = sp3; */
	}

	void setType(int t)
	{
		type = t; 
	}

};
//...............................
struct OUTPUT {
	vector<int> metabolitNamesIn;
	vector<int>	metabolitNamesOut;
	int		reacNum;
	float	Coef;
	bool	direction;

	int				type; // added for passOne --> which type of in/out is in a path 

	void setParams(vector<int> m, vector<int> s, int r, float c, bool d, int t/*, bitset<bitSetNumOfMeta> sp, bitset<bitSetNumOfMeta*3> sp3*/)
	{
		metabolitNamesIn	= m;
		metabolitNamesOut	= s; 
		reacNum		= r; 
		Coef		= c;
		direction	= d; 
		type		= t;
		/*SetPoints = sp;
		SetPoints3 = sp3;*/
	}

	void setType(int t)
	{
		type = t; 
	};

};
//...............................
struct NODE {
	vector<INPUT>	inputs;
	vector<OUTPUT>	outputs;
	//	vector<PATH>	paths;
	int				NodeName;

	bool			FWorBK; // for passTwo ==> TRUE: when we want to go through this node again for finding flux values 

	void Node(vector<INPUT> i, vector<OUTPUT> o, int nn, bool FB)
	{
		inputs		= i;
		outputs		= o; 
		NodeName	= nn; 
		FWorBK		= FB;
	}

	void setFWorBK(bool flag)
	{
		FWorBK = flag; 
	}

	bool getFWorBK(void)
	{
		return FWorBK; 
	}

};

//...............................
//...............................
typedef struct INPUT_ARRAY {
	int metabolitNamesIn [MaxNumOfMetInOut];
	int	metabolitNamesOut [MaxNumOfMetInOut];	

	int		reacNum;
	float	Coef;

	int		RECstatus;

	int		notGoodPrimaryCandida;  // 1:NOKay  0:OKay
	int		reactionNameWeGotBackTo; 

	//INPUT_ARRAY();//: metabolitNamesIn(&), metabolitNamesOut(&) {} ;

	//bool	direction;
	//int		type;  // added for passOne --> which type of in/out is in a path 
}INPUT_ARRAY;//INPUTS,*inputs;
//...............................
typedef struct OUTPUT_ARRAY {
	int metabolitNamesIn [MaxNumOfMetInOut];
	//size_t metabolitNamesIn_size;
	int	metabolitNamesOut [MaxNumOfMetInOut];	

	int		reacNum;
	float	Coef;

	int		RECstatus;

	int		notGoodPrimaryCandida;  // 1:NOKay  0:OKay
	int		reactionNameWeGotBackTo; 

	//OUTPUT_ARRAY();//: metabolitNamesIn(&), metabolitNamesOut(&) {}

	//bool	direction;
	//int		type;  // added for passOne --> which type of in/out is in a path 
}OUTPUT_ARRAY;//OUTPUTS, *outputs;
//...............................
typedef struct NODE_ARRAY {
	INPUT_ARRAY		inputs	[MaxNumOfRecInOut];
	OUTPUT_ARRAY	outputs [MaxNumOfRecInOut];
	//size_t INPUT_ARRAY_SIZE;
	//	vector<PATH>	paths;
	//int				NodeName;

	//bool			FWorBK; // for passTwo ==> TRUE: when we want to go through this node again for finding flux values 
	int METstatus;
	
//	int reacCover [NumberOfREACTIONSs];

	//NODE_ARRAY():outputs(new OUTPUT_ARRAY) {}
	
	//NODE_ARRAY(): outputs(&OUTPUTS), inputs(&INPUTS) {}
	
	//inputs(&INPUTS){}			
//	NODE_ARRAY():pINPUT_ARRAY(new INPUT_ARRAY) {}
}NODE_ARRAY;

//...............................
typedef struct NET {

	NODE_ARRAY net [NumberOfMETABOLITEs];

}NET;

//...............................
typedef struct oneEFM {

	//int recName [NumberOfREACTIONSs+1];

//	int		statusM					[NumberOfMETABOLITEs];
//	int		statusMetRecIn			[NumberOfMETABOLITEs]; 
//	int		statusMetRecOt			[NumberOfMETABOLITEs];

//	int		count					[NumberOfMETABOLITEs];
	//int		statusMetRecIn_BW		[NumberOfMETABOLITEs]; 
	//int		statusMetRecOt_BW		[NumberOfMETABOLITEs]; 

//	int		statusR					[NumberOfREACTIONSs];		//just the first bit is used 				
//	float	recFlux					[NumberOfREACTIONSs+1];

//	int		FREEZE;	

	//#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#**#

	int		metVisited			[NumberOfMETABOLITEs];	// bool -- set firstTime
//	int		dirChanged			[NumberOfMETABOLITEs];	// bool -- wentFromOutput changToInput ViseVersa

	int		loop				[NumberOfMETABOLITEs];	// number of the times you pass a node 

	int		primaryInput		[NumberOfMETABOLITEs];	// keep the index
	int		primaryOutput		[NumberOfMETABOLITEs];  // keep the index

	int		AllStable			[NumberOfMETABOLITEs];	// 0 or 1
	int		AllStableCount;

	//int		statusMetRecIn			[NumberOfMETABOLITEs];	//if for each non_primary reaction reached MAX(3) :: if (!dirChanged) change direction else nonEFM::TRUE
	//int		statusMetRecOt			[NumberOfMETABOLITEs];

	float	recFlux				[NumberOfREACTIONSs+1];	//lastIndex ==> 0:NotDone 1:DoneEFM 2:DoneNonEFM    if AllStable=0 AllStableCount++  (AllStableCount=TARGET) ==> DoneEFM
		
	int		isUpdate			[NumberOfREACTIONSs];	// for each reaction set 1 when u update the flux, reset to 0 when u 			

	// nonEFM situations ::: 
	// set a count on function changeFlux ==> after a TARGET set as nonEFM (for avoiding loops)

}oneEFM;
//---------------------------------------- Global Variables 

vector<char> IfBeingReversible;
int NumbersOfMetabolits=0, NumbersOfReactions=0;

oneEFM EFMs [NumberOfCandidates]; 

NET network [1];				//just to keep it like EFMs

//NODE_ARRAY* network [NumberOfMETABOLITEs-1];

vector<NODE> net; 
int sizeCount_glob = 0;

//vector<vector<FLUX>> EFMs; 
//---------------------------------------- Functions 

 vector<vector <float>> readFileS (istream & fin); 
 vector<vector<int>> FindSource(vector<vector<float>>& S, int j, float temp);
 vector<NODE> CreateNetwork(vector<vector<float>> S/*, vector<NODE>& network*/);
/* NODE_ARRAY **/ void CreateNetwork_array(vector<vector<float>> S/*, vector<NODE>& network*/);
 void printEFMs(ofstream &fout);
 /*NODE_ARRAY*/void convertNODEtoNODEArray(NODE node, int num);
//----------------------------------------

#endif /*STRUCTS_H*/