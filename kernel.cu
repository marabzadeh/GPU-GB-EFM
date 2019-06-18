//----------------------	EFM Finder simplifized architecture 
//----------------------	Mona Arabzadeh
//----------------------	Dey 1395  
//----------------------	PhD -- part3

//auto complete ==> ctlr + space

//nvcc -deviceemu	

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
#include <vector>
#include <algorithm>
#include <bitset>

#include <iterator> 

#include "cuPrintf.cuh"
#include "cuPrintf.cu"

#include <cuda_profiler_api.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "cuda.h"
#include "curand.h"
#include "curand_kernel.h"

/*#include <thrust/device_vector.h>
#include <thrust/copy.h>*/

#include "structs.h"
#include "funcs.cpp"
#include "lock.h"


using namespace std;  //introduces namespace std

//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
// read network from matrix 
// --- statics
// network is a set of NODEs 
// NODEs are consist of INPUTs and OUTPUTs 
// NODE dynamics: seen, FWorBW

// --- dynamics 
// PATH : flux, history(meta,pathnum), direction,    
//-------------------------------------------------------
__device__ unsigned int RNG()
{   
    unsigned int m_w = 150;
    unsigned int m_z = 40;

	unsigned int res = 0;

   // for(int i=0; i < 100; i++)
    //{
        m_z = 36969 * (m_z & 65535) + (m_z >> 16);
        m_w = 18000 * (m_w & 65535) + (m_w >> 16);

 //       cout <<(m_z << 16) + m_w << endl;  /* 32-bit result */
    //}

	res = ((m_z << 16) + m_w);
	return res;
}
//-------------------------------------------------------
__device__ void convDecimalToBase(int decimalNum, int BN, int* BS, int BSSize)
{
	//int* b = new int[BN];
	//memset(b, 0, BN);

	int x = decimalNum;
	int y = BN; 
	int point = BSSize - 1;

	if (x<y)
	{
		BS[point] = x;
	}
	else{
		while(x >= y)
		{
			BS[point] = x%y;
			x /= y; 
			point--;
		}

		BS[point] = x;
	}

	return; 
}
//-------------------------------------------------------
__device__ int findReacOrder(int reactionName, INPUT_ARRAY* inp)
{
	int ord = -1;

	for (int i=0; i<MaxNumOfRecInOut; i++)
	{
		if(inp[i].reacNum == reactionName)
		{
			ord = i;
			break;
		}
	}

	return ord;
}
//-------------------------------------------------------
__global__ void	depthFUNC(NET *d_network)
{
	int iMETA	=	blockIdx.x;					// on metabolites
	int	jThread	=	threadIdx.x;	

	int  baseStr[MaxDepth] = {0}; 
	int  baseNum = MaxNumOfRecInOut;

	convDecimalToBase(jThread, baseNum, baseStr, MaxDepth);

	int isDone		=	0;		//set if reach level = 10 || see iMETA 

	int firstRec	= baseStr[0];

	int currentNODE		[STACKSIZE];
	int currentLevel	[STACKSIZE];

	int stackSize		= 0;			//keep the size of elements which are needed to analyze 
	int stackPointer	= 0;			//just go forward when an element added

	currentNODE[0]		= iMETA;
	currentLevel[0]		= 1;
	stackSize++;

	int recSizeOut		=   (d_network[0].net[iMETA].METstatus & 0x000000F0) >> 4;
	int recSizeIn		=   (d_network[0].net[iMETA].METstatus & 0x0000000F);
	int recSizeInOut	=	0;

	if ((recSizeOut == 1) && (recSizeIn ==1)) //if it has only one input and only one output
		recSizeInOut = 1;

	while (stackSize>0)
	{
		int popCurrentNode  = currentNODE[stackPointer];
		int popCurrentLevel = currentLevel[stackPointer];

		if (stackPointer == 99)
			return;
		//cuPrintf("here");

		stackPointer++;
		stackSize--;

		int currentREC	= baseStr[popCurrentLevel-1]; // (jThread%MaxNumOfRecInOut);
		
		int ifLastMet	=   (d_network[0].net[popCurrentNode].outputs[currentREC].RECstatus & 0x0000FF00) >> 8;
		recSizeOut		=   (d_network[0].net[popCurrentNode].METstatus & 0x000000F0) >> 4;


		if(ifLastMet>0)
		{ //----------------------------------------------------------------------------------last metabolite :: the output reaction has no metabolite  

		//if(jThread==8)
		//	cuPrintf("CuuRec = %d;recSizeOut = %d\n\r", currentREC,recSizeOut);

			if (currentREC < recSizeOut) //if valid
			{
			//	cuPrintf("CuuRec = %d;recSizeOut = %d\n\r", currentREC,recSizeOut);
				int metSize = ifLastMet;//(d_network[0].net[popCurrentNode].outputs[currentREC].RECstatus & 0x0000FF00) >> 8; //??? check
				for (int i=0; i</*MaxNumOfMetInOut*/metSize; i++)
				{
					if (d_network[0].net[popCurrentNode].outputs[currentREC].metabolitNamesOut[i] == iMETA)
					{
						//set what is needed to be set
						isDone = 1;
						int ord = findReacOrder(d_network[0].net[popCurrentNode].outputs[currentREC].reacNum, d_network[0].net[iMETA].inputs);

						if (recSizeInOut == 0) {
							d_network[0].net[iMETA].inputs[ord].notGoodPrimaryCandida			= 1;
							d_network[0].net[iMETA].inputs[ord].reactionNameWeGotBackTo			= d_network[0].net[iMETA].outputs[firstRec].reacNum;
							//cuPrintf("META = %d;thread = %d; %d\n\r", iMETA,jThread,d_network[0].net[iMETA].outputs[firstRec].reacNum);

							d_network[0].net[iMETA].outputs[firstRec].notGoodPrimaryCandida		= 1;
							d_network[0].net[iMETA].outputs[firstRec].reactionNameWeGotBackTo	= d_network[0].net[iMETA].inputs[ord].reacNum;
							//cuPrintf("META = %d;thread = %d; %d\n\n\r", iMETA,jThread,d_network[0].net[iMETA].inputs[ord].reacNum);
						}
					}//endif we are done
					else if (popCurrentLevel <= MaxDepth) // full stack while  (level <= MaxDepth)
					{
						currentNODE	[stackPointer+stackSize] = d_network[0].net[popCurrentNode].outputs[currentREC].metabolitNamesOut[i];
						currentLevel[stackPointer+stackSize] = popCurrentLevel+1;
						stackSize++;
					}// end else push back nodeName nd level to stack

				}//endfor i:MaxNumOfMetInOut

				if(isDone == 1)
					break; //break while  if Done

			}// end if
			else
				break; // break the while -- not a valid reaction

		}//-----------------------------------------------------------------------------end if las metabolite (no reaction)
	}//end while


	//cuPrintf("Block = %d;thread = %d\n\r", iMETA,jThread);

}// end depthFUNC
//-------------------------------------------------------
//-------------------------------------------------------
__global__ void MetaINIT(NET *d_network, oneEFM *d_EFMs)
{
	int j = blockIdx.x;   // on candidates
	int k = threadIdx.x;  // on reactions

	//for (int j=0; j<NumberOfCandidates; j++)
	//{
		d_EFMs[j].recFlux[NumberOfREACTIONSsPlus-1] = 0;  // as an initial value in initialize ==> 0 :: notDone

		// set the first reaction flux and trig output and inputs 
		int firstReaction	= 0;									//sample :: 15			S2::0		S2::4
		//int firstMetabolite = 0;									//sample :: 15			S2::2		S2::2
		//for (int k=0; k<NumberOfREACTIONSs; k++)
		//{
			d_EFMs[j].recFlux[k] = 0;
		//}

		d_EFMs[j].recFlux[firstReaction] = 1;

		/*firstReaction	= 3;	
		d_EFMs[j].recFlux[firstReaction] = 1;

		firstReaction	= 8;	
		d_EFMs[j].recFlux[firstReaction] = 1;*/

	//}//endForThread
}
//---------------------------------------------------------
__device__ void wait()
{
	for (int i=0; i<10000; i++){}	
}
//----------------------------------------------------------
//----------------------------------------------------------
//----------------------------------------------------------
__global__ void METAx(NET *d_network, oneEFM *d_EFMs, int *parVal) // input :: the whole network, index ==> i :: Executation of METAx on i
{
	
	if (threadIdx.x == 0)
		__shared__ oneEFM test;// = d_EFMs[blockIdx.x];	//for sharedMEmory test ==> LATER

	while(1)
	{
		//...............................................................................
		int iMET=-1,jBLOCK=-1; 
		float extraFL = 0;
		//...............................................................................

		iMET	= threadIdx.x;					// on metabolites
		jBLOCK	= blockIdx.x;					// on blocks
		//...............................................................................
		int numberOfInputs		=	 d_network[0].net[iMET].METstatus & 0x0000000F;
		int numberOfOutputs		=   (d_network[0].net[iMET].METstatus & 0x000000F0) >> 4;
		//...............................................................................
		//if(iMET == 7)
		//	cuPrintf("exFlux = %f\n\r",extraFL);

		 //float ffll = 0;
		 //float ccff = 0; 
		int activeInpCnt = 0, activeOutCnt = 0;

calAgain:for (int fin=0; fin < numberOfInputs; fin++)
		 {
			float fl = d_EFMs[jBLOCK].recFlux[d_network[0].net[iMET].inputs[fin].reacNum];
			if (fl != 0)
				activeInpCnt++;
			 extraFL += fl/*d_EFMs[jBLOCK].recFlux[d_network[0].net[iMET].inputs[fin].reacNum]*/ * d_network[0].net[iMET].inputs[fin].Coef;
			 d_EFMs[jBLOCK].isUpdate[d_network[0].net[iMET].inputs[fin].reacNum] = 0;
		 }
		 for (int fin=0; fin < numberOfOutputs; fin++)
		 {
			 float fl = d_EFMs[jBLOCK].recFlux[d_network[0].net[iMET].outputs[fin].reacNum];
			 if (fl != 0)
				activeOutCnt++;
			 extraFL += fl/*d_EFMs[jBLOCK].recFlux[d_network[0].net[iMET].outputs[fin].reacNum]*/ * d_network[0].net[iMET].outputs[fin].Coef;
			 d_EFMs[jBLOCK].isUpdate[d_network[0].net[iMET].outputs[fin].reacNum] = 0;
		 }

		 //................................................................................
		 if (extraFL == 0)		
		 {
			 d_EFMs[jBLOCK].AllStable[iMET] = 0;
		 }//................................/endelse extraFL=0/............................
		 else if (extraFL > 0)
		 {
			 /*cuPrintf("extraFL > 0 \r\n");
			 parVal[0]++;
			 cuPrintf("step=%d", parVal[0]);*/

			 //cuPrintf("iMET = %d..... extraFL > 0", iMET);

			 d_EFMs[jBLOCK].AllStable[iMET] = 1;
			 //if nonVisited :: checkVisited go from Primaries
			 if(d_EFMs[jBLOCK].metVisited[iMET] == 0)
			 {
				 int tId = blockIdx.x + threadIdx.x;

				 int primaryCandidateInp = -1;
				 for (int fin=0; fin < numberOfInputs; fin++)
				 { if (d_EFMs[jBLOCK].recFlux[d_network[0].net[iMET].inputs[fin].reacNum] != 0) {primaryCandidateInp = fin; break;} }

				 //-------------------------------- go randomly out 
				 //unsigned int rnd = RNG();
				 
				 //commented for test
				 int outputITis		=tId %numberOfOutputs;//tId*(numberOfOutputs-1)%numberOfOutputs;//tId*(numberOfOutputs-1)%numberOfOutputs;
//				 int outputITis = 0; 
				 //if (iMET == 17)
				 //	 outputITis = 1;
				 int outputITisName	= d_network[0].net[iMET].outputs[outputITis].reacNum;

				  //for debug
				 /*if (iMET == 0)
					 outputITisName = 1;
				 else if(iMET == 1)
					 outputITisName = 3;
				 else if (iMET == 2)
					 outputITisName = 4;
				 else if (iMET == 4)
					 outputITisName = 7;
				 else if (iMET == 5)
					 outputITisName = 8;
				 else if (iMET == 6)
					 outputITisName = 9;
				 else if (iMET == 7)
					 outputITisName = 10;
				 else if (iMET == 8)
					 outputITisName = 11;
				 else if (iMET == 10)
					 outputITisName = 8;

				 for (int q=0; q<numberOfOutputs; q++)
				 { if (outputITisName == (d_network[0].net[iMET].outputs[q].reacNum))	{outputITis = q; break; }  }*/
				 //for debug


				 float newFlux = abs(extraFL/d_network[0].net[iMET].outputs[outputITis].Coef); 
				 /*if (d_EFMs[jBLOCK].recFlux[outputITisName] != 0 ){
				 if (d_EFMs[jBLOCK].recFlux[outputITisName] !=  newFlux)
				 {d_EFMs[jBLOCK].recFlux[NumberOfREACTIONSsPlus-1] = 2;}	//-- check as nonEFM	:: fluxes dont match	
				 }
				 else*/ 
				 if (d_EFMs[jBLOCK].isUpdate[outputITisName] == 0){
					 d_EFMs[jBLOCK].recFlux[outputITisName] = newFlux;
					 //cuPrintf("recName = %d;recFlux = %f\n\r", outputITisName,newFlux);
					 d_EFMs[jBLOCK].isUpdate[outputITisName] = 1;
				 }
				 else
				 {
					  //cuPrintf("washere");
					 goto calAgain;
				 }
				 //wait();
				 //-------------------------------- set primary input and output (primaryInput is a challenge => check one of them)
				 d_EFMs[jBLOCK].primaryInput[iMET]  = primaryCandidateInp;
				 d_EFMs[jBLOCK].primaryOutput[iMET] = outputITis;
				 //--------------------------------
				 d_EFMs[jBLOCK].metVisited[iMET] = 1;
			 }//.........................................................end if nonVisited

			 else if (d_EFMs[jBLOCK].metVisited[iMET] == 1)                           //   ==0 or the same as primary reaction  :: nonEFM
			 {
				 if(((d_network[0].net[iMET].outputs[d_EFMs[jBLOCK].primaryOutput[iMET]].notGoodPrimaryCandida == 0) && (d_EFMs[jBLOCK].loop[iMET] < LOOPLimit)) || (activeInpCnt == 1)/*d_EFMs[jBLOCK].loop[iMET] == 0*/) //from primaries - out
				 {
					 int primryOutName = d_network[0].net[iMET].outputs[d_EFMs[jBLOCK].primaryOutput[iMET]].reacNum;
					 float oldFlux = d_EFMs[jBLOCK].recFlux[primryOutName];
					 float c		  = abs(d_network[0].net[iMET].outputs[d_EFMs[jBLOCK].primaryOutput[iMET]].Coef);
					 if (d_EFMs[jBLOCK].isUpdate[primryOutName]==0){
						 d_EFMs[jBLOCK].recFlux[primryOutName] = ((oldFlux*c) + extraFL)/c;
						 //cuPrintf("recName = %d;recFlux = %f\n\r", primryOutName, d_EFMs[jBLOCK].recFlux[primryOutName]);
						 d_EFMs[jBLOCK].isUpdate[primryOutName] = 1;
					 }
					 else 
					 {
						 //cuPrintf("washere");
						 goto calAgain;
					 }
					 //wait();
					 d_EFMs[jBLOCK].loop[iMET]++; 
				 }
				 else if((d_network[0].net[iMET].outputs[d_EFMs[jBLOCK].primaryOutput[iMET]].notGoodPrimaryCandida == 1) && (d_EFMs[jBLOCK].loop[iMET] < LOOPLimit)/*(d_EFMs[jBLOCK].loop[iMET] > 0) && (d_EFMs[jBLOCK].loop[iMET] <= 2)*/)	//change direction
				 {
					 int primryInpName = d_network[0].net[iMET].inputs[d_EFMs[jBLOCK].primaryInput[iMET]].reacNum;
					 float oldFlux = d_EFMs[jBLOCK].recFlux[primryInpName];
					 float c		  = abs(d_network[0].net[iMET].inputs[d_EFMs[jBLOCK].primaryInput[iMET]].Coef);

					 float newFlux = abs(((oldFlux*c) - extraFL)/c);
					 /*if ((newFlux == oldFlux) || (newFlux == 0))
					 d_EFMs[jBLOCK].recFlux[NumberOfREACTIONSsPlus-1] = 2;
					 else*/
					 if(d_EFMs[jBLOCK].isUpdate[primryInpName] ==0){
						 d_EFMs[jBLOCK].recFlux[primryInpName] = newFlux;
						  //cuPrintf("recName = %d;recFlux = %f\n\r", primryInpName, newFlux);
						 d_EFMs[jBLOCK].isUpdate[primryInpName] = 1;
					 }
					 else
					 {
						 //cuPrintf("washere");
						 goto calAgain;
					 }
					 //wait();	
					 d_EFMs[jBLOCK].loop[iMET]++;
				 }
				 else if(d_EFMs[jBLOCK].loop[iMET] >= LOOPLimit)	//direction has changed before -- check as nonEFM (loop detected)
				 {
					 d_EFMs[jBLOCK].recFlux[NumberOfREACTIONSsPlus-1] = 2;
				 }

			 }//.........................................................end else visited
			 /*for (int test=0; test<NumberOfREACTIONSs; test++)
			 {
				 cuPrintf("%f\t",d_EFMs[jBLOCK].recFlux[test]);
			 }
			 cuPrintf("\r\n");
			 */
		 }//................................/endelse extraFL>0/............................
		 else if (extraFL < 0)
		 {
			 /*cuPrintf("extraFL < 0\r\n");
			 parVal[0]++;
			 cuPrintf("step=%d", parVal[0]);*/

			 extraFL = -extraFL;
			 d_EFMs[jBLOCK].AllStable[iMET] = 1;
			 //if nonVisited :: checkVisited go from Primaries
			 if(d_EFMs[jBLOCK].metVisited[iMET] == 0)
			 {
				 int tId = blockIdx.x + threadIdx.x;

				 int primaryCandidateOut = -1;
				 for (int fin=0; fin < numberOfOutputs; fin++)
				 { if (d_EFMs[jBLOCK].recFlux[d_network[0].net[iMET].outputs[fin].reacNum] != 0) {primaryCandidateOut = fin; break;} }

				 //-------------------------------- go randomly in
				 //unsigned int rnd = RNG();
				 int inputITis		= tId%numberOfInputs;//tId*(numberOfInputs-1)%numberOfInputs;//tId*(numberOfInputs-1)%numberOfInputs;
				 int inputITisName	= d_network[0].net[iMET].outputs[inputITis].reacNum;


				 //for debug
				/* if (iMET == 0)
					 inputITisName = 0;
				 else if(iMET == 1)
					 inputITisName = 1;
				 else if (iMET == 2)
					 inputITisName = 3; //1or3
				 else if (iMET == 4)
					 inputITisName = 4;
				 else if (iMET == 5)
					 inputITisName = 7;
				 else if (iMET == 6)
					 inputITisName = 8;
				 else if (iMET == 7)
					 inputITisName = 9;
				 else if (iMET == 8)
					 inputITisName = 10;
				 else if (iMET == 10)
					 inputITisName = 11;

				 for (int q=0; q<numberOfInputs; q++)
				 { if (inputITisName == (d_network[0].net[iMET].inputs[q].reacNum))	{inputITis = q; break; }  }*/

				 //for debug


				 float newFlux = abs(extraFL/d_network[0].net[iMET].outputs[inputITis].Coef); 
				 /*if (d_EFMs[jBLOCK].recFlux[inputITisName] != 0 ){
				 if (d_EFMs[jBLOCK].recFlux[inputITisName] !=  newFlux)
				 {d_EFMs[jBLOCK].recFlux[NumberOfREACTIONSsPlus-1] = 2;}	//-- check as nonEFM	:: fluxes dont match	
				 }
				 else */
				 if (d_EFMs[jBLOCK].isUpdate[inputITisName]==0){
					 d_EFMs[jBLOCK].recFlux[inputITisName] = newFlux;
					 //cuPrintf("recName = %d;recFlux = %f\n\r", inputITisName, newFlux);
					 d_EFMs[jBLOCK].isUpdate[inputITisName] = 1;
				 }
				 else 
				 {
					  //cuPrintf("washere");
					 goto calAgain;
				 }
				 //wait();
				 //-------------------------------- set primary input and output (primaryOutput is a challenge => check one of them)
				 d_EFMs[jBLOCK].primaryInput[iMET]  = inputITis;
				 d_EFMs[jBLOCK].primaryOutput[iMET] = primaryCandidateOut;
				 //--------------------------------
				 d_EFMs[jBLOCK].metVisited[iMET] = 1;
			 }//.........................................................end if nonVisited

			 else if (d_EFMs[jBLOCK].metVisited[iMET] == 1)                           //   ==0 or the same as primary reaction  :: nonEFM
			 {
				 if(((d_network[0].net[iMET].inputs[d_EFMs[jBLOCK].primaryInput[iMET]].notGoodPrimaryCandida==0) && (d_EFMs[jBLOCK].loop[iMET] < LOOPLimit)) || (activeOutCnt == 1)/*d_EFMs[jBLOCK].loop[iMET] == 0*/)		//from primaries - in
				 {
					 int primryInpName = d_network[0].net[iMET].inputs[d_EFMs[jBLOCK].primaryInput[iMET]].reacNum;

					 float oldFlux = d_EFMs[jBLOCK].recFlux[primryInpName];
					 float c		  = abs(d_network[0].net[iMET].inputs[d_EFMs[jBLOCK].primaryInput[iMET]].Coef);
					 if (d_EFMs[jBLOCK].isUpdate[primryInpName] == 0)
					 {
						 d_EFMs[jBLOCK].recFlux[primryInpName] = ((oldFlux*c) + extraFL)/c;
						 //cuPrintf("recName = %d;recFlux = %f\n\r", primryInpName,  d_EFMs[jBLOCK].recFlux[primryInpName]);
						 d_EFMs[jBLOCK].isUpdate[primryInpName] = 1 ;
					 }
					 else{
						//cuPrintf("washere");
						goto calAgain;
					 }
					 //wait();
					 d_EFMs[jBLOCK].loop[iMET]++; 
				 }
				 else if((d_network[0].net[iMET].inputs[d_EFMs[jBLOCK].primaryInput[iMET]].notGoodPrimaryCandida==1) && (d_EFMs[jBLOCK].loop[iMET] < LOOPLimit))/*(d_EFMs[jBLOCK].loop[iMET] > 0) && (d_EFMs[jBLOCK].loop[iMET] <= 2)*/	//change direction
				 {
					 int primryOutName = d_network[0].net[iMET].outputs[d_EFMs[jBLOCK].primaryOutput[iMET]].reacNum;

					 float oldFlux = d_EFMs[jBLOCK].recFlux[primryOutName];
					 float c		  = abs(d_network[0].net[iMET].outputs[d_EFMs[jBLOCK].primaryOutput[iMET]].Coef);
					 float newFlux = abs(((oldFlux*c) - extraFL)/c);
					 /*if ((newFlux == oldFlux) || (newFlux == 0))
					 d_EFMs[jBLOCK].recFlux[NumberOfREACTIONSsPlus-1] = 2;
					 else*/
					 if (d_EFMs[jBLOCK].isUpdate[primryOutName] == 0)
					 {
						 d_EFMs[jBLOCK].recFlux[primryOutName] = newFlux;
						 //cuPrintf("recName = %d;recFlux = %f\n\r", primryOutName, newFlux);
						 d_EFMs[jBLOCK].isUpdate[primryOutName] = 1;
					 }
					 else 
					 {
						// cuPrintf("washere");
						 goto calAgain;
					 }
					 //wait();

					 d_EFMs[jBLOCK].loop[iMET]++;
				 }
				 else if(d_EFMs[jBLOCK].loop[iMET] >= LOOPLimit)	//direction has changed before -- check as nonEFM (loop detected)
				 {
					 d_EFMs[jBLOCK].recFlux[NumberOfREACTIONSsPlus-1] = 2;
				 }

			 }//.........................................................end else visited
			 /*for (int test=0; test<NumberOfREACTIONSs; test++)
			 {
				 cuPrintf("%f\t",d_EFMs[jBLOCK].recFlux[test]);
			 }
			  cuPrintf("\r\n");*/
		 }//................................/endelse extraFL<0/............................
		 //................................................................................

		 //................................................................................
		 //...........................//All Threads Done Check//...........................
		 int AllEFMsDone = 1;

		 int check = 1;
		 if ((d_EFMs[jBLOCK].recFlux[NumberOfREACTIONSsPlus-1] == 0) && (d_EFMs[jBLOCK].recFlux[NumberOfREACTIONSsPlus-1] != 2))
		 {
			 for (int m=0; m<NumberOfMETABOLITEs; m++)
			 {
				 if (d_EFMs[jBLOCK].AllStable[m] == 1)
					 check = 0;
			 }
		 }
		 if (check == 1){
			 d_EFMs[jBLOCK].AllStableCount++;
			 //d_EFMs[jBLOCK].recFlux[NumberOfREACTIONSsPlus-1] = 1;
		 }
		 if(d_EFMs[jBLOCK].AllStableCount==8)
			 d_EFMs[jBLOCK].recFlux[NumberOfREACTIONSsPlus-1] = 1;

		 //--------------
		 for (int j=0; j<NumberOfCandidates; j++)
		 {
			 if (d_EFMs[j].recFlux[NumberOfREACTIONSsPlus-1] == 0)
				 AllEFMsDone = 0;
		 }
		 //---
		 if (AllEFMsDone) // all threads done
			 break;		// break while ==> copyback EFMs to main Memory

		 //__syncthreads();
		 //................................................................................
	}//end while 1

} //end METAx
//------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------
cudaError_t EFM()
{
	cudaError_t cudaStatus;

	oneEFM		*d_EFMs;		//[NumberOfThreads];
	NET			*d_network;

	int			*d_parVal;
	int			parVal[1]	= {0};

	// define device variables  

	cudaDeviceReset();
	// Choose which GPU to run on, change this on a multi-GPU system.
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}
	//-----------------------------------------------------------------------------
	//=========================================================================== show memory usage of GPU
	size_t free_byte;
	size_t total_byte;
	cudaStatus = cudaMemGetInfo(&free_byte, &total_byte);

	if ( cudaSuccess != cudaStatus ){
		printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cudaStatus) );
		exit(1);
	}

	double free_db = (double)free_byte; 
	double total_db = (double)total_byte;
	double used_db = total_db - free_db;

	printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n",
		used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);
	//==========================================================================
	//-----------------------------------------------------------------------------  copy memory to gpu
		//---
	cudaStatus = cudaMalloc((void**)&d_parVal, sizeof(int));	
	if (cudaStatus != cudaSuccess) {
	fprintf(stderr, "cudaMalloc failed!");
	goto Error;
	}

	cudaStatus = cudaMemcpy(d_parVal, &parVal, sizeof(int), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	//----


	int sizeE = sizeof EFMs;									//for test -- passed the test
	cudaStatus = cudaMalloc((void**)&d_EFMs, sizeE);			//for test -- passed the test

	int sizeN = sizeof network;
	cudaStatus = cudaMalloc((void**)&d_network, sizeN);	
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}


	cudaStatus = cudaMemcpy(d_network, &network, sizeN, cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}

	/*cudaStatus = cudaMemcpy(d_DONE, &DONE, sizeof(int), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
	fprintf(stderr, "cudaMemcpy failed!");
	goto Error;
	}*/

	//=========================================================================== show memory usage of GPU
	//   size_t free_byte;
	//   size_t total_byte;
	cudaStatus = cudaMemGetInfo(&free_byte, &total_byte);

	if ( cudaSuccess != cudaStatus ){
		printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cudaStatus) );
		exit(1);
	}

	free_db = (double)free_byte; 
	total_db = (double)total_byte;
	used_db = total_db - free_db;

	printf("GPU memory usage: used = %f, free = %f MB, total = %f MB\n",
		used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);
	//==========================================================================
	//-----------------------------------------------------------------------------
	//-----------------------------------------------------------------------------
	//-----------------------------------------------------------------------------
	//-----------------------------------------------------------------------------
	//-------------------------calculate related reactions  
	cudaPrintfInit ();

	//dim3 dimGrid0(NumberOfMETABOLITEs);		//each block == one metbolite
	//float x = MaxNumOfRecInOut; float y = MaxDepth; 
	//int maxPathNumber = pow(x, y);
	//dim3 dimBlock0(maxPathNumber);						//for each metbolite create all poibilities on threads ==> how many threads we nead -- how to calc indexes
	
	//depthFUNC<<<dimGrid0, dimBlock0>>> (d_network);

	//-------------------------
	//-------------------------


	dim3 dimGrid1(NumberOfCandidates);		
	dim3 dimBlock1(NumberOfREACTIONSs);

	MetaINIT<<<dimGrid1, dimBlock1>>>(d_network, d_EFMs);	

	cudaProfilerStart(); 
	
	// number of blocks == number of metabolits :: chejuri ru block haye mokhtalef yek function ro active konam?
	dim3 dimGrid2(NumberOfCandidates);		//NumberOfMETABOLITEs = NumberOfBlocks
	dim3 dimBlock2(NumberOfMETABOLITEs);		//NumberOfThreads

	METAx<<<dimGrid2, dimBlock2>>>(d_network, d_EFMs, d_parVal/*d_EFMsss,*/ /*,d_hora*/);	//Numofblocks = NumberOfMETABOLITEs	BlockSize = numOfThreads = 10 (masalan)

	cudaProfilerStop();

	cudaPrintfDisplay (stdout, true);
	cudaPrintfEnd ();

	cudaDeviceSynchronize();
	//-----------------------------------------------------------------------------
	//-----------------------------------------------------------------------------
	//-----------------------------------------------------------------------------
	//-----------------------------------------------------------------------------  copy_back memory from gpu 
	cudaStatus = cudaMemcpy(&EFMs[0], &d_EFMs[0], sizeE/*(size_t)sizeof(*d_EFMs)*/, cudaMemcpyDeviceToHost);	
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	//----------------------------------------------------------------------------
Error:
	cudaFree(d_network);			
	cudaFree(d_EFMs);

	return cudaStatus;
}


//------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------

int main()
{

	//randPathStab();

	cudaError_t cudaStatus;

	char timeStr [9];

	string inFile1, outFile1;
	ifstream *fin;
	ofstream *fout;

	//-------------------- 

	vector<vector<float>> S;
	//vector<NODE> network;
	//NODE_ARRAY network[NumberOfMETABOLITEs];


	_strtime_s(timeStr);
	printf( "\nThe current time is %s \n\n", timeStr);

	// mona
	inFile1 = "inputs\\cho.txt";		//cho	//revTest_1000i				//iAF1260
	outFile1 = "outputs\\cho.txt";				//revTest_1000i

	fin = new ifstream(inFile1.c_str());
	fout = new ofstream(outFile1.c_str());

	//-------------------------------------------------------- create network
	S = readFileS(*fin);

	//network = CreateNetwork(S/*, network*/);
	/*network =*/ CreateNetwork_array(S/*, network*/);

	//--------------------------------------------------------
	cudaStatus = EFM();

	//--------------------------------------------------------
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "EFMFinder failed!");
		return 1;
	}

	//-------------------------------------------------------- results
	//printf();
	//eliminateSameEFMs();
	printEFMs(*fout); 

	// cudaDeviceReset must be called before exiting in order for profiling and
	// tracing tools such as Nsight and Visual Profiler to show complete traces.
	cudaStatus = cudaDeviceReset();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceReset failed!");
		return 1;
	}

	_strtime_s(timeStr);
	printf( "\nThe current time is %s \n\n", timeStr);

	//delete[] network; 

	printf("Press any key to continue...");
	_getch();

	return 0;
}

//----------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------
// desiging pseudo-recursive functions killing meeeee .. true story :| 
// blew my mind; 
// ---------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------

