#include "structs.h"
#include <cstdlib>


//***********************************************************************************************************
 vector <vector <float>> readFileS (istream & fin)
{
	vector <vector <float>> S;
	vector<float> temp;
	float n = 0;
	char c;

	fin>>NumbersOfMetabolits;
	fin>>NumbersOfReactions;

	for (int l=0 ; l<NumbersOfReactions; l++){
		fin >> c;
		IfBeingReversible.push_back(c);
	}

	for (int k=0 ; k<NumbersOfMetabolits; k++){
		temp.clear();
		for (int h=0 ; h<NumbersOfReactions; h++){
			fin >> n;
			temp.push_back(n);
		}
		S.push_back(temp);
	}
	return S;
}
//*********************************************************************************************************** 
vector<vector<int>> FindSource(vector<vector<float>>& S, int j, float temp)
{
	vector<vector<int>> back; 
	vector<int> source;
	vector<int> sibling;

	if (temp > 0)
	{
		for (int h=0; h < NumbersOfMetabolits; h++)
		{
			if (S.at(h).at(j) < 0)
			{
				source.push_back(h);
			}//end if
			else if ((S.at(h).at(j) > 0))
			{
				sibling.push_back(h);
			}

		}

	}// end if temp

	//-------------
	else if (temp < 0)
	{
		for (int h=0; h < NumbersOfMetabolits; h++)
		{
			if (S.at(h).at(j) > 0)
			{
				source.push_back(h);
			}//end if
			else if ((S.at(h).at(j) < 0))
			{
				sibling.push_back(h);
			}

		}

	}// end else temp

	back.push_back(source);
	back.push_back(sibling);

	return back;
}
//***********************************************************************************************************
 vector<NODE> CreateNetwork(vector<vector<float>> S/*, vector<NODE>& network*/){

	vector<NODE> network;
	vector<vector<int>> source;

	for (int i=0; i < NumbersOfMetabolits; i++)
	{
		NODE node;
		node.NodeName = i;
		for (int j=0; j < NumbersOfReactions; j++)
		{
			vector<float> temp = S.at(i);
			source = FindSource(S, j, temp.at(j));
			// reversibility check 
			if (IfBeingReversible.at(j) == '0') //if (IfBeingReversible.at(j) == 'i')			<<<<<<< irreversible == 0 >>>>>
			{
				if (temp.at(j)>0)
				{
					//source = FindSource(S, j, temp.at(j));
					INPUT input;
					input.setParams(source.at(0), source.at(1), j, temp.at(j), false, 0/*, false*/);
					node.inputs.push_back(input);
				}
				else if (temp.at(j)<0)
				{
					//source = FindSource(S, j, temp.at(j));
					OUTPUT output;
					output.setParams(source.at(1), source.at(0), j, temp.at(j), false, 0/*, false*/);
					node.outputs.push_back(output);
				}
			} //end if not reversible
			else  // if (IfBeingReversible.at(j) == 'r')										 <<<<< reversible == 1 >>>>>
			{
				if (temp.at(j)>0)
				{
					//source = FindSource(S, j, temp.at(j));
					INPUT input;
					input.setParams(source.at(0), source.at(1), j, temp.at(j), false, 0/*, false*/);
					node.inputs.push_back(input);

					OUTPUT output;
					output.setParams(source.at(1), source.at(0), j, /*-*/temp.at(j), true, 0/*, false*/);		// if reaction reversible ::: if coefficient > (like input) ::::: but if it's output
					node.outputs.push_back(output);
				}
				else if (temp.at(j)<0)
				{
					//source = FindSource(S, j, temp.at(j));
					INPUT input;
					input.setParams(source.at(0), source.at(1), j, /*-*/temp.at(j), true, 0/*, false*/);		// if reaction reversible ::: if coefficient < (like output) ::::: but if it's input 
					node.inputs.push_back(input);
					OUTPUT output;
					output.setParams(source.at(1), source.at(0), j, temp.at(j), false, 0/*, false*/);
					node.outputs.push_back(output);
				}

			}// end if reversible

		}// for j

		//node.ConsumeProduct = 0;

		network.push_back(node);
	}// for i

	return network;
}
//***********************************************************************************************************
/*NODE_ARRAY*/ /*void convertNODEtoNODEArray(NODE node, int num)
 {
	 NODE_ARRAY NA;

	 //std::vector<double> v;
	 //double* a = &v[0];

	 NA.METstatus = 0;

	 vector<INPUT_ARRAY> temp3;

	 for (int i=0; i<(int)node.inputs.size(); i++)
	 {
		INPUT_ARRAY temp; 
		//(&NA.inputs[i])->Coef = node.inputs.at(i).Coef;   
		//(&NA.inputs[i])->reacNum = node.inputs.at(i).reacNum; 
	
		temp.Coef = node.inputs.at(i).Coef;   
		temp.reacNum = node.inputs.at(i).reacNum; 
		temp.RECstatus = 0; 

		//for (int j=0; j<(int)node.inputs.at(i).metabolitNamesIn.size(); j++)
		//{
		vector<int> temp2; 
			//((&NA.inputs[i])->metabolitNamesIn)[j] = node.inputs.at(i).metabolitNamesIn.at(j);
		if (!node.inputs.at(i).metabolitNamesIn.empty()){ 
			//temp2 = node.inputs.at(i).metabolitNamesIn;
			temp.metabolitNamesIn = &(node.inputs.at(i).metabolitNamesIn)[0];
		}

		//}

		if (!node.inputs.at(i).metabolitNamesOut.empty()){ 
			//temp2 = node.inputs.at(i).metabolitNamesOut;
			temp.metabolitNamesOut = &(node.inputs.at(i).metabolitNamesOut)[0];
		}


		temp3.push_back(temp);
		//for (int j=0; j<(int)node.inputs.at(i).metabolitNamesOut.size(); j++)
		//{
		//	((&NA.inputs[i])->metabolitNamesOut)[j] = node.inputs.at(i).metabolitNamesOut.at(j); 
		//}
	 }

	 
	 if (!temp3.empty())
	 {
		NA.inputs = &temp3[0];
	 }

	 //---------------------------------------------
	 vector<OUTPUT_ARRAY> temp4;
	 for (int i=0; i<(int)node.outputs.size(); i++)
	 {
		OUTPUT_ARRAY temp; 

		temp.Coef = node.outputs.at(i).Coef;   
		temp.reacNum = node.outputs.at(i).reacNum; 
		temp.RECstatus = 0; 

		vector<int> temp2; 

		if (!node.outputs.at(i).metabolitNamesIn.empty()){ 
			//temp2 = node.outputs.at(i).metabolitNamesIn;
			temp.metabolitNamesIn = &(node.outputs.at(i).metabolitNamesIn)[0];
		}

		if (!node.outputs.at(i).metabolitNamesOut.empty()){ 
			//temp2 = node.outputs.at(i).metabolitNamesOut;
			temp.metabolitNamesOut = &(node.outputs.at(i).metabolitNamesOut)[0]; 
		}


		temp4.push_back(temp);

	 }

	 if (!temp4.empty())
	 {
		 NA.outputs = &temp4[0];
	 }
	// int alki = NA.inputs->metabolitNamesIn[0];

	 network[num] = &NA; 

	 return;
 }*/
//***********************************************************************************************************
 /*NODE_ARRAY * */ void CreateNetwork_array(vector<vector<float>> S/*, vector<NODE>& network*/)
 {
	 //NODE_ARRAY network [NumberOfMETABOLITEs];
	 //vector<NODE_ARRAY> net;

	 vector<vector<int>> source;

	for (int i=0; i < NumbersOfMetabolits; i++)
	{
		NODE node;
		node.NodeName = i;
		for (int j=0; j < NumbersOfReactions; j++)
		{
			vector<float> temp = S.at(i);
			source = FindSource(S, j, temp.at(j));
			// reversibility check 
			if (IfBeingReversible.at(j) == '0') //if (IfBeingReversible.at(j) == 'i')			<<<<<<< irreversible == 0 >>>>>
			{
				if (temp.at(j)>0)
				{
					//source = FindSource(S, j, temp.at(j));
					INPUT input;
					input.setParams(source.at(0), source.at(1), j, temp.at(j), false, 0/*, false*/);
					node.inputs.push_back(input);
				}
				else if (temp.at(j)<0)
				{
					//source = FindSource(S, j, temp.at(j));
					OUTPUT output;
					output.setParams(source.at(1), source.at(0), j, temp.at(j), false, 0/*, false*/);
					node.outputs.push_back(output);
				}
			} //end if not reversible
			else  // if (IfBeingReversible.at(j) == 'r')										 <<<<< reversible == 1 >>>>>
			{
				if (temp.at(j)>0)
				{
					//source = FindSource(S, j, temp.at(j));
					INPUT input;
					input.setParams(source.at(0), source.at(1), j, temp.at(j), false, 0/*, false*/);
					node.inputs.push_back(input);

					OUTPUT output;
					output.setParams(source.at(1), source.at(0), j, /*-*/temp.at(j), true, 0/*, false*/);		// if reaction reversible ::: if coefficient > (like input) ::::: but if it's output
					node.outputs.push_back(output);
				}
				else if (temp.at(j)<0)
				{
					//source = FindSource(S, j, temp.at(j));
					INPUT input;
					input.setParams(source.at(0), source.at(1), j, /*-*/temp.at(j), true, 0/*, false*/);		// if reaction reversible ::: if coefficient < (like output) ::::: but if it's input 
					node.inputs.push_back(input);
					OUTPUT output;
					output.setParams(source.at(1), source.at(0), j, temp.at(j), false, 0/*, false*/);
					node.outputs.push_back(output);
				}

			}// end if reversible

		}// for j

		//node.ConsumeProduct = 0;

		net.push_back(node);
		/*NODE_ARRAY NA =  *///convertNODEtoNODEArray(node, i);
		//net.push_back(NA);
		//network[i] =  NA;

		//---------------------------------------------------------------------------------------------
		//---------------------------------------------------------------------------------------------
		NODE_ARRAY NA;// = new NODE_ARRAY();
		//INPUT_ARRAY *iinnppss = new INPUT_ARRAY[3];

		 NA.METstatus = 0;

		 int inoutSize = (int)node.inputs.size();		//inoutSize &= x000F; 	  
		 NA.METstatus |= inoutSize;

		 inoutSize = (int)node.outputs.size() << 4;		//inoutSize &= x00F0;
		 NA.METstatus |= inoutSize;

		
		// sizeCount_glob += sizeof(int);				// for metStatus
		

		 //vector<INPUT_ARRAY> *temp3;

		 //NA->inputs = new INPUT_ARRAY[];	//(int)node.inputs.size()
		  //NA->inputs = new INPUT_ARRAY[(int)node.inputs.size()];

		 for (int j=0; j<(int)node.inputs.size(); j++)
		 {
			//INPUT_ARRAY temp; 
			//INPUT_ARRAY *temp = new INPUT_ARRAY();
	
			//NA->inputs[i].Coef = node.inputs.at(i).Coef;
			
			// NA->inputs = new INPUT_ARRAY();

			NA.inputs[j].Coef = node.inputs.at(j).Coef;   
			NA.inputs[j].reacNum = node.inputs.at(j).reacNum; 
			//NA.reacCover[NA.inputs[j].reacNum] = 1;
			NA.inputs[j].RECstatus = 0; 
			NA.inputs[j].notGoodPrimaryCandida = 0; 

			inoutSize = (int)node.inputs.at(j).metabolitNamesIn.size();
			NA.inputs[j].RECstatus |= inoutSize;

			
			inoutSize = (int)node.inputs.at(j).metabolitNamesOut.size() << 8;
			NA.inputs[j].RECstatus |= inoutSize;

			//sizeCount_glob += 3*sizeof(int);											// for recStatus, Coef, reacNum		

			//NA->inputs[j].metabolitNamesIn = new int[(int)node.inputs.at(j).metabolitNamesIn.size()]; //
			for (int k=0; k<(int)node.inputs.at(j).metabolitNamesIn.size(); k++){
			
			//vector<int> temp2; 
				//((&NA.inputs[i])->metabolitNamesIn)[j] = node.inputs.at(i).metabolitNamesIn.at(j);
			//if (!node.inputs.at(i).metabolitNamesIn.empty()){ 
				//temp2 = node.inputs.at(i).metabolitNamesIn;
				
				NA.inputs[j].metabolitNamesIn[k] = node.inputs.at(j).metabolitNamesIn[k];
			}

			//sizeCount_glob += node.inputs.at(j).metabolitNamesIn.size()*sizeof(int);
			

			//}
			//NA->inputs[j].metabolitNamesOut = new int[(int)node.inputs.at(j).metabolitNamesOut.size()]; //
			//if (!node.inputs.at(i).metabolitNamesOut.empty()){ 
			for (int k=0; k<(int)node.inputs.at(j).metabolitNamesOut.size(); k++){
				//temp2 = node.inputs.at(i).metabolitNamesOut;
				//NA->inputs[i].metabolitNamesOut = &(node.inputs.at(i).metabolitNamesOut)[0];
				
				NA.inputs[j].metabolitNamesOut[k] = node.inputs.at(j).metabolitNamesOut.at(k);
			}

		    //sizeCount_glob += node.inputs.at(j).metabolitNamesOut.size()*sizeof(int);	

			//(&temp3)push_back(temp);  
		    //(NA->inputs) = (temp); 
			//(NA->inputs)++;
		 }
		 
	 
/*		 if (!temp3->empty())
		 {
			NA->inputs = temp3[0]; 
		 }*/

		 //---------------------------------------------
		 //vector<OUTPUT_ARRAY> temp4;

		 //NA->outputs = new OUTPUT_ARRAY[(int)node.outputs.size()];//

		 for (int j=0; j<(int)node.outputs.size(); j++)
		 {
			// OUTPUT_ARRAY temp; 
			// OUTPUT_ARRAY *temp = new OUTPUT_ARRAY();
			// OUTPUT_ARRAY *temp = new OUTPUT_ARRAY();

			//NA->outputs = new OUTPUT_ARRAY();

			NA.outputs[j].Coef = node.outputs.at(j).Coef;   
			NA.outputs[j].reacNum = node.outputs.at(j).reacNum; 
			//NA.reacCover[NA.outputs[j].reacNum] = 1;
			NA.outputs[j].RECstatus = 0; 
			NA.outputs[j].notGoodPrimaryCandida = 0; 

			inoutSize = (int)node.outputs.at(j).metabolitNamesIn.size();
			NA.outputs[j].RECstatus |= inoutSize;

			
			inoutSize = (int)node.outputs.at(j).metabolitNamesOut.size() << 8;
			NA.outputs[j].RECstatus |= inoutSize;


			//sizeCount_glob += 3*sizeof(int);											// for recStatus, Coef, reacNum	


			//vector<int> temp2; 

			//NA->outputs[j].metabolitNamesIn = new int[(int)node.outputs.at(j).metabolitNamesIn.size()]; //
			//if (!node.outputs.at(i).metabolitNamesIn.empty()){ 
			for (int k=0; k<(int)node.outputs.at(j).metabolitNamesIn.size(); k++){
				//temp2 = node.outputs.at(i).metabolitNamesIn;
				//NA->outputs[i].metabolitNamesIn = &(node.outputs.at(i).metabolitNamesIn)[0];
				
				NA.outputs[j].metabolitNamesIn[k] = node.outputs.at(j).metabolitNamesIn.at(k);
			}
			//sizeCount_glob += node.outputs.at(j).metabolitNamesIn.size()*sizeof(int);

			//NA->outputs[j].metabolitNamesIn_size = (int)node.outputs.at(j).metabolitNamesOut.size();
			//NA->outputs[j].metabolitNamesOut = new int[(int)node.outputs.at(j).metabolitNamesOut.size()]; //NA->outputs[j].metabolitNamesIn_size
			for (int k=0; k<(int)node.outputs.at(j).metabolitNamesOut.size(); k++){  //NA->outputs[j].metabolitNamesIn_size
			//if (!node.outputs.at(i).metabolitNamesOut.empty()){ 
				//temp2 = node.outputs.at(i).metabolitNamesOut;
				//NA->outputs[i].metabolitNamesOut = &(node.outputs.at(i).metabolitNamesOut)[0]; 
				
				NA.outputs[j].metabolitNamesOut[k] = node.outputs.at(j).metabolitNamesOut.at(k);
			}

			//sizeCount_glob += node.outputs.at(j).metabolitNamesOut.size()*sizeof(int);
			  
			//temp4.push_back(*temp);
			//NA->outputs = temp; 
			//(NA->outputs)++;

		 }

		 /*if (!temp4.empty())
		 {
			 NA->outputs = &temp4[0];
		 }*/

		 network[0].net[i] = NA; 	

		//---------------------------------------------------------------------------------------------
		//---------------------------------------------------------------------------------------------

	}// for i

	//network = &net[0];
	//int alaki = network[0]->inputs[0].metabolitNamesIn[1];
		
	// return network;  
 }
//***********************************************************************************************************
 void printEFMs(ofstream &fout) // check also for EFMs with "the same subset of reactions" ==> Just for test
{

	fout << NumberOfCandidates <<"\n";			//check the similarity of EFMs 
	fout << NumbersOfReactions <<"\n";

	for (int i=0; i< (int)NumberOfCandidates; i++)
	{
		//-------------
		//if (EFMs[i].recFlux[NumberOfREACTIONSsPlus-1]==1){
			for (int j = 0; j < NumberOfREACTIONSsPlus; j++)
			{
				fout << EFMs[i].recFlux[j] << "\t";

				if (j==(NumberOfREACTIONSsPlus-2))
					fout << "\t";
			}
			fout << "\n";
		//}
		//------------- 
	}// end for number of EFMs
}
//***********************************************************************************************************
//***********************************************************************************************************
//***********************************************************************************************************
//***********************************************************************************************************
//***********************************************************************************************************
//***********************************************************************************************************
 int setMet (int v, vector<int> F)
 {
	 int o = 0; 
	//if (v==0)
		 //

	if(v==1)
		 o = F.at(6) - F.at(7); 

	else if(v==2)
		 o = F.at(15) - F.at(6); 

	else if(v==3)
		 o = F.at(14) - F.at(15);

	else if(v==4)
		 o = F.at(5) - F.at(14); 

	else if(v==5)
		 o = F.at(26) - F.at(5); 

	else if(v==6)
		o = F.at(11) - (F.at(5)+ F.at(26)); 

	else if(v==7)
		 o =  F.at(8) - F.at(11); 

	else if(v==8)
		 o =  F.at(9) - F.at(8); 

	else if(v==9)
		 o =  (F.at(4) + F.at(17)) - (F.at(6) + F.at(9)); 

	else if(v==10)
		o = (F.at(6) + F.at(16)) - F.at(17); 

	else if(v==11)
		o =   F.at(22) - F.at(16); 

	else if(v==12)
		o = F.at(23) - (F.at(13)+F.at(25)); 

	else if(v==13)
		o =   F.at(13) - F.at(23); 

	else if(v==14)
		o =   F.at(12) - F.at(13); 

	//else if(v==15)
		 //

	else if(v==16)
		o = (F.at(22) + F.at(25)) - F.at(23); 

	else if(v==17)
		o =   F.at(25) - F.at(20); 

	else if(v==18)
		o =   F.at(20) - F.at(3); 

	else if(v==19)
		o =   F.at(20) - F.at(19);

	else if(v==20)
		o =   F.at(18) - F.at(25); 

	else if(v==21)
		o =   F.at(2) - F.at(18); 

	else if(v==22)
		o =   F.at(19) - F.at(18); 

	else if(v==23)
		o =   F.at(1) - F.at(2); 

	else if(v==24)
		o =   F.at(10) - F.at(1); 

	else if(v==25)
		o =   F.at(24) - F.at(10); 

	//else if(v==26)
		 //

	else if(v==27)
		o =   F.at(3) - F.at(4); 

	else if(v==28)
		o =   F.at(4) - F.at(21); 

	else if(v==29)
		o =   F.at(21) - F.at(0); 

	else if(v==30)
		o =   F.at(0) - F.at(22);

	return o;
 }
 //----
 int ifAny(vector<int> m)
 {
	 int o = 0;
	 for (int i=0; i<(int)m.size(); i++)
	 {
		if (m.at(i)!=0)
		{
			o = 1;
			return o;
		}
	 }
	 return o;
 }
 //----
 vector<int> adopt(int flag, vector<int> Fluxes, int j, int jVal)
 {
	 int jValOrg = jVal;
	 jVal = abs(jVal);
	if (flag==0)
	{

		if (jValOrg > 0) //backward
		{
			if (j==1)
			{ 	Fluxes.at(7) += jVal;		}
			else if(j==2)
			{ 	Fluxes.at(6) += jVal;		}
			else if(j==3)
			{ 	Fluxes.at(15) += jVal;		}
			else if(j==4)
			{ 	Fluxes.at(14) += jVal;		}
			else if(j==5)
			{ 	Fluxes.at(5) += jVal;		}
			else if(j==6)
			{ 	Fluxes.at(26) += jVal;		}
			else if(j==7)
			{ 	Fluxes.at(11) += jVal;		}
			else if(j==8)
			{ 	Fluxes.at(8) += jVal;		}
			else if(j==9)
			{ 	Fluxes.at(9) += jVal;		}
			else if(j==10)
			{ 	Fluxes.at(17) += jVal;		}
			else if(j==11)
			{ 	Fluxes.at(16) += jVal;		}
			else if(j==12)
			{ 	Fluxes.at(22) += jVal;		}
			else if(j==13)
			{ 	Fluxes.at(23) += jVal;		}
			else if(j==14)
			{ 	Fluxes.at(13) += jVal;		}
			else if(j==16)
			{ 	Fluxes.at(23) += jVal;		}
			else if(j==17)
			{ 	Fluxes.at(20) += jVal;		}
			else if(j==18)
			{ 	Fluxes.at(3) += jVal;		}
			else if(j==19)
			{ 	Fluxes.at(19) += jVal;		}
			else if(j==20)
			{ 	Fluxes.at(25) += jVal;		}
			else if(j==21)
			{ 	Fluxes.at(18) += jVal;		}
			else if(j==22)
			{ 	Fluxes.at(18) += jVal;		}
			else if(j==23)
			{ 	Fluxes.at(2) += jVal;		}
			else if(j==24)
			{ 	Fluxes.at(1) += jVal;		}
			else if(j==25)
			{ 	Fluxes.at(10) += jVal;		}
			else if(j==27)
			{ 	Fluxes.at(4) += jVal;		}
			else if(j==28)
			{ 	Fluxes.at(21) += jVal;		}
			else if(j==29)
			{ 	Fluxes.at(0) += jVal;		}
			else if(j==30)
			{ 	Fluxes.at(22) += jVal;		}

		}
		else if (jValOrg < 0) //forward 
		{
			if (j==1)
			{ 	Fluxes.at(6) += jVal;		}
			else if(j==2)
			{ 	Fluxes.at(15) += jVal;		}
			else if(j==3)
			{ 	Fluxes.at(14) += jVal;		}
			else if(j==4)
			{ 	Fluxes.at(5) += jVal;		}
			else if(j==5)
			{ 	Fluxes.at(26) += jVal;		}
			else if(j==6)
			{ 	Fluxes.at(11) += jVal;		}
			else if(j==7)
			{ 	Fluxes.at(8) += jVal;		}
			else if(j==8)
			{ 	Fluxes.at(9) += jVal;		}
			else if(j==9)
			{ 	Fluxes.at(4) += jVal;		}
			else if(j==10)
			{ 	Fluxes.at(16) += jVal;		}
			else if(j==11)
			{ 	Fluxes.at(22) += jVal;		}
			else if(j==12)
			{ 	Fluxes.at(23) += jVal;		}
			else if(j==13)
			{ 	Fluxes.at(13) += jVal;		}
			else if(j==14)
			{ 	Fluxes.at(12) += jVal;		}
			else if(j==16)
			{ 	Fluxes.at(25) += jVal;		}
			else if(j==17)
			{ 	Fluxes.at(25) += jVal;		}
			else if(j==18)
			{ 	Fluxes.at(20) += jVal;		}
			else if(j==19)
			{ 	Fluxes.at(20) += jVal;		}
			else if(j==20)
			{ 	Fluxes.at(18) += jVal;		}
			else if(j==21)
			{ 	Fluxes.at(2) += jVal;		}
			else if(j==22)
			{ 	Fluxes.at(19) += jVal;		}
			else if(j==23)
			{ 	Fluxes.at(1) += jVal;		}
			else if(j==24)
			{ 	Fluxes.at(10) += jVal;		}
			else if(j==25)
			{ 	Fluxes.at(24) += jVal;		}
			else if(j==27)
			{ 	Fluxes.at(3) += jVal;		}
			else if(j==28)
			{ 	Fluxes.at(4) += jVal;		}
			else if(j==29)
			{ 	Fluxes.at(21) += jVal;		}
			else if(j==30)
			{ 	Fluxes.at(0) += jVal;		}

		}
	}//end if flag

	//.....................Visited ::: select a random input/output 
	else 
	{
		int rd = rand()%2; 
		if ((rd==0)&& (jValOrg>0)){
			if (j==1)
			{ 	Fluxes.at(7) += jVal;		}
			else if(j==2)
			{ 	Fluxes.at(6) += jVal;		}
			else if(j==3)
			{ 	Fluxes.at(15) += jVal;		}
			else if(j==4)
			{ 	Fluxes.at(14) += jVal;		}
			else if(j==5)
			{ 	Fluxes.at(5) += jVal;		}
			else if(j==6)
			{ 	Fluxes.at(26) += jVal;		}
			else if(j==7)
			{ 	Fluxes.at(11) += jVal;		}
			else if(j==8)
			{ 	Fluxes.at(8) += jVal;		}
			else if(j==9)
			{ 	Fluxes.at(9) += jVal;		}
			else if(j==10)
			{ 	Fluxes.at(17) += jVal;		}
			else if(j==11)
			{ 	Fluxes.at(16) += jVal;		}
			else if(j==12)
			{ 	Fluxes.at(22) += jVal;		}
			else if(j==13)
			{ 	Fluxes.at(23) += jVal;		}
			else if(j==14)
			{ 	Fluxes.at(13) += jVal;		}
			else if(j==16)
			{ 	Fluxes.at(23) += jVal;		}
			else if(j==17)
			{ 	Fluxes.at(20) += jVal;		}
			else if(j==18)
			{ 	Fluxes.at(3) += jVal;		}
			else if(j==19)
			{ 	Fluxes.at(19) += jVal;		}
			else if(j==20)
			{ 	Fluxes.at(25) += jVal;		}
			else if(j==21)
			{ 	Fluxes.at(18) += jVal;		}
			else if(j==22)
			{ 	Fluxes.at(18) += jVal;		}
			else if(j==23)
			{ 	Fluxes.at(2) += jVal;		}
			else if(j==24)
			{ 	Fluxes.at(1) += jVal;		}
			else if(j==25)
			{ 	Fluxes.at(10) += jVal;		}
			else if(j==27)
			{ 	Fluxes.at(4) += jVal;		}
			else if(j==28)
			{ 	Fluxes.at(21) += jVal;		}
			else if(j==29)
			{ 	Fluxes.at(0) += jVal;		}
			else if(j==30)
			{ 	Fluxes.at(22) += jVal;		}
		}
		else if ((rd==1) && (jValOrg<0))
		{
			if (j==1)
			{ 	Fluxes.at(6) += jVal;		}
			else if(j==2)
			{ 	Fluxes.at(15) += jVal;		}
			else if(j==3)
			{ 	Fluxes.at(14) += jVal;		}
			else if(j==4)
			{ 	Fluxes.at(5) += jVal;		}
			else if(j==5)
			{ 	Fluxes.at(26) += jVal;		}
			else if(j==6)
			{ 	Fluxes.at(11) += jVal;		}
			else if(j==7)
			{ 	Fluxes.at(8) += jVal;		}
			else if(j==8)
			{ 	Fluxes.at(9) += jVal;		}
			else if(j==9)
			{ 	Fluxes.at(4) += jVal;		}
			else if(j==10)
			{ 	Fluxes.at(16) += jVal;		}
			else if(j==11)
			{ 	Fluxes.at(22) += jVal;		}
			else if(j==12)
			{ 	Fluxes.at(23) += jVal;		}
			else if(j==13)
			{ 	Fluxes.at(13) += jVal;		}
			else if(j==14)
			{ 	Fluxes.at(12) += jVal;		}
			else if(j==16)
			{ 	Fluxes.at(25) += jVal;		}
			else if(j==17)
			{ 	Fluxes.at(25) += jVal;		}
			else if(j==18)
			{ 	Fluxes.at(20) += jVal;		}
			else if(j==19)
			{ 	Fluxes.at(20) += jVal;		}
			else if(j==20)
			{ 	Fluxes.at(18) += jVal;		}
			else if(j==21)
			{ 	Fluxes.at(2) += jVal;		}
			else if(j==22)
			{ 	Fluxes.at(19) += jVal;		}
			else if(j==23)
			{ 	Fluxes.at(1) += jVal;		}
			else if(j==24)
			{ 	Fluxes.at(10) += jVal;		}
			else if(j==25)
			{ 	Fluxes.at(24) += jVal;		}
			else if(j==27)
			{ 	Fluxes.at(3) += jVal;		}
			else if(j==28)
			{ 	Fluxes.at(4) += jVal;		}
			else if(j==29)
			{ 	Fluxes.at(21) += jVal;		}
			else if(j==30)
			{ 	Fluxes.at(0) += jVal;		}
		}
		else if ((rd==0)&& (jValOrg<0)){
			if (j==1)
			{ 	Fluxes.at(7) -= jVal;		}
			else if(j==2)
			{ 	Fluxes.at(6) -= jVal;		}
			else if(j==3)
			{ 	Fluxes.at(15) -= jVal;		}
			else if(j==4)
			{ 	Fluxes.at(14) -= jVal;		}
			else if(j==5)
			{ 	Fluxes.at(5) -= jVal;		}
			else if(j==6)
			{ 	Fluxes.at(26) -= jVal;		}
			else if(j==7)
			{ 	Fluxes.at(11) -= jVal;		}
			else if(j==8)
			{ 	Fluxes.at(8) -= jVal;		}
			else if(j==9)
			{ 	Fluxes.at(9) -= jVal;		}
			else if(j==10)
			{ 	Fluxes.at(17) -= jVal;		}
			else if(j==11)
			{ 	Fluxes.at(16) -= jVal;		}
			else if(j==12)
			{ 	Fluxes.at(22) -= jVal;		}
			else if(j==13)
			{ 	Fluxes.at(23) -= jVal;		}
			else if(j==14)
			{ 	Fluxes.at(13) -= jVal;		}
			else if(j==16)
			{ 	Fluxes.at(23) -= jVal;		}
			else if(j==17)
			{ 	Fluxes.at(20) -= jVal;		}
			else if(j==18)
			{ 	Fluxes.at(3) -= jVal;		}
			else if(j==19)
			{ 	Fluxes.at(19) -= jVal;		}
			else if(j==20)
			{ 	Fluxes.at(25) -= jVal;		}
			else if(j==21)
			{ 	Fluxes.at(18) -= jVal;		}
			else if(j==22)
			{ 	Fluxes.at(18) -= jVal;		}
			else if(j==23)
			{ 	Fluxes.at(2) -= jVal;		}
			else if(j==24)
			{ 	Fluxes.at(1) -= jVal;		}
			else if(j==25)
			{ 	Fluxes.at(10) -= jVal;		}
			else if(j==27)
			{ 	Fluxes.at(4) -= jVal;		}
			else if(j==28)
			{ 	Fluxes.at(21) -= jVal;		}
			else if(j==29)
			{ 	Fluxes.at(0) -= jVal;		}
			else if(j==30)
			{ 	Fluxes.at(22) -= jVal;		}
		}
		else if ((rd==1) && (jValOrg>0))
		{
			if (j==1)
			{ 	Fluxes.at(6) -= jVal;		}
			else if(j==2)
			{ 	Fluxes.at(15) -= jVal;		}
			else if(j==3)
			{ 	Fluxes.at(14) -= jVal;		}
			else if(j==4)
			{ 	Fluxes.at(5) -= jVal;		}
			else if(j==5)
			{ 	Fluxes.at(26) -= jVal;		}
			else if(j==6)
			{ 	Fluxes.at(11) -= jVal;		}
			else if(j==7)
			{ 	Fluxes.at(8) -= jVal;		}
			else if(j==8)
			{ 	Fluxes.at(9) -= jVal;		}
			else if(j==9)
			{ 	Fluxes.at(4) -= jVal;		}
			else if(j==10)
			{ 	Fluxes.at(16) -= jVal;		}
			else if(j==11)
			{ 	Fluxes.at(22) -= jVal;		}
			else if(j==12)
			{ 	Fluxes.at(23) -= jVal;		}
			else if(j==13)
			{ 	Fluxes.at(13) -= jVal;		}
			else if(j==14)
			{ 	Fluxes.at(12) -= jVal;		}
			else if(j==16)
			{ 	Fluxes.at(25) -= jVal;		}
			else if(j==17)
			{ 	Fluxes.at(25) -= jVal;		}
			else if(j==18)
			{ 	Fluxes.at(20) -= jVal;		}
			else if(j==19)
			{ 	Fluxes.at(20) -= jVal;		}
			else if(j==20)
			{ 	Fluxes.at(18) -= jVal;		}
			else if(j==21)
			{ 	Fluxes.at(2) -= jVal;		}
			else if(j==22)
			{ 	Fluxes.at(19) -= jVal;		}
			else if(j==23)
			{ 	Fluxes.at(1) -= jVal;		}
			else if(j==24)
			{ 	Fluxes.at(10) -= jVal;		}
			else if(j==25)
			{ 	Fluxes.at(24) -= jVal;		}
			else if(j==27)
			{ 	Fluxes.at(3) -= jVal;		}
			else if(j==28)
			{ 	Fluxes.at(4) -= jVal;		}
			else if(j==29)
			{ 	Fluxes.at(21) -= jVal;		}
			else if(j==30)
			{ 	Fluxes.at(0) -= jVal;		}
		}
	}//end else flag 

	return Fluxes;
 }

 //***********************************************************************************************************
void randPathStab (void)
{
	vector<int>		Mets;
	vector<int>		Fluxes;
	vector<bool>	vis;

	for (int i=0; i<27; i++)
	{
		Fluxes.push_back(0);
	}

	for (int i=0; i<31; i++)
	{
		Mets.push_back(0);
	}

	for (int i=0; i<31; i++)
	{
		vis.push_back(0);
	}

	//set up 
	Fluxes.at(7) = 1;

	for (int k=0; k<31; k++)
	{
		Mets.at(k) = setMet(k, Fluxes);
	}

	//while the sum of Mets are ~= 0 
	int count =0; 
	while (ifAny(Mets) != 0) 
	{
		cout << count++ << "\t"; 
		for (int j=0; j<31; j++)
		{
			if (Mets.at(j)!=0)
			{
				int flag;
				if (vis.at(j) == 0)
				{vis.at(j)=1;	flag=0;}	//nonVisited
				else	{flag = 1;}			//visited

				//change input or output fluxes with adopt 
				Fluxes = adopt(flag, Fluxes, j, Mets.at(j));		//update fluxes
				for (int k=0; k<31; k++)
				{ 	Mets.at(k) = setMet(k, Fluxes); }				//update Mets accordingly 


			}//end if 
		}//end for

	}//end while

	for (int h=0; h<Fluxes.size();h++)
	{
		cout << Fluxes.at(h) << "\t";
	}
	cout << "\n\n";

}
 //***********************************************************************************************************