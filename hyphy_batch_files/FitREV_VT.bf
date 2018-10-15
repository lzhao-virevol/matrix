SetDialogPrompt ("Alignment file list:");
fscanf 			(PROMPT_FOR_FILE,"Lines",inFiles);
filesRead		= Columns (inFiles);
fprintf			(stdout, "Read ", filesRead, " files\n");

jointFreqs  	= {20,1};
jointChar   	= 0;
treeStrings 	= {filesRead,1}; 
lfDef			= "";

fprintf (stdout, "[PHASE 1.] Reading and checking sequence alignments...\n");

COUNT_GAPS_IN_FREQUENCIES = 0;
for (k=0; k<filesRead; k=k+1)
{
	ExecuteCommands ("DataSet ds_" + k + " = ReadDataFile (\"" + inFiles[k] + "\");\n");
	if (!IS_TREE_PRESENT_IN_DATA)
	{
		fprintf (stdout, "[***ERROR***] No tree read with alignment from file ", inFiles[k], ". Terminating...\n");
		return 0;
	}
	ExecuteCommands ("spc = ds_" + k + ".species;stc = ds_" + k + ".sites;usc=ds_"+k+".unique_sites;");
	ExecuteCommands ("DataSetFilter dsf_" + k + " = CreateFilter (ds_"+k+",1);\n");
	ExecuteCommands ("treeStrings["+k+"] = DATAFILE_TREE\n;");
	ExecuteCommands ("HarvestFrequencies (freqs, ds_" + k + ",1,1,1);");
	fprintf (stdout, "\tRead file ", inFiles[k], " with ", spc, " sequences, ", stc, " residues and ",usc," unique alignment columns.\n");
	jointChar = jointChar + spc*stc;
	jointFreqs = jointFreqs + freqs*spc*stc;
	if (k)
	{
		lfDef = lfDef + ",";
	}
	lfDef = lfDef + "dsf_" + k + ",T_" + k;
}
COUNT_GAPS_IN_FREQUENCIES = 1;
vectorOfFrequencies = jointFreqs * (1/jointChar);

ExecuteAFile 		(HYPHY_BASE_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "TemplateModels" + DIRECTORY_SEPARATOR + "defineGamma.mdl");
SetDialogPrompt 	("Save the resulting matrix to:");
fprintf 			(PROMPT_FOR_FILE, CLEAR_FILE);
matrixOutFile	= LAST_FILE_PATH;
lfFitFile		= matrixOutFile + ".likelihood";

/*-------------------------------------------------------------*/

allowedAACharacters = "ACDEFGHIKLMNPQRSTVWY";
fprintf (stdout, "\n[EMPIRICAL BASE COMPOSITION]\n");
for (k=0; k<20; k=k+1)
{
	fprintf (stdout, "\t", allowedAACharacters[k], ":", Format (vectorOfFrequencies[k]*100,8,3), "%\n");
}

PopulateModelMatrixWAG ("WAGMatrixF", vectorOfFrequencies);
Model WAG_MF = (WAGMatrixF, vectorOfFrequencies, -1);

fprintf (stdout, "[PHASE 2.] Fitting the WAG model to approximate branch lengths...\n");

overallCMX = {20,20};
approxL    = 0;

for (fileID =0; fileID<filesRead; fileID=fileID+1)
{
	ExecuteCommands 	("Tree WAG_Tree_"+fileID+ "= treeStrings[fileID];");
	ExecuteCommands 	("LikelihoodFunction WAG_LF = (dsf_"+fileID+",WAG_Tree_"+fileID+");");
	Optimize 			(res, WAG_LF);
	ExecuteCommands		("tbl = BranchLength (WAG_Tree_"+fileID+",-1);");
	
	tbl = (tbl*Transpose (tbl["1"]))[0];
	fprintf (stdout, "\tFinished with dataset ", fileID+1, ". Log(L) = ", Format (res[1][0],10,3), " Estimateed tree length (subs/site): ", 
					   tbl, "\n\t\tTabulating ML-based substitution counts...\n");
					   
	approxL = approxL + res[1][0];
	localCMX = countSubstitutionTypes (fileID);
	overallCMX = overallCMX + localCMX;
	prof	   = reportSubstitutionMatrix (localCMX,0);
	fprintf (stdout, "\t\t", Format (prof[0], 3,0), "/190 possible substitution kinds, with at least ", Format (prof[1], 8,0), " substitutions\n");
}


fprintf (stdout, "[APPROXIMATE JOINT LOG(L)]\n\t", Format (approxL, 8, 3), "\n");
fprintf (stdout, "[OVERALL SUBSTITUTION PROFILE]\n");
prof = reportSubstitutionMatrix (overallCMX,1);
fprintf (stdout, "\n\t", Format (prof[0], 3,0), "/190 possible substitution kinds, with at least ", Format (prof[1], 8,0), " substitutions\n");
fprintf (stdout, "\n[PHASE 3.] Fitting the REV model - this is going to take a while, most likely.\n");

/*-------------------------------------------------------------*/

aaNames = {{"Alanine",
"Cysteine",
"Aspartic_Acid",
"Glutamic_Acid",
"Phenylalanine",
"Glycine",
"Histidine",
"Isoleucine",
"Lysine",
"Leucine",
"Methionine",
"Asparagine",
"Proline",
"Glutamine",
"Arginine",
"Serine",
"Threonine",
"Valine",
"Tryptophan",
"Tyrosine"}};

t   = 1;
dMx = WAGMatrixF;

PopulateModelMatrix ("REVQ",vectorOfFrequencies);
Model AA_REV 	  = (REVQ, vectorOfFrequencies, 1);

for (k=0; k<filesRead; k=k+1)
{
	ExecuteCommands ("Tree T_" + k + " = " + treeStrings[k] + ";ReplicateConstraint(\"this1.?.t:=this2.?.t\", T_"+k+",WAG_Tree_"+k+"); ClearConstraints (T_"+k+");\n");
}

ExecuteCommands ("LikelihoodFunction lf = (" + lfDef + ");\n");
OPTIMIZATION_METHOD     = 0;
USE_LAST_RESULTS	= 1;
VERBOSITY_LEVEL         = 1;
Optimize (res,lf);
USE_LAST_RESULTS	= 0;
OPTIMIZATION_METHOD     = 4;

fprintf (stdout, lf);
saveLFO = LIKELIHOOD_FUNCTION_OUTPUT;
LIKELIHOOD_FUNCTION_OUTPUT = 6;
fprintf (lfFitFile, CLEAR_FILE, lf);
LIKELIHOOD_FUNCTION_OUTPUT = saveLFO;


fprintf (stdout, "\n[PHASE 4.] Imputing zero rates\n");

scalingFactor = 0;
for (fileID =0; fileID<filesRead; fileID=fileID+1)
{
	ExecuteCommands		("tbl = BranchLength (T_"+fileID+",-1);");
	
	tbl = (tbl*Transpose (tbl["1"]))[0];
	ExecuteCommands 	("scalingFactor=scalingFactor+tbl*(dsf_"+fileID+".sites+1);");
}


t 			 = 1;
REVQ_Numeric = REVQ;

mxScaler	 = 0;
for (k=0; k<19; k=k+1)
{
	for (k2 = k+1; k2 < 20; k2=k2+1)
	{
		mxScaler = mxScaler + (REVQ_Numeric[k][k2]+REVQ_Numeric[k2][k])*vectorOfFrequencies[k]*vectorOfFrequencies[k2];
	}
}    

REVQ_Numeric = REVQ_Numeric * (1/mxScaler);

for (k=0; k<19; k=k+1)
{
	for (k2 = k+1; k2 < 20; k2=k2+1)
	{
		if (REVQ[k][k2] == 0)
		{
			thisS = 1/((vectorOfFrequencies[k]+vectorOfFrequencies[k2])*scalingFactor);
			REVQ_Numeric[k][k2] = thisS;
			REVQ_Numeric[k2][k] = thisS;
			fprintf (stdout, "\t[", allowedAACharacters[k], "<->", allowedAACharacters[k2], "] ", Format (thisS, 8, 0),"\n");
		}	
	}
}    

mxScaler2	 = 0;
for (k=0; k<19; k=k+1)
{
	for (k2 = k+1; k2 < 20; k2=k2+1)
	{
		mxScaler2 = mxScaler2 + (REVQ_Numeric[k][k2]+REVQ_Numeric[k2][k])*vectorOfFrequencies[k]*vectorOfFrequencies[k2];
	}
}    

REVQ_Numeric = REVQ_Numeric * (1/mxScaler2);
fprintf (matrixOutFile, CLEAR_FILE, allowedAACharacters, "\n", REVQ_Numeric, "\n", vectorOfFrequencies);

/*-------------------------------------------------------------*/

function PopulateModelMatrix (ModelMatrixName&, EFV)
{

	ModelMatrixName = {20,20};
	if (categoriesUsed)
	{
		for (k=0; k<19; k=k+1)
		{
			for (k2 = k+1; k2 < 20; k2=k2+1)
			{
				subString = aaNames[k]+aaNames[k2];
				ExecuteCommands("global "+subString +"=" + dMx[k][k2]/dMx[7][9]+ "; ModelMatrixName[k][k2]:= c*t*"+subString+
								";ModelMatrixName[k2][k]:= c*t*"+subString+";");
			}
		}       
	}
	else
	{
		for (k=0; k<19; k=k+1)
		{
			for (k2 = k+1; k2 < 20; k2=k2+1)
			{
				subString = aaNames[k]+aaNames[k2];
				ExecuteCommands("global "+subString +"=" + dMx[k][k2]/dMx[7][9]+ "; ModelMatrixName[k][k2]:= t*"+subString+
								";ModelMatrixName[k2][k]:= t*"+subString+";");
			}
		}       
	}
	
	subString = aaNames[7]+aaNames[9];
	ExecuteCommands(subString +":=1;");
	return 1;
}

/*-------------------------------------------------------------*/
/* VT Matrix */

function PopulateModelMatrixWAG (ModelMatrixName&, EFV)
{
	ModelMatrixName = {20,20};
	ModelMatrixName[14][0]:=t*0.233108;
	ModelMatrixName[0][14]:=t*0.233108;
	ModelMatrixName[11][0]:=t*0.199097;
	ModelMatrixName[0][11]:=t*0.199097;
	ModelMatrixName[11][14]:=t*0.210797;
	ModelMatrixName[14][11]:=t*0.210797;
	ModelMatrixName[2][0]:=t*0.265145;
	ModelMatrixName[0][2]:=t*0.265145;
	ModelMatrixName[2][14]:=t*0.105191;
	ModelMatrixName[14][2]:=t*0.105191;
	ModelMatrixName[2][11]:=t*0.883422;
	ModelMatrixName[11][2]:=t*0.883422;
	ModelMatrixName[1][0]:=t*0.227333;
	ModelMatrixName[0][1]:=t*0.227333;
	ModelMatrixName[1][14]:=t*0.031726;
	ModelMatrixName[14][1]:=t*0.031726;
	ModelMatrixName[1][11]:=t*0.027495;
	ModelMatrixName[11][1]:=t*0.027495;
	ModelMatrixName[1][2]:=t*0.010313;
	ModelMatrixName[2][1]:=t*0.010313;
	ModelMatrixName[13][0]:=t*0.310084;
	ModelMatrixName[0][13]:=t*0.310084;
	ModelMatrixName[13][14]:=t*0.493763;
	ModelMatrixName[14][13]:=t*0.493763;
	ModelMatrixName[13][11]:=t*0.2757;
	ModelMatrixName[11][13]:=t*0.2757;
	ModelMatrixName[13][2]:=t*0.205842;
	ModelMatrixName[2][13]:=t*0.205842;
	ModelMatrixName[13][1]:=t*0.004315;
	ModelMatrixName[1][13]:=t*0.004315;
	ModelMatrixName[3][0]:=t*0.567957;
	ModelMatrixName[0][3]:=t*0.567957;
	ModelMatrixName[3][14]:=t*0.25524;
	ModelMatrixName[14][3]:=t*0.25524;
	ModelMatrixName[3][11]:=t*0.270417;
	ModelMatrixName[11][3]:=t*0.270417;
	ModelMatrixName[3][2]:=t*1.599461;
	ModelMatrixName[2][3]:=t*1.599461;
	ModelMatrixName[3][1]:=t*0.005321;
	ModelMatrixName[1][3]:=t*0.005321;
	ModelMatrixName[3][13]:=t*0.960976;
	ModelMatrixName[13][3]:=t*0.960976;
	ModelMatrixName[5][0]:=t*0.876213;
	ModelMatrixName[0][5]:=t*0.876213;
	ModelMatrixName[5][14]:=t*0.156945;
	ModelMatrixName[14][5]:=t*0.156945;
	ModelMatrixName[5][11]:=t*0.362028;
	ModelMatrixName[11][5]:=t*0.362028;
	ModelMatrixName[5][2]:=t*0.311718;
	ModelMatrixName[2][5]:=t*0.311718;
	ModelMatrixName[5][1]:=t*0.050876;
	ModelMatrixName[1][5]:=t*0.050876;
	ModelMatrixName[5][13]:=t*0.12866;
	ModelMatrixName[13][5]:=t*0.12866;
	ModelMatrixName[5][3]:=t*0.250447;
	ModelMatrixName[3][5]:=t*0.250447;
	ModelMatrixName[6][0]:=t*0.078692;
	ModelMatrixName[0][6]:=t*0.078692;
	ModelMatrixName[6][14]:=t*0.213164;
	ModelMatrixName[14][6]:=t*0.213164;
	ModelMatrixName[6][11]:=t*0.290006;
	ModelMatrixName[11][6]:=t*0.290006;
	ModelMatrixName[6][2]:=t*0.134252;
	ModelMatrixName[2][6]:=t*0.134252;
	ModelMatrixName[6][1]:=t*0.016695;
	ModelMatrixName[1][6]:=t*0.016695;
	ModelMatrixName[6][13]:=t*0.315521;
	ModelMatrixName[13][6]:=t*0.315521;
	ModelMatrixName[6][3]:=t*0.104458;
	ModelMatrixName[3][6]:=t*0.104458;
	ModelMatrixName[6][5]:=t*0.058131;
	ModelMatrixName[5][6]:=t*0.058131;
	ModelMatrixName[7][0]:=t*0.222972;
	ModelMatrixName[0][7]:=t*0.222972;
	ModelMatrixName[7][14]:=t*0.08151;
	ModelMatrixName[14][7]:=t*0.08151;
	ModelMatrixName[7][11]:=t*0.087225;
	ModelMatrixName[11][7]:=t*0.087225;
	ModelMatrixName[7][2]:=t*0.01172;
	ModelMatrixName[2][7]:=t*0.01172;
	ModelMatrixName[7][1]:=t*0.046398;
	ModelMatrixName[1][7]:=t*0.046398;
	ModelMatrixName[7][13]:=t*0.054602;
	ModelMatrixName[13][7]:=t*0.054602;
	ModelMatrixName[7][3]:=t*0.046589;
	ModelMatrixName[3][7]:=t*0.046589;
	ModelMatrixName[7][5]:=t*0.051089;
	ModelMatrixName[5][7]:=t*0.051089;
	ModelMatrixName[7][6]:=t*0.020039;
	ModelMatrixName[6][7]:=t*0.020039;
	ModelMatrixName[9][0]:=t*0.42463;
	ModelMatrixName[0][9]:=t*0.42463;
	ModelMatrixName[9][14]:=t*0.192364;
	ModelMatrixName[14][9]:=t*0.192364;
	ModelMatrixName[9][11]:=t*0.069245;
	ModelMatrixName[11][9]:=t*0.069245;
	ModelMatrixName[9][2]:=t*0.060863;
	ModelMatrixName[2][9]:=t*0.060863;
	ModelMatrixName[9][1]:=t*0.091709;
	ModelMatrixName[1][9]:=t*0.091709;
	ModelMatrixName[9][13]:=t*0.24353;
	ModelMatrixName[13][9]:=t*0.24353;
	ModelMatrixName[9][3]:=t*0.151924;
	ModelMatrixName[3][9]:=t*0.151924;
	ModelMatrixName[9][5]:=t*0.087056;
	ModelMatrixName[5][9]:=t*0.087056;
	ModelMatrixName[9][6]:=t*0.103552;
	ModelMatrixName[6][9]:=t*0.103552;
	ModelMatrixName[9][7]:=t*2.08989;
	ModelMatrixName[7][9]:=t*2.08989;
	ModelMatrixName[8][0]:=t*0.393245;
	ModelMatrixName[0][8]:=t*0.393245;
	ModelMatrixName[8][14]:=t*1.755838;
	ModelMatrixName[14][8]:=t*1.755838;
	ModelMatrixName[8][11]:=t*0.50306;
	ModelMatrixName[11][8]:=t*0.50306;
	ModelMatrixName[8][2]:=t*0.261101;
	ModelMatrixName[2][8]:=t*0.261101;
	ModelMatrixName[8][1]:=t*0.004067;
	ModelMatrixName[1][8]:=t*0.004067;
	ModelMatrixName[8][13]:=t*0.738208;
	ModelMatrixName[13][8]:=t*0.738208;
	ModelMatrixName[8][3]:=t*0.88863;
	ModelMatrixName[3][8]:=t*0.88863;
	ModelMatrixName[8][5]:=t*0.193243;
	ModelMatrixName[5][8]:=t*0.193243;
	ModelMatrixName[8][6]:=t*0.153323;
	ModelMatrixName[6][8]:=t*0.153323;
	ModelMatrixName[8][7]:=t*0.093181;
	ModelMatrixName[7][8]:=t*0.093181;
	ModelMatrixName[8][9]:=t*0.201204;
	ModelMatrixName[9][8]:=t*0.201204;
	ModelMatrixName[10][0]:=t*0.21155;
	ModelMatrixName[0][10]:=t*0.21155;
	ModelMatrixName[10][14]:=t*0.08793;
	ModelMatrixName[14][10]:=t*0.08793;
	ModelMatrixName[10][11]:=t*0.05742;
	ModelMatrixName[11][10]:=t*0.05742;
	ModelMatrixName[10][2]:=t*0.012182;
	ModelMatrixName[2][10]:=t*0.012182;
	ModelMatrixName[10][1]:=t*0.02369;
	ModelMatrixName[1][10]:=t*0.02369;
	ModelMatrixName[10][13]:=t*0.120801;
	ModelMatrixName[13][10]:=t*0.120801;
	ModelMatrixName[10][3]:=t*0.058643;
	ModelMatrixName[3][10]:=t*0.058643;
	ModelMatrixName[10][5]:=t*0.04656;
	ModelMatrixName[5][10]:=t*0.04656;
	ModelMatrixName[10][6]:=t*0.021157;
	ModelMatrixName[6][10]:=t*0.021157;
	ModelMatrixName[10][7]:=t*0.493845;
	ModelMatrixName[7][10]:=t*0.493845;
	ModelMatrixName[10][9]:=t*1.105667;
	ModelMatrixName[9][10]:=t*1.105667;
	ModelMatrixName[10][8]:=t*0.096474;
	ModelMatrixName[8][10]:=t*0.096474;
	ModelMatrixName[4][0]:=t*0.116646;
	ModelMatrixName[0][4]:=t*0.116646;
	ModelMatrixName[4][14]:=t*0.042569;
	ModelMatrixName[14][4]:=t*0.042569;
	ModelMatrixName[4][11]:=t*0.039769;
	ModelMatrixName[11][4]:=t*0.039769;
	ModelMatrixName[4][2]:=t*0.016577;
	ModelMatrixName[2][4]:=t*0.016577;
	ModelMatrixName[4][1]:=t*0.051127;
	ModelMatrixName[1][4]:=t*0.051127;
	ModelMatrixName[4][13]:=t*0.026235;
	ModelMatrixName[13][4]:=t*0.026235;
	ModelMatrixName[4][3]:=t*0.028168;
	ModelMatrixName[3][4]:=t*0.028168;
	ModelMatrixName[4][5]:=t*0.050143;
	ModelMatrixName[5][4]:=t*0.050143;
	ModelMatrixName[4][6]:=t*0.079807;
	ModelMatrixName[6][4]:=t*0.079807;
	ModelMatrixName[4][7]:=t*0.32102;
	ModelMatrixName[7][4]:=t*0.32102;
	ModelMatrixName[4][9]:=t*0.946499;
	ModelMatrixName[9][4]:=t*0.946499;
	ModelMatrixName[4][8]:=t*0.038261;
	ModelMatrixName[8][4]:=t*0.038261;
	ModelMatrixName[4][10]:=t*0.173052;
	ModelMatrixName[10][4]:=t*0.173052;
	ModelMatrixName[12][0]:=t*0.399143;
	ModelMatrixName[0][12]:=t*0.399143;
	ModelMatrixName[12][14]:=t*0.12848;
	ModelMatrixName[14][12]:=t*0.12848;
	ModelMatrixName[12][11]:=t*0.083956;
	ModelMatrixName[11][12]:=t*0.083956;
	ModelMatrixName[12][2]:=t*0.160063;
	ModelMatrixName[2][12]:=t*0.160063;
	ModelMatrixName[12][1]:=t*0.011137;
	ModelMatrixName[1][12]:=t*0.011137;
	ModelMatrixName[12][13]:=t*0.15657;
	ModelMatrixName[13][12]:=t*0.15657;
	ModelMatrixName[12][3]:=t*0.205134;
	ModelMatrixName[3][12]:=t*0.205134;
	ModelMatrixName[12][5]:=t*0.124492;
	ModelMatrixName[5][12]:=t*0.124492;
	ModelMatrixName[12][6]:=t*0.078892;
	ModelMatrixName[6][12]:=t*0.078892;
	ModelMatrixName[12][7]:=t*0.054797;
	ModelMatrixName[7][12]:=t*0.054797;
	ModelMatrixName[12][9]:=t*0.169784;
	ModelMatrixName[9][12]:=t*0.169784;
	ModelMatrixName[12][8]:=t*0.212302;
	ModelMatrixName[8][12]:=t*0.212302;
	ModelMatrixName[12][10]:=t*0.010363;
	ModelMatrixName[10][12]:=t*0.010363;
	ModelMatrixName[12][4]:=t*0.042564;
	ModelMatrixName[4][12]:=t*0.042564;
	ModelMatrixName[15][0]:=t*1.817198;
	ModelMatrixName[0][15]:=t*1.817198;
	ModelMatrixName[15][14]:=t*0.292327;
	ModelMatrixName[14][15]:=t*0.292327;
	ModelMatrixName[15][11]:=t*0.847049;
	ModelMatrixName[11][15]:=t*0.847049;
	ModelMatrixName[15][2]:=t*0.461519;
	ModelMatrixName[2][15]:=t*0.461519;
	ModelMatrixName[15][1]:=t*0.17527;
	ModelMatrixName[1][15]:=t*0.17527;
	ModelMatrixName[15][13]:=t*0.358017;
	ModelMatrixName[13][15]:=t*0.358017;
	ModelMatrixName[15][3]:=t*0.406035;
	ModelMatrixName[3][15]:=t*0.406035;
	ModelMatrixName[15][5]:=t*0.612843;
	ModelMatrixName[5][15]:=t*0.612843;
	ModelMatrixName[15][6]:=t*0.167406;
	ModelMatrixName[6][15]:=t*0.167406;
	ModelMatrixName[15][7]:=t*0.081567;
	ModelMatrixName[7][15]:=t*0.081567;
	ModelMatrixName[15][9]:=t*0.214977;
	ModelMatrixName[9][15]:=t*0.214977;
	ModelMatrixName[15][8]:=t*0.400072;
	ModelMatrixName[8][15]:=t*0.400072;
	ModelMatrixName[15][10]:=t*0.090515;
	ModelMatrixName[10][15]:=t*0.090515;
	ModelMatrixName[15][4]:=t*0.138119;
	ModelMatrixName[4][15]:=t*0.138119;
	ModelMatrixName[15][12]:=t*0.430431;
	ModelMatrixName[12][15]:=t*0.430431;
	ModelMatrixName[16][0]:=t*0.877877;
	ModelMatrixName[0][16]:=t*0.877877;
	ModelMatrixName[16][14]:=t*0.204109;
	ModelMatrixName[14][16]:=t*0.204109;
	ModelMatrixName[16][11]:=t*0.471268;
	ModelMatrixName[11][16]:=t*0.471268;
	ModelMatrixName[16][2]:=t*0.178197;
	ModelMatrixName[2][16]:=t*0.178197;
	ModelMatrixName[16][1]:=t*0.079511;
	ModelMatrixName[1][16]:=t*0.079511;
	ModelMatrixName[16][13]:=t*0.248992;
	ModelMatrixName[13][16]:=t*0.248992;
	ModelMatrixName[16][3]:=t*0.321028;
	ModelMatrixName[3][16]:=t*0.321028;
	ModelMatrixName[16][5]:=t*0.136266;
	ModelMatrixName[5][16]:=t*0.136266;
	ModelMatrixName[16][6]:=t*0.101117;
	ModelMatrixName[6][16]:=t*0.101117;
	ModelMatrixName[16][7]:=t*0.376588;
	ModelMatrixName[7][16]:=t*0.376588;
	ModelMatrixName[16][9]:=t*0.243227;
	ModelMatrixName[9][16]:=t*0.243227;
	ModelMatrixName[16][8]:=t*0.446646;
	ModelMatrixName[8][16]:=t*0.446646;
	ModelMatrixName[16][10]:=t*0.184609;
	ModelMatrixName[10][16]:=t*0.184609;
	ModelMatrixName[16][4]:=t*0.08587;
	ModelMatrixName[4][16]:=t*0.08587;
	ModelMatrixName[16][12]:=t*0.207143;
	ModelMatrixName[12][16]:=t*0.207143;
	ModelMatrixName[16][15]:=t*1.767766;
	ModelMatrixName[15][16]:=t*1.767766;
	ModelMatrixName[18][0]:=t*0.030309;
	ModelMatrixName[0][18]:=t*0.030309;
	ModelMatrixName[18][14]:=t*0.046417;
	ModelMatrixName[14][18]:=t*0.046417;
	ModelMatrixName[18][11]:=t*0.010459;
	ModelMatrixName[11][18]:=t*0.010459;
	ModelMatrixName[18][2]:=t*0.011393;
	ModelMatrixName[2][18]:=t*0.011393;
	ModelMatrixName[18][1]:=t*0.007732;
	ModelMatrixName[1][18]:=t*0.007732;
	ModelMatrixName[18][13]:=t*0.021248;
	ModelMatrixName[13][18]:=t*0.021248;
	ModelMatrixName[18][3]:=t*0.018844;
	ModelMatrixName[3][18]:=t*0.018844;
	ModelMatrixName[18][5]:=t*0.02399;
	ModelMatrixName[5][18]:=t*0.02399;
	ModelMatrixName[18][6]:=t*0.020009;
	ModelMatrixName[6][18]:=t*0.020009;
	ModelMatrixName[18][7]:=t*0.034954;
	ModelMatrixName[7][18]:=t*0.034954;
	ModelMatrixName[18][9]:=t*0.083439;
	ModelMatrixName[9][18]:=t*0.083439;
	ModelMatrixName[18][8]:=t*0.023321;
	ModelMatrixName[8][18]:=t*0.023321;
	ModelMatrixName[18][10]:=t*0.022019;
	ModelMatrixName[10][18]:=t*0.022019;
	ModelMatrixName[18][4]:=t*0.12805;
	ModelMatrixName[4][18]:=t*0.12805;
	ModelMatrixName[18][12]:=t*0.014584;
	ModelMatrixName[12][18]:=t*0.014584;
	ModelMatrixName[18][15]:=t*0.035933;
	ModelMatrixName[15][18]:=t*0.035933;
	ModelMatrixName[18][16]:=t*0.020437;
	ModelMatrixName[16][18]:=t*0.020437;
	ModelMatrixName[19][0]:=t*0.087061;
	ModelMatrixName[0][19]:=t*0.087061;
	ModelMatrixName[19][14]:=t*0.09701;
	ModelMatrixName[14][19]:=t*0.09701;
	ModelMatrixName[19][11]:=t*0.093268;
	ModelMatrixName[11][19]:=t*0.093268;
	ModelMatrixName[19][2]:=t*0.051664;
	ModelMatrixName[2][19]:=t*0.051664;
	ModelMatrixName[19][1]:=t*0.042823;
	ModelMatrixName[1][19]:=t*0.042823;
	ModelMatrixName[19][13]:=t*0.062544;
	ModelMatrixName[13][19]:=t*0.062544;
	ModelMatrixName[19][3]:=t*0.0552;
	ModelMatrixName[3][19]:=t*0.0552;
	ModelMatrixName[19][5]:=t*0.037568;
	ModelMatrixName[5][19]:=t*0.037568;
	ModelMatrixName[19][6]:=t*0.286027;
	ModelMatrixName[6][19]:=t*0.286027;
	ModelMatrixName[19][7]:=t*0.086237;
	ModelMatrixName[7][19]:=t*0.086237;
	ModelMatrixName[19][9]:=t*0.189842;
	ModelMatrixName[9][19]:=t*0.189842;
	ModelMatrixName[19][8]:=t*0.068689;
	ModelMatrixName[8][19]:=t*0.068689;
	ModelMatrixName[19][10]:=t*0.073223;
	ModelMatrixName[10][19]:=t*0.073223;
	ModelMatrixName[19][4]:=t*0.898663;
	ModelMatrixName[4][19]:=t*0.898663;
	ModelMatrixName[19][12]:=t*0.032043;
	ModelMatrixName[12][19]:=t*0.032043;
	ModelMatrixName[19][15]:=t*0.121979;
	ModelMatrixName[15][19]:=t*0.121979;
	ModelMatrixName[19][16]:=t*0.094617;
	ModelMatrixName[16][19]:=t*0.094617;
	ModelMatrixName[19][18]:=t*0.124746;
	ModelMatrixName[18][19]:=t*0.124746;
	ModelMatrixName[17][0]:=t*1.230985;
	ModelMatrixName[0][17]:=t*1.230985;
	ModelMatrixName[17][14]:=t*0.113146;
	ModelMatrixName[14][17]:=t*0.113146;
	ModelMatrixName[17][11]:=t*0.049824;
	ModelMatrixName[11][17]:=t*0.049824;
	ModelMatrixName[17][2]:=t*0.048769;
	ModelMatrixName[2][17]:=t*0.048769;
	ModelMatrixName[17][1]:=t*0.163831;
	ModelMatrixName[1][17]:=t*0.163831;
	ModelMatrixName[17][13]:=t*0.112027;
	ModelMatrixName[13][17]:=t*0.112027;
	ModelMatrixName[17][3]:=t*0.205868;
	ModelMatrixName[3][17]:=t*0.205868;
	ModelMatrixName[17][5]:=t*0.082579;
	ModelMatrixName[5][17]:=t*0.082579;
	ModelMatrixName[17][6]:=t*0.068575;
	ModelMatrixName[6][17]:=t*0.068575;
	ModelMatrixName[17][7]:=t*3.65443;
	ModelMatrixName[7][17]:=t*3.65443;
	ModelMatrixName[17][9]:=t*1.337571;
	ModelMatrixName[9][17]:=t*1.337571;
	ModelMatrixName[17][8]:=t*0.144587;
	ModelMatrixName[8][17]:=t*0.144587;
	ModelMatrixName[17][10]:=t*0.307309;
	ModelMatrixName[10][17]:=t*0.307309;
	ModelMatrixName[17][4]:=t*0.247329;
	ModelMatrixName[4][17]:=t*0.247329;
	ModelMatrixName[17][12]:=t*0.129315;
	ModelMatrixName[12][17]:=t*0.129315;
	ModelMatrixName[17][15]:=t*0.1277;
	ModelMatrixName[15][17]:=t*0.1277;
	ModelMatrixName[17][16]:=t*0.740372;
	ModelMatrixName[16][17]:=t*0.740372;
	ModelMatrixName[17][18]:=t*0.022134;
	ModelMatrixName[18][17]:=t*0.022134;
	ModelMatrixName[17][19]:=t*0.125733;
	ModelMatrixName[19][17]:=t*0.125733;
	return 1;
}

/*-------------------------------------------------------------*/

function countSubstitutionTypes (dataSetNumber)
{
	DataSet dsA									= ReconstructAncestors (WAG_LF);
	DataSetFilter filteredDataA 				= CreateFilter  	   (dsA,1);
	ExecuteCommands	("DataSet dsJoint 			= Combine		       (dsA,ds_"+dataSetNumber+");");
	DataSetFilter filteredDataJ 				= CreateFilter         (dsJoint,1);

	observedCEFV = seqToBranchMap;

	ExecuteCommands ("branchNames = BranchName (WAG_Tree_"+dataSetNumber+",-1);");
	h = Columns (branchNames);

	seqToBranchMap 	= {h, 2};
	/* maps sequence names to branch order in column 1 
	   and the other way around in column 2 */
	   
	ExecuteCommands ("spc=dsf_"+dataSetNumber+".species;usc=dsf_"+dataSetNumber+".unique_sites;tsc=dsf_"+dataSetNumber+".sites;");

	for (k=0; k<spc; k=k+1)
	{
		ExecuteCommands ("GetString (seqName, dsf_"+dataSetNumber+", k);");
		seqToBranchMap[k][0] = -1;
		for (v=0; v<h; v=v+1)
		{
			if (branchNames[v] % seqName)
			{
				seqToBranchMap[k][0] = v;
				seqToBranchMap[v][1] = k;
				break;
			}
		}
	}

	seqToBranchMap[spc][0] = h-1;
	seqToBranchMap[h-1][1] = spc;

	for (k=1; k<filteredDataA.species; k=k+1)
	{
		GetString (seqName, filteredDataA, k);
		seqToBranchMap[spc+k][0] = -1;
		for (v=0; v<h; v=v+1)
		{
			if (branchNames[v] % seqName)
			{
				seqToBranchMap[k+spc][0] = v;
				seqToBranchMap[v][1] = k+spc;
				break;
			}
		}
	}


	/* get residue matrix */

	codonInfo  = {spc, usc};
	codonInfo2 = {filteredDataA.species, filteredDataA.unique_sites};

	ExecuteCommands ("GetDataInfo    (dupInfo, dsf_" + dataSetNumber +");");
	GetDataInfo	    (dupInfoA, filteredDataA);

	matrixTrick  = {1,20}["_MATRIX_ELEMENT_COLUMN_"];
	matrixTrick2 = {1,20}["1"];

	ExecuteCommands ("flatTreeRep	  = Abs (WAG_Tree_"+dataSetNumber+");");

	for (v=0; v<usc;v=v+1)
	{
		for (h=0; h<spc;h=h+1)
		{
			ExecuteCommands ("GetDataInfo (siteInfo, dsf_"+ dataSetNumber+", h, v);");
			_SITE_ES_COUNT = matrixTrick2 * siteInfo; 
			if (_SITE_ES_COUNT[0] == 1)
			{
				siteInfo = matrixTrick * siteInfo;
				codonInfo[h][v] = siteInfo[0];
			}
			else
			{
				codonInfo[h][v] = -1;
			}
		}
	}

	for (v=0; v<filteredDataA.unique_sites;v=v+1)
	{
		for (h=0; h<filteredDataA.species;h=h+1)
		{
			GetDataInfo (siteInfo, filteredDataA, h, v);
			siteInfo2 = matrixTrick2 * siteInfo;
			if (siteInfo2[0] == 1)
			{
				siteInfo = matrixTrick * siteInfo;
				codonInfo2[h][v] = siteInfo[0];
			}
			else
			{
				codonInfo2[h][v] = -1;
			}
		}
	}

	aaSubsCount     = {20,20};

	for (v=0; v<tsc;v=v+1)
	{
		k = spc+1;

		c1 = dupInfoA[v];
		for (h=1; h<filteredDataA.species; h=h+1)
		{
			p1 = seqToBranchMap[k][0];
			p2 = flatTreeRep[p1];
			p2 = seqToBranchMap[p2][1]-spc;
			
			cd1 = codonInfo2[h] [c1];
			cd2 = codonInfo2[p2][c1];
			
			if (cd1 >= 0 && cd2 >= 0)
			{
				aaSubsCount[cd1][cd2] = aaSubsCount[cd1][cd2] + 1;		
			}
			k=k+1;
		}
		
		/* now do the leaves */
		
		observedCEFV = {{0}};
		
		c2 = dupInfo[v];
		for (h=0; h<filteredData.species; h=h+1)
		{
			p1 = seqToBranchMap[h][0];
			p2 = flatTreeRep[p1];
			p2 = seqToBranchMap[p2][1]-spc;
			
			cd2 = codonInfo2[p2][c1];
			cd1 = codonInfo[h] [c2];
			
			if (cd1>=0)
			/* no ambiguities */
			{
				aaSubsCount[cd1][cd2] = aaSubsCount[cd1][cd2] + 1;		
			}	
			else
			/* ambiguities here */
			{
				ExecuteCommands ("GetDataInfo    (ambInfo, dsf_" + dataSetNumber + ", h, c2);");	
				if (Rows(observedCEFV) == 1)
				{
					siteFilter = ""+v;
					DataSetFilter filteredDataSite = CreateFilter (dsJoint,1,siteFilter,"");
					HarvestFrequencies			  (observedCEFV,filteredDataSite,1,1,0);
				}
				
				weightFactor = matrixTrick2*ambInfo;
				if (weightFactor[0]<20)
				{
					ambInfo  	 = ambInfo$observedCEFV;
					
					weightFactor = 0;
					tempMx = -1;
					for (k=0; k<20; k=k+1)
					{
						if (ambInfo[k]>weightFactor)
						{
							weightFactor = ambInfo[k];
							tempMx = k;
						}
					}
					if (tempMx>=0)
					{
						aaSubsCount[tempMx][cd2] = aaSubsCount[tempMx][cd2] + 1;		
					}
				}
			}
		}
	}

	return aaSubsCount;
}

/*-------------------------------------------------------------*/

function reportSubstitutionMatrix (subCountsMatrix, verb)
{
	totalNZ = 0;
	totalSC = 0;
	for (k1 = 0; k1 < 20; k1 = k1 + 1)
	{
		for (k2 = k1 + 1; k2 < 20; k2 = k2 + 1)
		{
			pairCount = subCountsMatrix[k1][k2] + subCountsMatrix[k2][k1];
			if (pairCount > 0)
			{
				totalNZ = totalNZ + 1;
				totalSC = totalSC + pairCount;
				if (verb)
				{
					fprintf (stdout, "\t[", allowedAACharacters[k1], "<->", allowedAACharacters[k2], "] ", Format (pairCount, 8, 0)," substitutions \n");
				}
			}
		}
	}
	return {{totalNZ__,totalSC__}};
}