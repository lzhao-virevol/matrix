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
/* WAG Matrix */

function PopulateModelMatrixWAG (ModelMatrixName&, EFV)
{
	ModelMatrixName = {20,20};
	ModelMatrixName[14][0]:=t*0.425093;
	ModelMatrixName[0][14]:=t*0.425093;
	ModelMatrixName[11][0]:=t*0.276818;
	ModelMatrixName[0][11]:=t*0.276818;
	ModelMatrixName[11][14]:=t*0.751878;
	ModelMatrixName[14][11]:=t*0.751878;
	ModelMatrixName[2][0]:=t*0.395144;
	ModelMatrixName[0][2]:=t*0.395144;
	ModelMatrixName[2][14]:=t*0.123954;
	ModelMatrixName[14][2]:=t*0.123954;
	ModelMatrixName[2][11]:=t*5.076149;
	ModelMatrixName[11][2]:=t*5.076149;
	ModelMatrixName[1][0]:=t*2.489084;
	ModelMatrixName[0][1]:=t*2.489084;
	ModelMatrixName[1][14]:=t*0.534551;
	ModelMatrixName[14][1]:=t*0.534551;
	ModelMatrixName[1][11]:=t*0.528768;
	ModelMatrixName[11][1]:=t*0.528768;
	ModelMatrixName[1][2]:=t*0.062556;
	ModelMatrixName[2][1]:=t*0.062556;
	ModelMatrixName[13][0]:=t*0.969894;
	ModelMatrixName[0][13]:=t*0.969894;
	ModelMatrixName[13][14]:=t*2.807908;
	ModelMatrixName[14][13]:=t*2.807908;
	ModelMatrixName[13][11]:=t*1.695752;
	ModelMatrixName[11][13]:=t*1.695752;
	ModelMatrixName[13][2]:=t*0.523386;
	ModelMatrixName[2][13]:=t*0.523386;
	ModelMatrixName[13][1]:=t*0.084808;
	ModelMatrixName[1][13]:=t*0.084808;
	ModelMatrixName[3][0]:=t*1.038545;
	ModelMatrixName[0][3]:=t*1.038545;
	ModelMatrixName[3][14]:=t*0.363970;
	ModelMatrixName[14][3]:=t*0.363970;
	ModelMatrixName[3][11]:=t*0.541712;
	ModelMatrixName[11][3]:=t*0.541712;
	ModelMatrixName[3][2]:=t*5.243870;
	ModelMatrixName[2][3]:=t*5.243870;
	ModelMatrixName[3][1]:=t*0.003499;
	ModelMatrixName[1][3]:=t*0.003499;
	ModelMatrixName[3][13]:=t*4.128591;
	ModelMatrixName[13][3]:=t*4.128591;
	ModelMatrixName[5][0]:=t*2.066040;
	ModelMatrixName[0][5]:=t*2.066040;
	ModelMatrixName[5][14]:=t*0.390192;
	ModelMatrixName[14][5]:=t*0.390192;
	ModelMatrixName[5][11]:=t*1.437645;
	ModelMatrixName[11][5]:=t*1.437645;
	ModelMatrixName[5][2]:=t*0.844926;
	ModelMatrixName[2][5]:=t*0.844926;
	ModelMatrixName[5][1]:=t*0.569265;
	ModelMatrixName[1][5]:=t*0.569265;
	ModelMatrixName[5][13]:=t*0.267959;
	ModelMatrixName[13][5]:=t*0.267959;
	ModelMatrixName[5][3]:=t*0.348847;
	ModelMatrixName[3][5]:=t*0.348847;
	ModelMatrixName[6][0]:=t*0.358858;
	ModelMatrixName[0][6]:=t*0.358858;
	ModelMatrixName[6][14]:=t*2.426601;
	ModelMatrixName[14][6]:=t*2.426601;
	ModelMatrixName[6][11]:=t*4.509238;
	ModelMatrixName[11][6]:=t*4.509238;
	ModelMatrixName[6][2]:=t*0.927114;
	ModelMatrixName[2][6]:=t*0.927114;
	ModelMatrixName[6][1]:=t*0.640543;
	ModelMatrixName[1][6]:=t*0.640543;
	ModelMatrixName[6][13]:=t*4.813505;
	ModelMatrixName[13][6]:=t*4.813505;
	ModelMatrixName[6][3]:=t*0.423881;
	ModelMatrixName[3][6]:=t*0.423881;
	ModelMatrixName[6][5]:=t*0.311484;
	ModelMatrixName[5][6]:=t*0.311484;
	ModelMatrixName[7][0]:=t*0.149830;
	ModelMatrixName[0][7]:=t*0.149830;
	ModelMatrixName[7][14]:=t*0.126991;
	ModelMatrixName[14][7]:=t*0.126991;
	ModelMatrixName[7][11]:=t*0.191503;
	ModelMatrixName[11][7]:=t*0.191503;
	ModelMatrixName[7][2]:=t*0.010690;
	ModelMatrixName[2][7]:=t*0.010690;
	ModelMatrixName[7][1]:=t*0.320627;
	ModelMatrixName[1][7]:=t*0.320627;
	ModelMatrixName[7][13]:=t*0.072854;
	ModelMatrixName[13][7]:=t*0.072854;
	ModelMatrixName[7][3]:=t*0.044265;
	ModelMatrixName[3][7]:=t*0.044265;
	ModelMatrixName[7][5]:=t*0.008705;
	ModelMatrixName[5][7]:=t*0.008705;
	ModelMatrixName[7][6]:=t*0.108882;
	ModelMatrixName[6][7]:=t*0.108882;
	ModelMatrixName[9][0]:=t*0.395337;
	ModelMatrixName[0][9]:=t*0.395337;
	ModelMatrixName[9][14]:=t*0.301848;
	ModelMatrixName[14][9]:=t*0.301848;
	ModelMatrixName[9][11]:=t*0.068427;
	ModelMatrixName[11][9]:=t*0.068427;
	ModelMatrixName[9][2]:=t*0.015076;
	ModelMatrixName[2][9]:=t*0.015076;
	ModelMatrixName[9][1]:=t*0.594007;
	ModelMatrixName[1][9]:=t*0.594007;
	ModelMatrixName[9][13]:=t*0.582457;
	ModelMatrixName[13][9]:=t*0.582457;
	ModelMatrixName[9][3]:=t*0.069673;
	ModelMatrixName[3][9]:=t*0.069673;
	ModelMatrixName[9][5]:=t*0.044261;
	ModelMatrixName[5][9]:=t*0.044261;
	ModelMatrixName[9][6]:=t*0.366317;
	ModelMatrixName[6][9]:=t*0.366317;
	ModelMatrixName[9][7]:=t*4.145067;
	ModelMatrixName[7][9]:=t*4.145067;
	ModelMatrixName[8][0]:=t*0.536518;
	ModelMatrixName[0][8]:=t*0.536518;
	ModelMatrixName[8][14]:=t*6.326067;
	ModelMatrixName[14][8]:=t*6.326067;
	ModelMatrixName[8][11]:=t*2.145078;
	ModelMatrixName[11][8]:=t*2.145078;
	ModelMatrixName[8][2]:=t*0.282959;
	ModelMatrixName[2][8]:=t*0.282959;
	ModelMatrixName[8][1]:=t*0.013266;
	ModelMatrixName[1][8]:=t*0.013266;
	ModelMatrixName[8][13]:=t*3.234294;
	ModelMatrixName[13][8]:=t*3.234294;
	ModelMatrixName[8][3]:=t*1.807177;
	ModelMatrixName[3][8]:=t*1.807177;
	ModelMatrixName[8][5]:=t*0.296636;
	ModelMatrixName[5][8]:=t*0.296636;
	ModelMatrixName[8][6]:=t*0.697264;
	ModelMatrixName[6][8]:=t*0.697264;
	ModelMatrixName[8][7]:=t*0.159069;
	ModelMatrixName[7][8]:=t*0.159069;
	ModelMatrixName[8][9]:=t*0.137500;
	ModelMatrixName[9][8]:=t*0.137500;
	ModelMatrixName[10][0]:=t*1.124035;
	ModelMatrixName[0][10]:=t*1.124035;
	ModelMatrixName[10][14]:=t*0.484133;
	ModelMatrixName[14][10]:=t*0.484133;
	ModelMatrixName[10][11]:=t*0.371004;
	ModelMatrixName[11][10]:=t*0.371004;
	ModelMatrixName[10][2]:=t*0.025548;
	ModelMatrixName[2][10]:=t*0.025548;
	ModelMatrixName[10][1]:=t*0.893680;
	ModelMatrixName[1][10]:=t*0.893680;
	ModelMatrixName[10][13]:=t*1.672569;
	ModelMatrixName[13][10]:=t*1.672569;
	ModelMatrixName[10][3]:=t*0.173735;
	ModelMatrixName[3][10]:=t*0.173735;
	ModelMatrixName[10][5]:=t*0.139538;
	ModelMatrixName[5][10]:=t*0.139538;
	ModelMatrixName[10][6]:=t*0.442472;
	ModelMatrixName[6][10]:=t*0.442472;
	ModelMatrixName[10][7]:=t*4.273607;
	ModelMatrixName[7][10]:=t*4.273607;
	ModelMatrixName[10][9]:=t*6.312358;
	ModelMatrixName[9][10]:=t*6.312358;
	ModelMatrixName[10][8]:=t*0.656604;
	ModelMatrixName[8][10]:=t*0.656604;
	ModelMatrixName[4][0]:=t*0.253701;
	ModelMatrixName[0][4]:=t*0.253701;
	ModelMatrixName[4][14]:=t*0.052722;
	ModelMatrixName[14][4]:=t*0.052722;
	ModelMatrixName[4][11]:=t*0.089525;
	ModelMatrixName[11][4]:=t*0.089525;
	ModelMatrixName[4][2]:=t*0.017416;
	ModelMatrixName[2][4]:=t*0.017416;
	ModelMatrixName[4][1]:=t*1.105251;
	ModelMatrixName[1][4]:=t*1.105251;
	ModelMatrixName[4][13]:=t*0.035855;
	ModelMatrixName[13][4]:=t*0.035855;
	ModelMatrixName[4][3]:=t*0.018811;
	ModelMatrixName[3][4]:=t*0.018811;
	ModelMatrixName[4][5]:=t*0.089586;
	ModelMatrixName[5][4]:=t*0.089586;
	ModelMatrixName[4][6]:=t*0.682139;
	ModelMatrixName[6][4]:=t*0.682139;
	ModelMatrixName[4][7]:=t*1.112727;
	ModelMatrixName[7][4]:=t*1.112727;
	ModelMatrixName[4][9]:=t*2.592692;
	ModelMatrixName[9][4]:=t*2.592692;
	ModelMatrixName[4][8]:=t*0.023918;
	ModelMatrixName[8][4]:=t*0.023918;
	ModelMatrixName[4][10]:=t*1.798853;
	ModelMatrixName[10][4]:=t*1.798853;
	ModelMatrixName[12][0]:=t*1.177651;
	ModelMatrixName[0][12]:=t*1.177651;
	ModelMatrixName[12][14]:=t*0.332533;
	ModelMatrixName[14][12]:=t*0.332533;
	ModelMatrixName[12][11]:=t*0.161787;
	ModelMatrixName[11][12]:=t*0.161787;
	ModelMatrixName[12][2]:=t*0.394456;
	ModelMatrixName[2][12]:=t*0.394456;
	ModelMatrixName[12][1]:=t*0.075382;
	ModelMatrixName[1][12]:=t*0.075382;
	ModelMatrixName[12][13]:=t*0.624294;
	ModelMatrixName[13][12]:=t*0.624294;
	ModelMatrixName[12][3]:=t*0.419409;
	ModelMatrixName[3][12]:=t*0.419409;
	ModelMatrixName[12][5]:=t*0.196961;
	ModelMatrixName[5][12]:=t*0.196961;
	ModelMatrixName[12][6]:=t*0.508851;
	ModelMatrixName[6][12]:=t*0.508851;
	ModelMatrixName[12][7]:=t*0.078281;
	ModelMatrixName[7][12]:=t*0.078281;
	ModelMatrixName[12][9]:=t*0.249060;
	ModelMatrixName[9][12]:=t*0.249060;
	ModelMatrixName[12][8]:=t*0.390322;
	ModelMatrixName[8][12]:=t*0.390322;
	ModelMatrixName[12][10]:=t*0.099849;
	ModelMatrixName[10][12]:=t*0.099849;
	ModelMatrixName[12][4]:=t*0.094464;
	ModelMatrixName[4][12]:=t*0.094464;
	ModelMatrixName[15][0]:=t*4.727182;
	ModelMatrixName[0][15]:=t*4.727182;
	ModelMatrixName[15][14]:=t*0.858151;
	ModelMatrixName[14][15]:=t*0.858151;
	ModelMatrixName[15][11]:=t*4.008358;
	ModelMatrixName[11][15]:=t*4.008358;
	ModelMatrixName[15][2]:=t*1.240275;
	ModelMatrixName[2][15]:=t*1.240275;
	ModelMatrixName[15][1]:=t*2.784478;
	ModelMatrixName[1][15]:=t*2.784478;
	ModelMatrixName[15][13]:=t*1.223828;
	ModelMatrixName[13][15]:=t*1.223828;
	ModelMatrixName[15][3]:=t*0.611973;
	ModelMatrixName[3][15]:=t*0.611973;
	ModelMatrixName[15][5]:=t*1.739990;
	ModelMatrixName[5][15]:=t*1.739990;
	ModelMatrixName[15][6]:=t*0.990012;
	ModelMatrixName[6][15]:=t*0.990012;
	ModelMatrixName[15][7]:=t*0.064105;
	ModelMatrixName[7][15]:=t*0.064105;
	ModelMatrixName[15][9]:=t*0.182287;
	ModelMatrixName[9][15]:=t*0.182287;
	ModelMatrixName[15][8]:=t*0.748683;
	ModelMatrixName[8][15]:=t*0.748683;
	ModelMatrixName[15][10]:=t*0.346960;
	ModelMatrixName[10][15]:=t*0.346960;
	ModelMatrixName[15][4]:=t*0.361819;
	ModelMatrixName[4][15]:=t*0.361819;
	ModelMatrixName[15][12]:=t*1.338132;
	ModelMatrixName[12][15]:=t*1.338132;
	ModelMatrixName[16][0]:=t*2.139501;
	ModelMatrixName[0][16]:=t*2.139501;
	ModelMatrixName[16][14]:=t*0.578987;
	ModelMatrixName[14][16]:=t*0.578987;
	ModelMatrixName[16][11]:=t*2.000679;
	ModelMatrixName[11][16]:=t*2.000679;
	ModelMatrixName[16][2]:=t*0.425860;
	ModelMatrixName[2][16]:=t*0.425860;
	ModelMatrixName[16][1]:=t*1.143480;
	ModelMatrixName[1][16]:=t*1.143480;
	ModelMatrixName[16][13]:=t*1.080136;
	ModelMatrixName[13][16]:=t*1.080136;
	ModelMatrixName[16][3]:=t*0.604545;
	ModelMatrixName[3][16]:=t*0.604545;
	ModelMatrixName[16][5]:=t*0.129836;
	ModelMatrixName[5][16]:=t*0.129836;
	ModelMatrixName[16][6]:=t*0.584262;
	ModelMatrixName[6][16]:=t*0.584262;
	ModelMatrixName[16][7]:=t*1.033739;
	ModelMatrixName[7][16]:=t*1.033739;
	ModelMatrixName[16][9]:=t*0.302936;
	ModelMatrixName[9][16]:=t*0.302936;
	ModelMatrixName[16][8]:=t*1.136863;
	ModelMatrixName[8][16]:=t*1.136863;
	ModelMatrixName[16][10]:=t*2.020366;
	ModelMatrixName[10][16]:=t*2.020366;
	ModelMatrixName[16][4]:=t*0.165001;
	ModelMatrixName[4][16]:=t*0.165001;
	ModelMatrixName[16][12]:=t*0.571468;
	ModelMatrixName[12][16]:=t*0.571468;
	ModelMatrixName[16][15]:=t*6.472279;
	ModelMatrixName[15][16]:=t*6.472279;
	ModelMatrixName[18][0]:=t*0.180717;
	ModelMatrixName[0][18]:=t*0.180717;
	ModelMatrixName[18][14]:=t*0.593607;
	ModelMatrixName[14][18]:=t*0.593607;
	ModelMatrixName[18][11]:=t*0.045376;
	ModelMatrixName[11][18]:=t*0.045376;
	ModelMatrixName[18][2]:=t*0.029890;
	ModelMatrixName[2][18]:=t*0.029890;
	ModelMatrixName[18][1]:=t*0.670128;
	ModelMatrixName[1][18]:=t*0.670128;
	ModelMatrixName[18][13]:=t*0.236199;
	ModelMatrixName[13][18]:=t*0.236199;
	ModelMatrixName[18][3]:=t*0.077852;
	ModelMatrixName[3][18]:=t*0.077852;
	ModelMatrixName[18][5]:=t*0.268491;
	ModelMatrixName[5][18]:=t*0.268491;
	ModelMatrixName[18][6]:=t*0.597054;
	ModelMatrixName[6][18]:=t*0.597054;
	ModelMatrixName[18][7]:=t*0.111660;
	ModelMatrixName[7][18]:=t*0.111660;
	ModelMatrixName[18][9]:=t*0.619632;
	ModelMatrixName[9][18]:=t*0.619632;
	ModelMatrixName[18][8]:=t*0.049906;
	ModelMatrixName[8][18]:=t*0.049906;
	ModelMatrixName[18][10]:=t*0.696175;
	ModelMatrixName[10][18]:=t*0.696175;
	ModelMatrixName[18][4]:=t*2.457121;
	ModelMatrixName[4][18]:=t*2.457121;
	ModelMatrixName[18][12]:=t*0.095131;
	ModelMatrixName[12][18]:=t*0.095131;
	ModelMatrixName[18][15]:=t*0.248862;
	ModelMatrixName[15][18]:=t*0.248862;
	ModelMatrixName[18][16]:=t*0.140825;
	ModelMatrixName[16][18]:=t*0.140825;
	ModelMatrixName[19][0]:=t*0.218959;
	ModelMatrixName[0][19]:=t*0.218959;
	ModelMatrixName[19][14]:=t*0.314440;
	ModelMatrixName[14][19]:=t*0.314440;
	ModelMatrixName[19][11]:=t*0.612025;
	ModelMatrixName[11][19]:=t*0.612025;
	ModelMatrixName[19][2]:=t*0.135107;
	ModelMatrixName[2][19]:=t*0.135107;
	ModelMatrixName[19][1]:=t*1.165532;
	ModelMatrixName[1][19]:=t*1.165532;
	ModelMatrixName[19][13]:=t*0.257336;
	ModelMatrixName[13][19]:=t*0.257336;
	ModelMatrixName[19][3]:=t*0.120037;
	ModelMatrixName[3][19]:=t*0.120037;
	ModelMatrixName[19][5]:=t*0.054679;
	ModelMatrixName[5][19]:=t*0.054679;
	ModelMatrixName[19][6]:=t*5.306834;
	ModelMatrixName[6][19]:=t*5.306834;
	ModelMatrixName[19][7]:=t*0.232523;
	ModelMatrixName[7][19]:=t*0.232523;
	ModelMatrixName[19][9]:=t*0.299648;
	ModelMatrixName[9][19]:=t*0.299648;
	ModelMatrixName[19][8]:=t*0.131932;
	ModelMatrixName[8][19]:=t*0.131932;
	ModelMatrixName[19][10]:=t*0.481306;
	ModelMatrixName[10][19]:=t*0.481306;
	ModelMatrixName[19][4]:=t*7.803902;
	ModelMatrixName[4][19]:=t*7.803902;
	ModelMatrixName[19][12]:=t*0.089613;
	ModelMatrixName[12][19]:=t*0.089613;
	ModelMatrixName[19][15]:=t*0.400547;
	ModelMatrixName[15][19]:=t*0.400547;
	ModelMatrixName[19][16]:=t*0.245841;
	ModelMatrixName[16][19]:=t*0.245841;
	ModelMatrixName[19][18]:=t*3.151815;
	ModelMatrixName[18][19]:=t*3.151815;
	ModelMatrixName[17][0]:=t*2.547870;
	ModelMatrixName[0][17]:=t*2.547870;
	ModelMatrixName[17][14]:=t*0.170887;
	ModelMatrixName[14][17]:=t*0.170887;
	ModelMatrixName[17][11]:=t*0.083688;
	ModelMatrixName[11][17]:=t*0.083688;
	ModelMatrixName[17][2]:=t*0.037967;
	ModelMatrixName[2][17]:=t*0.037967;
	ModelMatrixName[17][1]:=t*1.959291;
	ModelMatrixName[1][17]:=t*1.959291;
	ModelMatrixName[17][13]:=t*0.210332;
	ModelMatrixName[13][17]:=t*0.210332;
	ModelMatrixName[17][3]:=t*0.245034;
	ModelMatrixName[3][17]:=t*0.245034;
	ModelMatrixName[17][5]:=t*0.076701;
	ModelMatrixName[5][17]:=t*0.076701;
	ModelMatrixName[17][6]:=t*0.119013;
	ModelMatrixName[6][17]:=t*0.119013;
	ModelMatrixName[17][7]:=t*10.649107;
	ModelMatrixName[7][17]:=t*10.649107;
	ModelMatrixName[17][9]:=t*1.702745;
	ModelMatrixName[9][17]:=t*1.702745;
	ModelMatrixName[17][8]:=t*0.185202;
	ModelMatrixName[8][17]:=t*0.185202;
	ModelMatrixName[17][10]:=t*1.898718;
	ModelMatrixName[10][17]:=t*1.898718;
	ModelMatrixName[17][4]:=t*0.654683;
	ModelMatrixName[4][17]:=t*0.654683;
	ModelMatrixName[17][12]:=t*0.296501;
	ModelMatrixName[12][17]:=t*0.296501;
	ModelMatrixName[17][15]:=t*0.098369;
	ModelMatrixName[15][17]:=t*0.098369;
	ModelMatrixName[17][16]:=t*2.188158;
	ModelMatrixName[16][17]:=t*2.188158;
	ModelMatrixName[17][18]:=t*0.189510;
	ModelMatrixName[18][17]:=t*0.189510;
	ModelMatrixName[17][19]:=t*0.249313;
	ModelMatrixName[19][17]:=t*0.249313;
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