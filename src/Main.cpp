/*
 * Main.cpp
 *
 *  Created on: Jan 14, 2017
 *      Author: Julia
 */
#include<iostream>
#include<dirent.h>
#include <sys/stat.h>
#include <RInside.h>
#include <fstream>
#include<cstring>
using namespace std;
#include "SEinventory.h"
#include "CSVRow.h"
#include "CSVIterator.h"
#include "Summator.h"
#include "Mean.h"
#include "Variance.h"
#include "Correlation.h"
#include "tTest.h"
#include "LinearRegression.h"
#include "ZScoreMethod.h"



#define help 911
#define bufOverflow 999
#define notFound 404
#define parseError 123
#define paramError 500

void help_message(int err_message)
{
	cerr<<endl;
	cerr<<endl;
	cerr<<endl;
	cerr<<"SE_Compare: Written by Julia D. Sommer"<<endl; 
	cerr<<endl;
	cerr<<"Usage: SE_Compare [Options]"<<endl;
	cerr<<endl;

	cerr<<"Options:"<<endl;
	cerr<<"--rTrain: [Required]: right side skeletal elements ref. population"<<endl;
	cerr<<"--lTrain: [Required]: left side skeletal elements ref. population"<<endl;
	cerr<<"--rClass: right side skeletal elements population to be tested"<<endl;
	cerr<<"--lClass: left side skeletal elements population to be tested"<<endl;
	cerr<<"--uweightedZ: The unweighted Zscore method for generating p-values is used"<<endl;
	cerr<<"--effectSizeZ: The weighted (effect size) Zscore method for generating p-values is used"<<endl;
	cerr<<"--standardErrorZ: The weighted (standard error) Zscore method for generating p-values is used"<<endl;
	cerr<<"--tTest: The t-test for summed measurements to generate p-values is used"<<endl;
	cerr<<"--abs_tTest: The absolute value t-test for summed measurements to generate p-values is used"<<endl;
	cerr<<"--tTest_wMean: The t-test (mean-centered) for summed measurements to generate p-values is used"<<endl;
	cerr<<"--LOOCV: Leave one out cross validation (Performed on training data)"<<endl;
	cerr<<"--descLen: Maximum skeletal element name length (default 50)"<<endl;
	cerr<<"--header: The first line of the input files is the header line (default true)"<<endl;	
	
	cerr<<endl;
        cerr<<endl;
        cerr<<endl;
	cerr<<"Example Usage:"<<endl;
	cerr<<endl;

	cerr<<"Perform Leave-One-Out Crossvalidation on Reference population set using unweighted Z method"<<endl; 
	cerr<<"./SE_Compare --rTrain fileR --lTrain fileL --uweightedZ TRUE"<<endl;
	cerr<<endl;	
	
	cerr<<"Use t-test method to obtain p-values for test data using a reference population set"<<endl;
	cerr<<"./SE_Compare --rTrain fileR --lTrain fileL --rClass file2R --lClass file2L --tTest TRUE"<<endl;

	
	exit(911);
}

int main(int argc, char * argv[])
{
	RInside R;

	//Param options and defaults

	string rightFileTrain = "noFile"; 	//right side skeletal elements ref. population
	string leftFileTrain = "noFile";	//left side skeletal elements ref. population
	string rightFileClass = "noFile";	//right side skeletal elements population to be tested
	string leftFileClass = "noFile";	//left side skeletal elements population to be tested


	bool unweightedZscore = false;		//The unweighted Zscore method for generating p-values is used.
	bool weightedEffectZscore = false;	//The weighted (effect size) Zscore method for generating p-values is used.
	bool weightedStdErrZscore = false;	//The weighted (standard error) Zscore method for generating p-values is used.
	bool sumtTest = false;			//The t-test for summed measurements to generate p-values is used.
	bool abstTest = false;			//The absolute value t-test for summed measurements to generate p-values is used.
	bool meanSumtTest = false; 		//The t-test (mean-centered) for summed measurements to generate p-values is used.

	int minSig = 1;

	bool LOOCV = false;			//Leave one out cross validation

	bool KFCV = false;			//K-fold cross validation.If both the training and classification
	//files are available than this parameter is disabled.

	//Other Parameters
	int seNameLen = 50;                    //Maximum skeletal element length
	bool header = TRUE;					 //Are there headers on the input files: default true


	//user options
	string rightFileTrainS = "--rTrain";
	string leftFileTrainS = "--lTrain";
	string rightFileClassS = "--rClass";
	string leftFileClassS = "--lClass";

	string unweightedZscoreS = "--uweightedZ";
	string weightedEffectZscoreS = "--effectSizeZ";
	string weightedStdErrZscoreS = "--standardErrorZ";
	string sumtTestS = "--tTest";
	string abstTestS = "--abs_tTest";
	string meanSumtTestS = "--tTest_wMean";

	string LOOCVS = "--LOOCV";
	string KFCVS = "--KFCV";

	string seNameLenS = "--descLen";
	string headerS = "--header";

	for(int i = 1; i < argc; i+=2)
	{

		//The input files
		if(argv[i] == rightFileTrainS)
		{
			rightFileTrain = argv[i+1];
		}

		if(argv[i] == leftFileTrainS)
		{
			leftFileTrain = argv[i+1];
		}

		if(argv[i] == rightFileClassS)
		{
			rightFileClass = argv[i+1];
		}

		if(argv[i] == leftFileClassS)
		{
			leftFileClass = argv[i+1];
		}


		//Methods
		if(argv[i] == unweightedZscoreS)
		{
			string TorF = argv[i+1];
			if(TorF.compare("FALSE") == 0 || TorF.compare("F") == 0)
			{
				unweightedZscore = false;
			}else{
				if(TorF.compare("TRUE") == 0 || TorF.compare("T") == 0)
				{
					unweightedZscore = true;
				}
			}
		}

		if(argv[i] == weightedEffectZscoreS)
		{
			string TorF = argv[i+1];
			if(TorF.compare("FALSE") == 0 || TorF.compare("F") == 0)
			{
				weightedEffectZscore = false;
			}else{
				if(TorF.compare("TRUE") == 0 || TorF.compare("T") == 0)
				{
					weightedEffectZscore = true;
				}
			}
		}

		if(argv[i] == weightedStdErrZscoreS)
		{
			string TorF = argv[i+1];
			if(TorF.compare("FALSE") == 0 || TorF.compare("F") == 0)
			{
				weightedStdErrZscore = false;
			}else{
				if(TorF.compare("TRUE") == 0 || TorF.compare("T") == 0)
				{
					weightedStdErrZscore = true;
				}
			}
		}

		if(argv[i] == sumtTestS)
		{
			string TorF = argv[i+1];
			if(TorF.compare("FALSE") == 0 || TorF.compare("F") == 0)
			{
				sumtTest = false;
			}else{
				if(TorF.compare("TRUE") == 0 || TorF.compare("T") == 0)
				{
					sumtTest = true;
				}
			}
		}

		if(argv[i] == abstTestS)
		{
			string TorF = argv[i+1];
			if(TorF.compare("FALSE") == 0 || TorF.compare("F") == 0)
			{
				abstTest = false;
			}else{
				if(TorF.compare("TRUE") == 0 || TorF.compare("T") == 0)
				{
					abstTest = true;
				}
			}
		}			


		if(argv[i] == meanSumtTestS)
		{
			string TorF = argv[i+1];
			if(TorF.compare("FALSE") == 0 || TorF.compare("F") == 0)
			{
				meanSumtTest = false;
			}else{
				if(TorF.compare("TRUE") == 0 || TorF.compare("T") == 0)
				{
					meanSumtTest = true;
				}
			}
		}

		//Cross-Validation methods
		if(argv[i] == LOOCVS)
		{
			string TorF = argv[i+1];
			if(TorF.compare("FALSE") == 0 || TorF.compare("F") == 0)
			{
				LOOCV = false;
			}else{
				if(TorF.compare("TRUE") == 0 || TorF.compare("T") == 0)
				{
					LOOCV = true;
				}
			}
		}

		if(argv[i] == KFCVS)
		{
			string TorF = argv[i+1];
			if(TorF.compare("FALSE") == 0 || TorF.compare("F") == 0)
			{
				KFCV = false;
			}else{
				if(TorF.compare("TRUE") == 0 || TorF.compare("T") == 0)
				{
					KFCV = true;
				}
			}
		}

		//Other Parameters
		if(argv[i] == seNameLenS)
		{
			seNameLen =  atoi(argv[i+1]);
		}

		if(argv[i] == headerS)
		{
			string TorF = argv[i+1];
			if(TorF.compare("FALSE") == 0 || TorF.compare("F") == 0)
			{
				header = false;
			}else{
				if(TorF.compare("TRUE") == 0 || TorF.compare("T") == 0)
				{
					header = true;
				}
			}
		}

	}

	//PARAMETER CHECKING//
	//User forgot a file, oops
	if(argc == 1 ||  rightFileTrain.compare("noFile") == 0 || leftFileTrain.compare("noFile") == 0)
	{
		help_message(notFound);
	}

	if(rightFileClass.compare("noFile") != 0 || leftFileClass.compare("noFile") != 0 )
	{
		LOOCV = false; KFCV = false;
	}

	if(argc == 1 ||  rightFileClass.compare("noFile") != 0 || leftFileClass.compare("noFile") == 0)
	{
		help_message(notFound);
	}

	if(argc == 1 ||  rightFileClass.compare("noFile") == 0 || leftFileClass.compare("noFile") != 0)
	{
		help_message(notFound);
	}


	ifstream iFile (rightFileTrain.c_str());

	if(header)
	{
		string line;
		getline(iFile, line);
	}

	int numRecordsRT = 0; int numMeasurementsRT = 0;
	string elementTypeRT = ""; string elementSideRT = "";
	for(CSVIterator loop(iFile); loop != CSVIterator(); ++loop)
	{
		numMeasurementsRT = (*loop).size() - 3;
		elementSideRT = (*loop)[1];
		elementTypeRT = (*loop)[2];

		numRecordsRT++;
	}

	iFile.close();

	SEinventory rightTrain(numMeasurementsRT, numRecordsRT, seNameLen, elementSideRT, elementTypeRT);

	float * measures = new float[numMeasurementsRT]; string id;

	iFile.open(rightFileTrain.c_str());

	if(header)
	{
		string line;
		getline(iFile, line);
	}

	for(CSVIterator loop(iFile); loop != CSVIterator(); ++loop)
	{
		for(int i = 0; i < (*loop).size()-3; i++)
		{
			measures[i] = atof((*loop)[i+3].c_str());
		}

		id = (*loop)[0];
		rightTrain.addMeasurements(id, measures);
	}

	iFile.close();

	iFile.open(leftFileTrain.c_str());

	if(header)
	{
		string line;
		getline(iFile, line);
	}

	int numRecordsLT = 0; int numMeasurementsLT = 0;
	string elementTypeLT = ""; string elementSideLT = "";
	for(CSVIterator loop(iFile); loop != CSVIterator(); ++loop)
	{
		numMeasurementsLT = (*loop).size() - 3;
		elementSideLT = (*loop)[1];
		elementTypeLT = (*loop)[2];
		numRecordsLT++;
	}

	iFile.close();

	if(numMeasurementsLT != numMeasurementsRT)
	{
		delete [] measures;
		cerr<<"Number of measurements in right and left ";
		cerr<<"training files do not match"<<endl;
		help_message(parseError);
	}

	if(numRecordsLT != numRecordsRT)
	{
		delete [] measures;

		cerr<<"Number of records in right and left ";
		cerr<<"training files do not match"<<endl;
		help_message(parseError);
	}


	SEinventory leftTrain(numMeasurementsLT, numRecordsLT, seNameLen, elementSideLT, elementTypeLT);

	iFile.open(leftFileTrain.c_str());

	if(header)
	{
		string line;
		getline(iFile, line);
	}

	for(CSVIterator loop(iFile); loop != CSVIterator(); ++loop)
	{
		for(int i = 0; i < (*loop).size()-3; i++)
		{
			measures[i] = atof((*loop)[i+3].c_str());
		}

		id = (*loop)[0];
		leftTrain.addMeasurements(id, measures);
	}

	iFile.close();

	int numRecordsTrain = (numRecordsRT+numRecordsLT)/2;
	int numMeasurementsTrain = (numMeasurementsRT + numMeasurementsLT)/2;

	double ** Dvals = new double * [numMeasurementsTrain];
	for(int i = 0; i < numMeasurementsTrain; i++)
	{
		Dvals[i] = new double[numRecordsTrain];
	}

	for(int i = 0; i < numMeasurementsTrain; i++)
	{
		for(int j = 0; j < numRecordsTrain; j++)
		{
			Dvals[i][j] = leftTrain.getRecordMeas(j,i) - rightTrain.getRecordMeas(j,i); //CHANGED
		}
	}

	double means[numMeasurementsTrain]; Mean<double> mean;
	for(int i = 0; i < numMeasurementsTrain; i++)
	{
		means[i] = mean.getMean(Dvals[i], numRecordsTrain);
	}

	double vars[numMeasurementsTrain]; Variance<double> var;
	for(int i = 0; i < numMeasurementsTrain; i++)
	{
		vars[i] = var.getVar(Dvals[i], numRecordsTrain, means[i]);
	}


	double ** cors = new double * [numMeasurementsTrain];
	for(int i = 0; i < numMeasurementsTrain; i++)
		cors[i] = new double[numMeasurementsTrain];


	if(rightFileClass.compare("noFile") != 0 || leftFileClass.compare("noFile") != 0 )
	{
		iFile.open(rightFileClass.c_str());

		if(header)
		{
			string line;
			getline(iFile, line);
		}

		int numRecordsRC = 0; int numMeasurementsRC = 0;
		string elementTypeRC = ""; string elementSideRC = "";
		for(CSVIterator loop(iFile); loop != CSVIterator(); ++loop)
		{
			numMeasurementsRC = (*loop).size() - 3;
			elementSideRC = (*loop)[1];
			elementTypeRC = (*loop)[2];
			numRecordsRC++;
		}

		iFile.close();

		SEinventory rightClass(numMeasurementsRC, numRecordsRC, seNameLen, elementSideRC, elementTypeRC);

		float * measures = new float[numMeasurementsRC]; string id;

		iFile.open(rightFileClass.c_str());

		if(header)
		{
			string line;
			getline(iFile, line);
		}

		for(CSVIterator loop(iFile); loop != CSVIterator(); ++loop)
		{
			for(int i = 0; i < (*loop).size()-3; i++)
			{
				measures[i] = atof((*loop)[i+3].c_str());
			}

			id = (*loop)[0];
			rightClass.addMeasurements(id, measures);
		}

		iFile.close();

		iFile.open(leftFileClass.c_str());

		if(header)
		{
			string line;
			getline(iFile, line);
		}

		int numRecordsLC = 0; int numMeasurementsLC = 0;
		string elementTypeLC = ""; string elementSideLC = "";
		for(CSVIterator loop(iFile); loop != CSVIterator(); ++loop)
		{
			numMeasurementsLC = (*loop).size() - 3;
			elementSideLC = (*loop)[1];
			elementTypeLC = (*loop)[2];
			numRecordsLC++;
		}

		iFile.close();

		if(numMeasurementsLC != numMeasurementsRC)
		{
			delete [] measures;
			cerr<<"Number of measurements in right and left ";
			cerr<<"classification files do not match"<<endl;
			help_message(parseError);
		}

		if(numRecordsLC != numRecordsRC)
		{
			delete [] measures;

			cerr<<"Number of records in right and left ";
			cerr<<"classification files do not match"<<endl;
			help_message(parseError);
		}

		SEinventory leftClass(numMeasurementsLC, numRecordsLC, seNameLen, elementSideLC, elementTypeLC);

		iFile.open(leftFileClass.c_str());

		if(header)
		{
			string line;
			getline(iFile, line);
		}

		for(CSVIterator loop(iFile); loop != CSVIterator(); ++loop)
		{
			for(int i = 0; i < (*loop).size()-3; i++)
			{
				measures[i] = atof((*loop)[i+3].c_str());
			}

			id = (*loop)[0];
			leftClass.addMeasurements(id, measures);
		}

		iFile.close();

		int numRecordsClass = (numRecordsRC+numRecordsLC)/2;
		int numMeasurementsClass = (numMeasurementsRC + numMeasurementsLC)/2;

		if(unweightedZscore && minSig == 1)
		{
			double zMean[numMeasurementsTrain];

			double ** crossDvals = new double * [numMeasurementsClass];
			double ** crossPvals = new double * [numMeasurementsClass];
			double * result = new double [numRecordsTrain];
			double ** Zscores = new double * [numMeasurementsTrain];
			double ** pVals = new double * [numMeasurementsTrain];

			for(int i = 0; i < numMeasurementsTrain; i++)
			{
				crossDvals[i] = new double [numRecordsClass];
				crossPvals[i] = new double [numRecordsClass];
				Zscores[i] = new double [numRecordsTrain];
				pVals[i] = new double  [numRecordsTrain];
			}


			for(int k = 0; k < numMeasurementsTrain; k++)
			{
				tTest<double> Tstats;
				Tstats.getPvals(Dvals[k], pVals[k], means[k], vars[k], numRecordsTrain);
			}

			ZScoreMethod zScoreMethod;
			zScoreMethod.OneTpValsToZ(pVals, Zscores, numRecordsTrain, numMeasurementsTrain);

			for(int k = 0; k < numMeasurementsTrain; k++)
			{
				zMean[k] = mean.getMean(Zscores[k], numRecordsTrain);
			}

			for(int k = 0; k < numMeasurementsTrain; k++)
			{
				for(int h = k+1; h < numMeasurementsTrain; h++)
				{
					Correlation<double> correlation;
					cors[k][h] = correlation.getCor(Zscores[k], Zscores[h], numRecordsTrain, zMean[k], zMean[h]);
				}
			}


			for(int i = 0; i < numRecordsClass; i++)
			{
				double crossZvals[numMeasurementsTrain]; double weights[numMeasurementsTrain];
				for(int k = 0; k < numMeasurementsTrain; k++)
				{
					crossDvals[k][i] = leftClass.getRecordMeas(i,k) - rightClass.getRecordMeas(i,k);

					tTest<double> Tstats;

					double tStat = Tstats.getTStatistic(crossDvals[k], i, means[k], vars[k]);
					crossPvals[k][i] = Tstats.twoTailedPVal(tStat, numRecordsTrain);

					crossZvals[k] = zScoreMethod.getZVal1T(crossPvals[k][i]); //zScoreMethod.getZVal2T(crossPvals[k][j]);
					weights[k] = 1;

				}

				double finalP = zScoreMethod.combineZ(crossZvals, weights, cors, numMeasurementsTrain);
				cout<<leftClass.getRecordId(i)<<","<<rightClass.getRecordId(i)<<","<<finalP<<endl;
			}
		}


		if(weightedEffectZscore && minSig == 1)
		{
			double zMean[numMeasurementsTrain];

			double ** crossDvals = new double * [numMeasurementsClass];
			double ** crossPvals = new double * [numMeasurementsClass];
			double * result = new double [numRecordsTrain];
			double ** Zscores = new double * [numMeasurementsTrain];
			double ** pVals = new double * [numMeasurementsTrain];

			for(int i = 0; i < numMeasurementsTrain; i++)
			{
				crossDvals[i] = new double [numRecordsClass];
				crossPvals[i] = new double [numRecordsClass];
				Zscores[i] = new double [numRecordsTrain];
				pVals[i] = new double  [numRecordsTrain];
			}


			for(int k = 0; k < numMeasurementsTrain; k++)
			{
				tTest<double> Tstats;
				Tstats.getPvals(Dvals[k], pVals[k], means[k], vars[k], numRecordsTrain);
			}

			ZScoreMethod zScoreMethod;
			zScoreMethod.OneTpValsToZ(pVals, Zscores, numRecordsTrain, numMeasurementsTrain);

			for(int k = 0; k < numMeasurementsTrain; k++)
			{
				zMean[k] = mean.getMean(Zscores[k], numRecordsTrain);
			}

			for(int k = 0; k < numMeasurementsTrain; k++)
			{
				for(int h = k+1; h < numMeasurementsTrain; h++)
				{
					Correlation<double> correlation;
					cors[k][h] = correlation.getCor(Zscores[k], Zscores[h], numRecordsTrain, zMean[k], zMean[h]);
				}
			}


			for(int i = 0; i < numRecordsClass; i++)
			{
				double crossZvals[numMeasurementsTrain]; double weights[numMeasurementsTrain];
				for(int k = 0; k < numMeasurementsTrain; k++)
				{
					crossDvals[k][i] = leftClass.getRecordMeas(i,k) - rightClass.getRecordMeas(i,k);

					tTest<double> Tstats;

					double tStat = Tstats.getTStatistic(crossDvals[k], i, means[k], vars[k]);
					crossPvals[k][i] = Tstats.twoTailedPVal(tStat, numRecordsTrain);

					crossZvals[k] = zScoreMethod.getZVal1T(crossPvals[k][i]); //zScoreMethod.getZVal2T(crossPvals[k][j]);
					weights[k] =  abs(tStat);

				}

				double finalP = zScoreMethod.combineZ(crossZvals, weights, cors, numMeasurementsTrain);
				cout<<leftClass.getRecordId(i)<<","<<rightClass.getRecordId(i)<<","<<finalP<<endl;
			}
		}

		if(weightedStdErrZscore && minSig == 1)
		{
			double zMean[numMeasurementsTrain];

			double ** crossDvals = new double * [numMeasurementsClass];
			double ** crossPvals = new double * [numMeasurementsClass];
			double * result = new double [numRecordsTrain];
			double ** Zscores = new double * [numMeasurementsTrain];
			double ** pVals = new double * [numMeasurementsTrain];

			for(int i = 0; i < numMeasurementsTrain; i++)
			{
				crossDvals[i] = new double [numRecordsClass];
				crossPvals[i] = new double [numRecordsClass];
				Zscores[i] = new double [numRecordsTrain];
				pVals[i] = new double  [numRecordsTrain];
			}


			for(int k = 0; k < numMeasurementsTrain; k++)
			{
				tTest<double> Tstats;
				Tstats.getPvals(Dvals[k], pVals[k], means[k], vars[k], numRecordsTrain);
			}

			ZScoreMethod zScoreMethod;
			zScoreMethod.OneTpValsToZ(pVals, Zscores, numRecordsTrain, numMeasurementsTrain);

			for(int k = 0; k < numMeasurementsTrain; k++)
			{
				zMean[k] = mean.getMean(Zscores[k], numRecordsTrain);
			}

			for(int k = 0; k < numMeasurementsTrain; k++)
			{
				for(int h = k+1; h < numMeasurementsTrain; h++)
				{
					Correlation<double> correlation;
					cors[k][h] = correlation.getCor(Zscores[k], Zscores[h], numRecordsTrain, zMean[k], zMean[h]);
				}
			}


			for(int i = 0; i < numRecordsClass; i++)
			{
				double crossZvals[numMeasurementsTrain]; double weights[numMeasurementsTrain];
				for(int k = 0; k < numMeasurementsTrain; k++)
				{
					crossDvals[k][i] = leftClass.getRecordMeas(i,k) - rightClass.getRecordMeas(i,k);

					tTest<double> Tstats;

					double tStat = Tstats.getTStatistic(crossDvals[k], i, means[k], vars[k]);
					crossPvals[k][i] = Tstats.twoTailedPVal(tStat, numRecordsTrain);

					crossZvals[k] = zScoreMethod.getZVal1T(crossPvals[k][i]); //zScoreMethod.getZVal2T(crossPvals[k][j]);
					weights[k] =  sqrt(vars[k]);

				}

				double finalP = zScoreMethod.combineZ(crossZvals, weights, cors, numMeasurementsTrain);
				cout<<leftClass.getRecordId(i)<<","<<rightClass.getRecordId(i)<<","<<finalP<<endl;
			}
		}


		if(sumtTest)
		{
			double * DvalSums = new double [numRecordsTrain];
			for(int i = 0; i < numRecordsTrain; i++)
			{
				DvalSums[i] = 0;
				for(int j = 0; j < numMeasurementsTrain; j++)
				{
					DvalSums[i]+=Dvals[j][i];
				}

			}
			double * crossDvalSums = new double [numRecordsTrain];


			double SumMean = mean.getMean(DvalSums, numRecordsTrain);
			double SumVar = var.getVar(DvalSums, numRecordsTrain, SumMean);
			tTest<double> Tstats;

			for(int i = 0; i < numRecordsClass; i++)
			{
				crossDvalSums[i] = 0;
				for(int k = 0; k < numMeasurementsClass; k++)
				{
					crossDvalSums[i]+=rightClass.getRecordMeas(i,k)-leftClass.getRecordMeas(i,k);
				}

				double tStat = Tstats.getTStatistic(crossDvalSums, i, 0, SumVar);
				double pVal = Tstats.twoTailedPVal(tStat, numRecordsTrain);

				cout<<leftClass.getRecordId(i)<<","<<rightClass.getRecordId(i)<<","<<pVal<<endl;
			}
		}


		if(meanSumtTest)
		{
			double * DvalSums = new double [numRecordsTrain];
			for(int i = 0; i < numRecordsTrain; i++)
			{
				DvalSums[i] = 0;
				for(int j = 0; j < numMeasurementsTrain; j++)
				{
					DvalSums[i]+=Dvals[j][i];
				}

			}
			double * crossDvalSums = new double [numRecordsTrain];


			double SumMean = mean.getMean(DvalSums, numRecordsTrain);
			double SumVar = var.getVar(DvalSums, numRecordsTrain, SumMean);
			tTest<double> Tstats;

			for(int i = 0; i < numRecordsClass; i++)
			{
				crossDvalSums[i] = 0;
				for(int k = 0; k < numMeasurementsClass; k++)
				{
					crossDvalSums[i]+=rightClass.getRecordMeas(i,k)-leftClass.getRecordMeas(i,k);
				}

				double tStat = Tstats.getTStatistic(crossDvalSums, i, SumMean, SumVar);
				double pVal = Tstats.twoTailedPVal(tStat, numRecordsTrain);

				cout<<leftClass.getRecordId(i)<<","<<rightClass.getRecordId(i)<<","<<pVal<<endl;
			}
		}


		if(abstTest)
		{
			double * DvalSums = new double [numRecordsTrain];
			for(int i = 0; i < numRecordsTrain; i++)
			{
				DvalSums[i] = 0;
				for(int j = 0; j < numMeasurementsTrain; j++)
				{
					DvalSums[i]+=abs(Dvals[j][i]);
				}

				DvalSums[i] = pow(DvalSums[i] + 0.00005, 0.33);
			}
			double * crossDvalSums = new double [numRecordsTrain];

			double SumMean = mean.getMean(DvalSums, numRecordsTrain);
			double SumVar = var.getVar(DvalSums, numRecordsTrain, SumMean);
			tTest<double> Tstats;

			for(int i = 0; i < numRecordsClass; i++)
			{
				crossDvalSums[i] = 0;
				for(int k = 0; k < numMeasurementsTrain; k++)
				{
					crossDvalSums[i]+=abs(rightTrain.getRecordMeas(i,k)-leftTrain.getRecordMeas(i,k));
				}

				crossDvalSums[i] = pow(crossDvalSums[i] + 0.00005, 0.33);

				double tStat = Tstats.getTStatistic(crossDvalSums, i, SumMean, SumVar);
				double pVal = Tstats.twoTailedPVal(tStat, numRecordsTrain);

				pVal=pVal/2;
				cout<<leftClass.getRecordId(i)<<","<<rightClass.getRecordId(i)<<","<<pVal<<endl;
			}

		}	
	}

	if(LOOCV && unweightedZscore && minSig == 1)
	{
		Summator<double> summator;
		double sums[numMeasurementsTrain];
		double sqSums[numMeasurementsTrain];

		double adjMean[numMeasurementsTrain];
		double adjVars[numMeasurementsTrain];

		double zMean[numMeasurementsTrain];

		for(int i = 0; i < numMeasurementsTrain; i++)
		{
			sums[i] = summator.sum(Dvals[i], numRecordsTrain);
			sqSums[i] = summator.sumXY(Dvals[i], Dvals[i], numRecordsTrain);
		}

		double ** crossDvals = new double * [numMeasurementsTrain];
		double ** crossPvals = new double * [numMeasurementsTrain];
		double * result = new double [numRecordsTrain];
		double ** Zscores = new double * [numMeasurementsTrain];
		double ** pVals = new double * [numMeasurementsTrain];

		for(int i = 0; i < numMeasurementsTrain; i++)
		{
			crossDvals[i] = new double [numRecordsTrain];
			crossPvals[i] = new double [numRecordsTrain];
			Zscores[i] = new double [numRecordsTrain];
			pVals[i] = new double  [numRecordsTrain];
		}


		int LO[2] = {0};
		for(int i = 0; i < numRecordsTrain; i++)
		{
			for(int j = 0; j < numRecordsTrain; j++)
			{
				LO[0] = i; LO[1] = j; int numLO = (i==j) ? 1:2;
				for(int k = 0; k < numMeasurementsTrain; k++)
				{
					adjMean[k] = mean.getMeanLOO(Dvals[k], numRecordsTrain, LO, numLO, sums[k]);
					adjVars[k] = var.getVarLOO(Dvals[k], numRecordsTrain, LO, numLO, sums[k], sqSums[k], adjMean[k]);

					tTest<double> Tstats;

					Tstats.getPvals(Dvals[k], pVals[k], adjMean[k], adjVars[k], numRecordsTrain);
				}

				ZScoreMethod zScoreMethod;
				zScoreMethod.OneTpValsToZ(pVals, Zscores, numRecordsTrain, numMeasurementsTrain);

				for(int k = 0; k < numMeasurementsTrain; k++)
				{
					double sum = summator.sum(Zscores[k], numRecordsTrain);
					zMean[k] = mean.getMeanLOO(Zscores[k], numRecordsTrain, LO, numLO, sum);
				}

				for(int k = 0; k < numMeasurementsTrain; k++)
				{
					for(int h = k+1; h < numMeasurementsTrain; h++)
					{
						Correlation<double> correlation;
						cors[k][h] = correlation.getCorLO(Zscores[k], Zscores[h], numRecordsTrain, LO, numLO, zMean[k], zMean[h]);
					}
				}

				double crossZvals[numMeasurementsTrain]; double weights[numMeasurementsTrain];
				for(int k = 0; k < numMeasurementsTrain; k++)
				{
					crossDvals[k][j] = leftTrain.getRecordMeas(j,k) - rightTrain.getRecordMeas(i,k);

					tTest<double> Tstats;

					double tStat = Tstats.getTStatistic(crossDvals[k], j, adjMean[k], adjVars[k]);
					crossPvals[k][j] = Tstats.twoTailedPVal(tStat, numRecordsTrain-numLO);

					crossZvals[k] = zScoreMethod.getZVal1T(crossPvals[k][j]); //zScoreMethod.getZVal2T(crossPvals[k][j]);
					weights[k] = 1; //abs(tStat); //sqrt(adjVars[k]);

				}

				double finalP = zScoreMethod.combineZ(crossZvals, weights, cors, numMeasurementsTrain); //zScoreMethod.combineZ2T(crossZvals, weights, cors, numMeasurements);
				cout<<leftTrain.getRecordId(j)<<","<<rightTrain.getRecordId(i)<<","<<finalP<<endl;
			}
		}

		for(int i = 0; i < numMeasurementsTrain; i++)
		{
			delete [] crossDvals[i];
			delete [] crossPvals[i];
			delete [] Zscores [i];
			delete [] pVals [i];
		}

		delete [] crossDvals;
		delete [] crossPvals;
		delete [] result;
		delete [] Zscores;
		delete [] pVals;
	}


	if(LOOCV && weightedEffectZscore && minSig == 1)
	{
		Summator<double> summator;
		double sums[numMeasurementsTrain];
		double sqSums[numMeasurementsTrain];

		double adjMean[numMeasurementsTrain];
		double adjVars[numMeasurementsTrain];

		double zMean[numMeasurementsTrain];

		for(int i = 0; i < numMeasurementsTrain; i++)
		{
			sums[i] = summator.sum(Dvals[i], numRecordsTrain);
			sqSums[i] = summator.sumXY(Dvals[i], Dvals[i], numRecordsTrain);
		}

		double ** crossDvals = new double * [numMeasurementsTrain];
		double ** crossPvals = new double * [numMeasurementsTrain];
		double * result = new double [numRecordsTrain];
		double ** Zscores = new double * [numMeasurementsTrain];
		double ** pVals = new double * [numMeasurementsTrain];

		for(int i = 0; i < numMeasurementsTrain; i++)
		{
			crossDvals[i] = new double [numRecordsTrain];
			crossPvals[i] = new double [numRecordsTrain];
			Zscores[i] = new double [numRecordsTrain];
			pVals[i] = new double  [numRecordsTrain];
		}


		int LO[2] = {0};
		for(int i = 0; i < numRecordsTrain; i++)
		{
			for(int j = 0; j < numRecordsTrain; j++)
			{
				LO[0] = i; LO[1] = j; int numLO = (i==j) ? 1:2;
				for(int k = 0; k < numMeasurementsTrain; k++)
				{
					adjMean[k] = mean.getMeanLOO(Dvals[k], numRecordsTrain, LO, numLO, sums[k]);
					adjVars[k] = var.getVarLOO(Dvals[k], numRecordsTrain, LO, numLO, sums[k], sqSums[k], adjMean[k]);

					tTest<double> Tstats;

					Tstats.getPvals(Dvals[k], pVals[k], adjMean[k], adjVars[k], numRecordsTrain);
				}

				ZScoreMethod zScoreMethod;
				zScoreMethod.OneTpValsToZ(pVals, Zscores, numRecordsTrain, numMeasurementsTrain);

				for(int k = 0; k < numMeasurementsTrain; k++)
				{
					double sum = summator.sum(Zscores[k], numRecordsTrain);
					zMean[k] = mean.getMeanLOO(Zscores[k], numRecordsTrain, LO, numLO, sum);
				}

				for(int k = 0; k < numMeasurementsTrain; k++)
				{
					for(int h = k+1; h < numMeasurementsTrain; h++)
					{
						Correlation<double> correlation;
						cors[k][h] = correlation.getCorLO(Zscores[k], Zscores[h], numRecordsTrain, LO, numLO, zMean[k], zMean[h]);
					}
				}

				double crossZvals[numMeasurementsTrain]; double weights[numMeasurementsTrain];
				for(int k = 0; k < numMeasurementsTrain; k++)
				{
					crossDvals[k][j] = leftTrain.getRecordMeas(j,k) - rightTrain.getRecordMeas(i,k);

					tTest<double> Tstats;

					double tStat = Tstats.getTStatistic(crossDvals[k], j, adjMean[k], adjVars[k]);
					crossPvals[k][j] = Tstats.twoTailedPVal(tStat, numRecordsTrain-numLO);

					crossZvals[k] = zScoreMethod.getZVal1T(crossPvals[k][j]); //zScoreMethod.getZVal2T(crossPvals[k][j]);
					weights[k] = abs(tStat); //sqrt(adjVars[k]); //1; 

				}

				double finalP = zScoreMethod.combineZ(crossZvals, weights, cors, numMeasurementsTrain); //zScoreMethod.combineZ2T(crossZvals, weights, cors, numMeasurements);
				cout<<leftTrain.getRecordId(j)<<","<<rightTrain.getRecordId(i)<<","<<finalP<<endl;
			}
		}

		for(int i = 0; i < numMeasurementsTrain; i++)
		{
			delete [] crossDvals[i];
			delete [] crossPvals[i];
			delete [] Zscores [i];
			delete [] pVals [i];
		}

		delete [] crossDvals;
		delete [] crossPvals;
		delete [] result;
		delete [] Zscores;
		delete [] pVals;
	}



	if(LOOCV && weightedStdErrZscore && minSig == 1)
	{
		Summator<double> summator;
		double sums[numMeasurementsTrain];
		double sqSums[numMeasurementsTrain];

		double adjMean[numMeasurementsTrain];
		double adjVars[numMeasurementsTrain];

		double zMean[numMeasurementsTrain];

		for(int i = 0; i < numMeasurementsTrain; i++)
		{
			sums[i] = summator.sum(Dvals[i], numRecordsTrain);
			sqSums[i] = summator.sumXY(Dvals[i], Dvals[i], numRecordsTrain);
		}

		double ** crossDvals = new double * [numMeasurementsTrain];
		double ** crossPvals = new double * [numMeasurementsTrain];
		double * result = new double [numRecordsTrain];
		double ** Zscores = new double * [numMeasurementsTrain];
		double ** pVals = new double * [numMeasurementsTrain];

		for(int i = 0; i < numMeasurementsTrain; i++)
		{
			crossDvals[i] = new double [numRecordsTrain];
			crossPvals[i] = new double [numRecordsTrain];
			Zscores[i] = new double [numRecordsTrain];
			pVals[i] = new double  [numRecordsTrain];
		}


		int LO[2] = {0};
		for(int i = 0; i < numRecordsTrain; i++)
		{
			for(int j = 0; j < numRecordsTrain; j++)
			{
				LO[0] = i; LO[1] = j; int numLO = (i==j) ? 1:2;
				for(int k = 0; k < numMeasurementsTrain; k++)
				{
					adjMean[k] = mean.getMeanLOO(Dvals[k], numRecordsTrain, LO, numLO, sums[k]);
					adjVars[k] = var.getVarLOO(Dvals[k], numRecordsTrain, LO, numLO, sums[k], sqSums[k], adjMean[k]);

					tTest<double> Tstats;

					Tstats.getPvals(Dvals[k], pVals[k], adjMean[k], adjVars[k], numRecordsTrain);
				}

				ZScoreMethod zScoreMethod;
				zScoreMethod.OneTpValsToZ(pVals, Zscores, numRecordsTrain, numMeasurementsTrain);

				for(int k = 0; k < numMeasurementsTrain; k++)
				{
					double sum = summator.sum(Zscores[k], numRecordsTrain);
					zMean[k] = mean.getMeanLOO(Zscores[k], numRecordsTrain, LO, numLO, sum);
				}

				for(int k = 0; k < numMeasurementsTrain; k++)
				{
					for(int h = k+1; h < numMeasurementsTrain; h++)
					{
						Correlation<double> correlation;
						cors[k][h] = correlation.getCorLO(Zscores[k], Zscores[h], numRecordsTrain, LO, numLO, zMean[k], zMean[h]);
					}
				}

				double crossZvals[numMeasurementsTrain]; double weights[numMeasurementsTrain];
				for(int k = 0; k < numMeasurementsTrain; k++)
				{
					crossDvals[k][j] = leftTrain.getRecordMeas(j,k) - rightTrain.getRecordMeas(i,k);

					tTest<double> Tstats;

					double tStat = Tstats.getTStatistic(crossDvals[k], j, adjMean[k], adjVars[k]);
					crossPvals[k][j] = Tstats.twoTailedPVal(tStat, numRecordsTrain-numLO);

					crossZvals[k] = zScoreMethod.getZVal1T(crossPvals[k][j]); //zScoreMethod.getZVal2T(crossPvals[k][j]);
					weights[k] = sqrt(adjVars[k]); //1; //abs(tStat); 

				}

				double finalP = zScoreMethod.combineZ(crossZvals, weights, cors, numMeasurementsTrain); //zScoreMethod.combineZ2T(crossZvals, weights, cors, numMeasurements);
				cout<<leftTrain.getRecordId(j)<<","<<rightTrain.getRecordId(i)<<","<<finalP<<endl;
			}
		}

		for(int i = 0; i < numMeasurementsTrain; i++)
		{
			delete [] crossDvals[i];
			delete [] crossPvals[i];
			delete [] Zscores [i];
			delete [] pVals [i];
		}

		delete [] crossDvals;
		delete [] crossPvals;
		delete [] result;
		delete [] Zscores;
		delete [] pVals;
	}



	if(LOOCV && minSig > 1) //at least U significant
	{
		Summator<double> summator;
		double sums[numMeasurementsTrain];
		double sqSums[numMeasurementsTrain];

		double adjMean[numMeasurementsTrain];
		double adjVars[numMeasurementsTrain];

		double zMean[numMeasurementsTrain];

		int  orderZ  [numMeasurementsTrain];

		for(int i = 0; i < numMeasurementsTrain; i++)
		{
			sums[i] = summator.sum(Dvals[i], numRecordsTrain);
			sqSums[i] = summator.sumXY(Dvals[i], Dvals[i], numRecordsTrain);
		}

		double ** crossDvals = new double * [numMeasurementsTrain];
		double ** crossPvals = new double * [numMeasurementsTrain];
		double * result = new double [numRecordsTrain];
		double ** Zscores = new double * [numMeasurementsTrain];
		double ** pVals = new double * [numMeasurementsTrain];

		for(int i = 0; i < numMeasurementsTrain; i++)
		{
			crossDvals[i] = new double [numRecordsTrain];
			crossPvals[i] = new double [numRecordsTrain];
			Zscores[i] = new double [numRecordsTrain];
			pVals[i] = new double  [numRecordsTrain];
		}

		int LO[2] = {0};
		for(int i = 0; i < numRecordsTrain; i++)
		{
			for(int j = 0; j < numRecordsTrain; j++)
			{
				LO[0] = i; LO[1] = j; int numLO = (i==j) ? 1:2;
				for(int k = 0; k < numMeasurementsTrain; k++)
				{
					adjMean[k] = mean.getMeanLOO(Dvals[k], numRecordsTrain, LO, numLO, sums[k]);
					adjVars[k] = var.getVarLOO(Dvals[k], numRecordsTrain, LO, numLO, sums[k], sqSums[k], adjMean[k]);
					tTest<double> Tstats;

					Tstats.getPvals(Dvals[k], pVals[k], adjMean[k], adjVars[k], numRecordsTrain);
					orderZ[k] = k;
				}

				ZScoreMethod zScoreMethod;
				zScoreMethod.TwoTpValsToZ(pVals, Zscores, numRecordsTrain, numMeasurementsTrain);
				for(int k = 0; k < numMeasurementsTrain; k++)
				{
					double sum = summator.sum(Zscores[k], numRecordsTrain);
					zMean[k] = mean.getMeanLOO(Zscores[k], numRecordsTrain, LO, numLO, sum);
				}

				for(int k = 0; k < numMeasurementsTrain; k++)
				{
					for(int h = 0; h < numMeasurementsTrain; h++)
					{
						Correlation<double> correlation;
						cors[k][h] = correlation.getCorLO(Zscores[k], Zscores[h], numRecordsTrain, LO, numLO, zMean[k], zMean[h]);
					}
				}

				double crossZvals[numMeasurementsTrain]; double weights[numMeasurementsTrain];
				for(int k = 0; k < numMeasurementsTrain; k++)
				{
					crossDvals[k][j] = rightTrain.getRecordMeas(i,k)-leftTrain.getRecordMeas(j,k);
					tTest<double> Tstats;

					double tStat = Tstats.getTStatistic(crossDvals[k], j, adjMean[k], adjVars[k]);
					crossPvals[k][j] = Tstats.twoTailedPVal(tStat, numRecordsTrain-numLO); //CHANGED

					crossZvals[k] = zScoreMethod.getZVal2T(crossPvals[k][j]);
					weights[k] = 1/sqrt(adjVars[k]);
				}

				int numTest = numMeasurementsTrain - minSig + 1;

				double finalP = zScoreMethod.combineZ2TOrdered(crossZvals, weights, cors, orderZ, numTest, numMeasurementsTrain);
				cout<<leftTrain.getRecordId(j)<<","<<rightTrain.getRecordId(i)<<","<<finalP<<endl;
			}
		}

		for(int i = 0; i < numMeasurementsTrain; i++)
		{
			delete [] crossDvals[i];
			delete [] crossPvals[i];
			delete [] Zscores [i];
			delete [] pVals [i];
		}

		delete [] crossDvals;
		delete [] crossPvals;
		delete [] result;
		delete [] Zscores;
		delete [] pVals;
	}


	if(LOOCV && sumtTest)
	{
		double * DvalSums = new double [numRecordsTrain];
		for(int i = 0; i < numRecordsTrain; i++)
		{
			DvalSums[i] = 0;
			for(int j = 0; j < numMeasurementsTrain; j++)
			{
				DvalSums[i]+=Dvals[j][i];
			}

		}
		double * crossDvalSums = new double [numRecordsTrain];

		Summator<double> summator;
		double sum = summator.sum(DvalSums, numRecordsTrain);
		double sqSum = summator.sumXY(DvalSums, DvalSums, numRecordsTrain);
		tTest<double> Tstats;

		int LO[2] = {0};
		for(int i = 0; i < numRecordsTrain; i++)
		{
			for(int j = 0; j < numRecordsTrain; j++)
			{
				LO[0] = i; LO[1] = j; int numLO = (i==j) ? 1:2;
				double adjMean = mean.getMeanLOO(DvalSums, numRecordsTrain, LO, numLO, sum);
				double adjVar = var.getVarLOO(DvalSums, numRecordsTrain, LO, numLO, sum, sqSum, adjMean);

				crossDvalSums[j] = 0;
				for(int k = 0; k < numMeasurementsTrain; k++)
				{
					crossDvalSums[j]+=(rightTrain.getRecordMeas(i,k)-leftTrain.getRecordMeas(j,k));
				}


				double tStat = Tstats.getTStatistic(crossDvalSums, j, 0, adjVar);
				double pVal = Tstats.twoTailedPVal(tStat, numRecordsTrain-numLO);

				cout<<leftTrain.getRecordId(j)<<","<<rightTrain.getRecordId(i)<<","<<pVal<<endl;

			}
		}

		delete [] crossDvalSums;
	}


	if(LOOCV && meanSumtTest)
	{
		double * DvalSums = new double [numRecordsTrain];
		for(int i = 0; i < numRecordsTrain; i++)
		{
			DvalSums[i] = 0;
			for(int j = 0; j < numMeasurementsTrain; j++)
			{
				DvalSums[i]+=Dvals[j][i];
			}

		}
		double * crossDvalSums = new double [numRecordsTrain];

		Summator<double> summator;
		double sum = summator.sum(DvalSums, numRecordsTrain);
		double sqSum = summator.sumXY(DvalSums, DvalSums, numRecordsTrain);
		tTest<double> Tstats;

		int LO[2] = {0};
		for(int i = 0; i < numRecordsTrain; i++)
		{
			for(int j = 0; j < numRecordsTrain; j++)
			{
				LO[0] = i; LO[1] = j; int numLO = (i==j) ? 1:2;
				double adjMean = mean.getMeanLOO(DvalSums, numRecordsTrain, LO, numLO, sum);
				double adjVar = var.getVarLOO(DvalSums, numRecordsTrain, LO, numLO, sum, sqSum, adjMean);

				crossDvalSums[j] = 0;
				for(int k = 0; k < numMeasurementsTrain; k++)
				{
					crossDvalSums[j]+=(rightTrain.getRecordMeas(i,k)-leftTrain.getRecordMeas(j,k));
				}

				double tStat = Tstats.getTStatistic(crossDvalSums, j, adjMean, adjVar);
				double pVal = Tstats.twoTailedPVal(tStat, numRecordsTrain-numLO);

				cout<<leftTrain.getRecordId(j)<<","<<rightTrain.getRecordId(i)<<","<<pVal<<endl;

			}
		}

		delete [] crossDvalSums;
	}


	if(LOOCV && abstTest)
	{
		double * DvalSums = new double [numRecordsTrain];
		for(int i = 0; i < numRecordsTrain; i++)
		{
			DvalSums[i] = 0;
			for(int j = 0; j < numMeasurementsTrain; j++)
			{
				DvalSums[i]+=abs(Dvals[j][i]);
			}

			DvalSums[i] = pow(DvalSums[i] + 0.00005, 0.33);
		}
		double * crossDvalSums = new double [numRecordsTrain];

		Summator<double> summator;
		double sum = summator.sum(DvalSums, numRecordsTrain);
		double sqSum = summator.sumXY(DvalSums, DvalSums, numRecordsTrain);
		tTest<double> Tstats;

		int LO[2] = {0};
		for(int i = 0; i < numRecordsTrain; i++)
		{
			for(int j = 0; j < numRecordsTrain; j++)
			{
				LO[0] = i; LO[1] = j; int numLO = (i==j) ? 1:2;
				double adjMean = mean.getMeanLOO(DvalSums, numRecordsTrain, LO, numLO, sum);
				double adjVar = var.getVarLOO(DvalSums, numRecordsTrain, LO, numLO, sum, sqSum, adjMean);

				crossDvalSums[j] = 0;
				for(int k = 0; k < numMeasurementsTrain; k++)
				{
					crossDvalSums[j]+=abs(rightTrain.getRecordMeas(i,k)-leftTrain.getRecordMeas(j,k));
				}

				crossDvalSums[j] = pow(crossDvalSums[j] + 0.00005, 0.33);

				double tStat = Tstats.getTStatistic(crossDvalSums, j, adjMean, adjVar);
				double pVal = Tstats.twoTailedPVal(tStat, numRecordsTrain-numLO);

				pVal = pVal/2;

				cout<<leftTrain.getRecordId(j)<<","<<rightTrain.getRecordId(i)<<","<<pVal<<endl;
			}
		}

		delete [] crossDvalSums;
	}



	for(int i = 0; i < numMeasurementsTrain; i++){
		delete [] Dvals [i];
		delete [] cors [i];
	}


	delete [] Dvals;
	delete [] cors;

	return 0;
}

