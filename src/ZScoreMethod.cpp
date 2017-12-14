/*
 * ZScoreMethod.cpp
 *
 *  Created on: Jan 16, 2017
 *      Author: Julia
 */
#include<iostream>
#include<math.h>
#include<Rcpp.h>
#include<RInside.h>
using namespace std;
#include "CompareClass.h"
#include "ZScoreMethod.h"
#include "Correlation.h"

ZScoreMethod::ZScoreMethod() {
}

ZScoreMethod::~ZScoreMethod() {
	// TODO Auto-generated destructor stub
}

void ZScoreMethod::TwoTpValsToZ(double ** pVals, double ** Z, const int numRecords, const int numMeasPerRecord) {

		for(int i = 0; i < numMeasPerRecord; i++)
		{
			for(int j = 0; j < numRecords; j++)
			{
				Z[i][j] = getZVal2T(pVals[i][j]);
			}
		}
}

void ZScoreMethod::OneTpValsToZ(double ** pVals, double ** Z, const int numRecords, const int numMeasPerRecord) {

		for(int i = 0; i < numMeasPerRecord; i++)
		{
			for(int j = 0; j < numRecords; j++)
			{
				Z[i][j] = getZVal1T(pVals[i][j]);
			}
		}
}

void ZScoreMethod::TwoTpValsToOneToZ(double ** pVals, double ** Z, const int numRecords, const int numMeasPerRecord) {

		for(int i = 0; i < numMeasPerRecord; i++)
		{
			for(int j = 0; j < numRecords; j++)
			{
				Z[i][j] = getZVal2T(pVals[i][j]);
			}
		}
}

void ZScoreMethod::combineZs(double ** Z, double * result, double * weights, double ** cor, const int numRecords, const int numMeasPerRecord){

	double denominator = getDenominator(weights, cor, numMeasPerRecord);

	for(int i = 0; i < numRecords; i++)
	{
		double numerator = getNumerator(Z, i, weights, numMeasPerRecord);

		double val1 = numerator/denominator;
		Rcpp::NumericVector x(1,val1);
		Rcpp::NumericVector Z_final = Rcpp::pnorm(x, 0, 1, TRUE, FALSE); //SECOND TRUE FROM FALSE

		double pZ = 1-Z_final[0];

		result[i] = pZ;
	}
}

double ZScoreMethod::combineZ(double Z [], double weights [], double ** cor, const int numMeasPerRecord){

		double denominator = getDenominator(weights, cor, numMeasPerRecord);
		double numerator = getNumerator(Z, weights, numMeasPerRecord);

		double val1 = numerator/denominator;
		Rcpp::NumericVector x(1,val1);
		Rcpp::NumericVector Z_final = Rcpp::pnorm(x, 0, 1, TRUE, FALSE); //SECOND TRUE FROM FALSE

		double pZ = 1-Z_final[0];

		return pZ;
}

void ZScoreMethod::combineZs2T(double ** Z, double * result, double * weights, double ** cor, const int numRecords, const int numMeasPerRecord){

	double denominator = getDenominator(weights, cor, numMeasPerRecord);
	for(int i = 0; i < numRecords; i++)
	{
		double numerator = getNumerator(Z, i, weights, numMeasPerRecord);

		double val1 = numerator/denominator;
		Rcpp::NumericVector x(1,val1);
		Rcpp::NumericVector Z = Rcpp::pnorm(x, 0, 1, TRUE, FALSE); //SECOND TRUE FROM FALSE

		double pZ = 1-Z[0];

		if(pZ < .5)
		{
			pZ=2*pZ;
		}else{
			pZ=2*(1-pZ);
		}
		result[i] = pZ;
	}
}

double ZScoreMethod::combineZ2T(double Z [], double weights [], double ** cor,  const int numMeasPerRecord){

	double denominator = getDenominator(weights, cor, numMeasPerRecord);
	double numerator = getNumerator(Z,  weights, numMeasPerRecord);

	double val1 = numerator/denominator;

	Rcpp::NumericVector x(1,val1);
	Rcpp::NumericVector Z_final = Rcpp::pnorm(x, 0, 1, TRUE, FALSE); //SECOND TRUE FROM FALSE

	double pZ = 1-Z_final[0];

	if(pZ < .5)
	{
		pZ=2*pZ;
	}else{
		pZ=2*(1-pZ);
	}
	return pZ;
}

double ZScoreMethod::combineZ2TOrdered(double Z [], double weights [], double ** cor, int orderZ [], int U, const int numMeasPerRecord){

	sort(orderZ, orderZ+numMeasPerRecord, CompareClass(Z));

	double denominator = getDenominatorOrdered(weights, cor, orderZ, U, numMeasPerRecord);
	double numerator = getNumeratorOrdered(Z,  weights, orderZ, U, numMeasPerRecord);

	double val1 = numerator/denominator;


	Rcpp::NumericVector x(1,val1);
	Rcpp::NumericVector Z_final = Rcpp::pnorm(x, 0, 1, TRUE, FALSE); //SECOND TRUE FROM FALSE

	double pZ = 1-Z_final[0];

	if(pZ < .5)
	{
		pZ=2*pZ;
	}else{
		pZ=2*(1-pZ);
	}
	return pZ;
}

double ZScoreMethod::getZVal2T(const double pVal) {

	long double val1 = (1-pVal/2);
	if(val1 == 1) val1 = .99999999;

	Rcpp::NumericVector x(1,val1);

	Rcpp::NumericVector Z = Rcpp::qnorm(x,0,1,TRUE,FALSE); //SECOND TRUE FROM FALSE

	return Z[0];
}

double ZScoreMethod::getZVal1T(const double pVal) {

	long double val1 = (1-pVal);
	if(val1 == 1) val1 = .99999999;
	if(val1 == 0) val1 = .00000001;
	Rcpp::NumericVector x(1,val1);

	Rcpp::NumericVector Z = Rcpp::qnorm(x,0,1,TRUE,FALSE); //SECOND TRUE FROM FALSE

	return Z[0];
}

double ZScoreMethod::getDenominator(double * weights, double ** cors, const int numPvals) {

	double sum1 = 0;
 	for(int i = 0; i < numPvals; i++)
	{
		sum1+=(pow(weights[i],2));
	}

 	double sum2 = 0;
 	for(int i = 0; i < numPvals; i++)
 	{
 		for(int j = i+1; j < numPvals; j++)
 		{
 			sum2+=(weights[i]*weights[j]*cors[i][j]);
 		}
 	}


 	double sum3 = sum1 + 2*sum2;
 	return sqrt(sum3);
}

double ZScoreMethod::getDenominatorOrdered(double * weights, double ** cors, int order [], int U, const int numPvals) {

	double sum1 = 0;
 	for(int i = 0; i < U; i++)
	{
		sum1+=(pow(weights[order[i]],2));
	}

 	double sum2 = 0;
 	for(int i = 0; i < U; i++)
 	{
 		for(int j = i+1; j < U; j++)
 		{
 			sum2+=(weights[order[i]]*weights[order[j]]*cors[order[i]][order[j]]);
 		}
 	}

 	double sum3 = sum1 + 2*sum2;
 	return sqrt(sum3);
}

double ZScoreMethod::getNumerator(double ** zScores, const int currRecord, double * weights, const int numPvals) {

		double sum = 0;
		for(int i = 0; i < numPvals; i++)
		{
			sum+=(weights[i]*zScores[i][currRecord]);
		}
		return sum;
}

double ZScoreMethod::getNumerator(double zScores [],  double * weights, const int numPvals) {

		double sum = 0;
		for(int i = 0; i < numPvals; i++)
		{
			sum+=(weights[i]*zScores[i]);
		}
		return sum;
}

double ZScoreMethod::getNumeratorOrdered(double zScores [],  double * weights, int order [], const int U, const int numPvals) {

		double sum = 0;
		for(int i = 0; i < U; i++)
		{
			sum+=(weights[order[i]]*zScores[order[i]]);
		}
		return sum;
}
