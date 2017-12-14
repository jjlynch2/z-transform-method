/*
 * LinearRegression.cpp
 *
 *  Created on: Feb 3, 2017
 *      Author: Julia
 */

#include <iostream>
#include<math.h>
using namespace std;
#include "Summator.h"
#include "LinearRegression.h"

LinearRegression::LinearRegression(double * X, double * Y, const int arrSize) {
			Summator<double>sum;

			XXsum = sum.sumXY(X, X, arrSize);
			Xsum = sum.sum(X, arrSize);
			Ysum = sum.sum(Y, arrSize);
			XYsum = sum.sumXY(X,Y, arrSize);
}

LinearRegression::~LinearRegression() {
	// TODO Auto-generated destructor stub
}

double LinearRegression::getM(const int arrSize) {

		double numerator = (arrSize * XYsum)-(Xsum*Ysum);
		double denominator = (arrSize * XXsum) - (Xsum * Xsum);

		return numerator/denominator;
}

double LinearRegression::getMLO(double * X, double * Y, const int arrSize, const int LO [], const int numLO) {

	double XXsumLO, XsumLO, YsumLO, XYsumLO = 0;

	for(int i = 0; i < numLO; i++)
	{
		XXsumLO = XXsum - (X[LO[i]] * X[LO[i]]);
		XsumLO = Xsum - X[LO[i]];
		YsumLO = Ysum - Y[LO[i]];
		XYsumLO = XYsum - X[LO[i]]*Y[LO[i]];
	}

	double numerator = ((arrSize-numLO) * XYsumLO)-(XsumLO*YsumLO);
	double denominator = ((arrSize-numLO) * XXsumLO) - (XsumLO * XsumLO);

	return numerator/denominator;
}

double LinearRegression::getB(const int arrSize) {

		double numerator = (XXsum*Ysum)-(Xsum*XYsum);
		double denominator = (arrSize * XXsum) - (Xsum * Xsum);

		return numerator/denominator;
}

double LinearRegression::getBLO(double * X, double * Y, const int arrSize, const int LO [], const int numLO) {

	double XXsumLO, XsumLO, YsumLO, XYsumLO = 0;

	for(int i = 0; i < numLO; i++)
	{
			XXsumLO = XXsum - (X[LO[i]] * X[LO[i]]);
			XsumLO = Xsum - X[LO[i]];
			YsumLO = Ysum - Y[LO[i]];
			XYsumLO = XYsum - X[LO[i]]*Y[LO[i]];
	}

	double numerator = (XXsumLO*YsumLO)-(XsumLO*XYsumLO);
	double denominator = ((arrSize-numLO) * XXsumLO) - (XsumLO * XsumLO);

	return numerator/denominator;
}

double LinearRegression::getSE(double * X, double * Y, const int arrSize) {

	double M = getM(arrSize); double B = getB(arrSize);

	double sum = 0;
	for(int i = 0; i < arrSize; i++)
	{
		double predY = M * X[i] + B;
		sum += pow(Y[i] - predY, 2);
	}

	double N = arrSize;
	return sqrt(sum/N);
}

double LinearRegression::getSELO(double * X, double * Y, const double SE, const int arrSize, const int LO [], const int numLO) {

	double M = getMLO(X, Y, arrSize, LO, numLO); double B = getBLO(X, Y, arrSize, LO, numLO);

	double sum = pow(SE,2) * arrSize;

	for(int i = 0; i < numLO; i++)
	{
		double predY = M * X[LO[i]] + B;
		sum-= pow(Y[LO[i]] - predY, 2);
	}

	double N = arrSize-numLO;
	return sqrt(sum/N);
}


double LinearRegression::getTStatistic(const double * X, const double * Y, const int & pos1, const int & pos2, const double & mean, const double & var, const double & SE, const double & arrSize) {

	double M = getM(arrSize); double B = getB(arrSize);
	double predY = M * X[pos1] + B;

	double sum1 = pow((X[pos1] - mean),2)/(arrSize * var);

	double numerator = predY - Y[pos2];
	double denominator = SE * sqrt(1 + 1/arrSize + sum1);

	return numerator/denominator;
}
