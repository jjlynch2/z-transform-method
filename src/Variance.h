/*
 * Variance.h
 *
 *  Created on: Jan 14, 2017
 *      Author: Julia
 */

#ifndef VARIANCE_H_
#define VARIANCE_H_

#include <iostream>
#include <math.h>
using namespace std;

template <class T>
class Variance {
public:
	Variance();
	virtual ~Variance();
	double getVar(const T *, const int, const double) const;
	double getVarLOO(const T *, const int, const int [], const int, const double, const double, const double);
};

template <class T>
Variance<T>::Variance() {};

template <class T>
Variance<T>::~Variance() {
}

template <class T>
double Variance<T>::getVar(const T * arr, const int arrSize, const double mean) const {

		double sum = 0;
		for(int i = 0; i < arrSize; i++)
		{
			sum+=(pow(arr[i] - mean, 2));
		}

		float size = arrSize;
		return sum/(size-1);
}

template <class T>
double Variance<T>::getVarLOO(const T * arr, const int arrSize, const int LO [], const int numLO, const double sum, const double sqSum, const double mean) {

	double firstTerm = sqSum;
	for(int i = 0; i < numLO; i++)
			firstTerm-=pow(arr[LO[i]],2);

	double squaredMean = pow(mean,2);
	double lastTerm = (arrSize-numLO) * squaredMean;

	double sumLOO = sum;
	for(int i = 0; i < numLO; i++)
		sumLOO-=arr[LO[i]];
	double middleTerm = -2*(mean * sumLOO);

	double denominator = arrSize-1-numLO;

	double newVar = (firstTerm + middleTerm + lastTerm)/denominator;

	return newVar;
}

#endif /* VARIANCE_H_ */
