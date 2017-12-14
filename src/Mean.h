/*
 * Mean.h
 *
 *  Created on: Jan 14, 2017
 *      Author: Julia
 */

#ifndef MEAN_H_
#define MEAN_H_

#include <iostream>
using namespace std;

template <class T>
class Mean {
public:
	Mean();
	virtual ~Mean();
	double getMean(const T *, const int) const;
	double getMeanLOO(const T *, const int, const int [], const int, const double) const;
};

template <class T>
Mean<T>::Mean() {};

template <class T>
Mean<T>::~Mean() {
}

template <class T>
double Mean<T>::getMean(const T * arr, const int arrSize) const {

		double sum = 0;
		for(int i = 0; i < arrSize; i++)
		{
			sum+=arr[i];
		}

		float size = arrSize;
		return sum/size;
}

template <class T>
double Mean<T>::getMeanLOO(const T * arr, const int arrSize, const int LO [], const int numLO, const double sum) const
{
	double numerator = sum;
	for(int i = 0; i < numLO; i++)
	numerator-=arr[LO[i]];

	double denominator = arrSize-numLO;

	return (numerator/denominator);
}

#endif /* MEAN_H_ */
