/*
 * Summator.h
 *
 *  Created on: Jan 15, 2017
 *      Author: Julia
 */

#ifndef SUMMATOR_H_
#define SUMMATOR_H_

template <class T>
class Summator {
public:
	Summator();
	virtual ~Summator();
	double sum(const T *, const int) const;
	double sumLOO(const T *, const int, const int [], const int) const;
	double sumXY(const T *, const T *, const int) const;
	double sumXYLOO(const T *, const T *, const int, const int [], const int) const;
};

template <class T>
Summator<T>::Summator() {};

template <class T>
Summator<T>::~Summator() {
}

template <class T>
double Summator<T>::sum(const T * arr, const int arrSize) const {

	double sum = 0;
	for(int i = 0; i < arrSize; i++)
		sum+=arr[i];

	return sum;
}

template <class T>
double Summator<T>::sumLOO(const T * arr, const int arrSize, const int LOO[], const int numLO) const {

	double sum = 0;
	for(int i = 0; i < arrSize; i++)
		sum+=arr[i];

	for(int i = 0; i < numLO; i++)
		sum-=arr[LOO[i]];
	return sum;
}

template <class T>
double Summator<T>::sumXY(const T * X, const T * Y, const int arrSize) const {

	double sum = 0;
	for(int i = 0; i < arrSize; i++)
		sum+=(X[i]*Y[i]);

	return sum;
}

template <class T>
double Summator<T>::sumXYLOO(const T * X, const T * Y, const int arrSize, const int LOO[], const int numLO) const {

	double sum = 0;
	for(int i = 0; i < arrSize; i++)
		sum+=(X[i]*Y[i]);

	for(int i = 0; i < numLO; i++)
			sum-=(X[LOO[i]] * Y[LOO[i]]);

	return sum;
}

#endif /* SUMMATOR_H_ */
