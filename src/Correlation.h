/*
 * Correlation.h
 *
 *  Created on: Jan 14, 2017
 *      Author: Julia
 */

#ifndef CORRELATION_H_
#define CORRELATION_H_

template <class T>
class Correlation {
public:
	Correlation();
	virtual ~Correlation();
	double getCor(const T *, const T *, const int, const double, const double) const;
	double getCorAdj(const double, const double) const;
	double getCorLO(const T *, const T *, const int,const int [], const int, const double, const double) const;

};

template <class T>
Correlation<T>::Correlation() {};

template <class T>
Correlation<T>::~Correlation() {
}

template <class T>
double Correlation<T>::getCor(const T * X, const T * Y, const int arrSize, const double meanX, const double meanY) const {

	double numerator = 0;  double sum1 = 0; double sum2 = 0;
	for(int i = 0; i < arrSize; i++)
	{
		sum1+=pow((X[i]-meanX), 2); sum2+=pow((Y[i]-meanY), 2);
		numerator+=((X[i]-meanX) * (Y[i]-meanY));
	}

	double denominator = (sqrt(sum1)*sqrt(sum2));
	return numerator/denominator;
}

template <class T>
double Correlation<T>::getCorLO(const T * X, const T * Y, const int arrSize, const int LO [], const int numLO, const double meanX, const double meanY) const {

		double numerator = 0;  double sum1 = 0; double sum2 = 0;
		for(int i = 0; i < arrSize; i++)
		{
			sum1+=pow((X[i]-meanX), 2); sum2+=pow((Y[i]-meanY), 2);
			numerator+=((X[i]-meanX) * (Y[i]-meanY));
		}

		for(int i = 0; i < numLO; i++)
		{
			sum1-=pow((X[i]-meanX), 2);
			sum2-=pow((Y[i]-meanY), 2);
			numerator-=((X[i]-meanX) * (Y[i]-meanY));
		}

		double denominator = (sqrt(sum1)*sqrt(sum2));
		return numerator/denominator;
}


template <class T>
double Correlation<T>::getCorAdj(const double r, const double size) const {
	return r*(1+(1-pow(r,2))/(2*size));
}

#endif /* CORRELATION_H_ */
