/*
 * tTest.h
 *
 *  Created on: Jan 15, 2017
 *      Author: Julia
 */

#ifndef TTEST_H_
#define TTEST_H_
#include<iostream>
#include<math.h>
#include<Rcpp.h>
#include<RInside.h>
using namespace std;

template <class T>
class tTest {
public:
	tTest();
	void getPvals(const T *, T *, const double &, const double &, const double &);
	double getTStatistic(const T *, const int &, const double &, const double &, const double &);
	double getTStatistic(const T *, const int &, const double &, const double &);
	double leftTailPVal(const double &, const int &);
	double rightTailPVal(const double &, const int &);
	double twoTailedPVal(const double &, const int &);
	virtual ~tTest();

};

template <class T>
tTest<T>::tTest() {};

template <class T>
tTest<T>::~tTest() {
}

template <class T>
double tTest<T>::getTStatistic(const T * arr, const int & pos, const double & mean, const double & var, const double & sampleSize) {
	return ((arr[pos]-mean)/(sqrt(var)/sampleSize));
}

template <class T>
double tTest<T>::getTStatistic(const T * arr, const int & pos, const double & mean, const double & var) {

	return ((arr[pos]-mean)/(sqrt(var)));
}



template <class T>
void tTest<T>::getPvals(const T * arr, T * pArr, const double & mean, const double & var, const double & sampleSize) {

		for(int i = 0; i < sampleSize; i++)
		{
			double Tstatistic = getTStatistic(arr, i, mean, var);
			//cout<<"TSTAT "<<Tstatistic<<endl;
			pArr[i] =  twoTailedPVal(Tstatistic, sampleSize-1);
			//cout<<pArr[i]<<" PVAL "<<endl;
		}
}


template <class T>
double tTest<T>::leftTailPVal(const double & tStat, const int & df) {

	Rcpp::NumericVector x(1,tStat);
	Rcpp::NumericVector Y = Rcpp::pt(x, df, true, false);

	return Y[0];
}

template <class T>
double tTest<T>::rightTailPVal(const double & tStat, const int & df) {

	Rcpp::NumericVector x(1,tStat);
	Rcpp::NumericVector Y = Rcpp::pt(x, df, false, false);

	return Y[0];
}

template <class T>
double tTest<T>::twoTailedPVal(const double & tStat, const int & df) {
	Rcpp::NumericVector x(1, tStat);
	Rcpp::NumericVector Y = 2 * Rcpp::pt(Rcpp::abs(x), df, false, false);

	return Y[0];
}


#endif /* TTEST_H_ */
