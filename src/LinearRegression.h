/*
 * LinearRegression.h
 *
 *  Created on: Feb 3, 2017
 *      Author: Julia
 */

#ifndef LINEARREGRESSION_H_
#define LINEARREGRESSION_H_

class LinearRegression {
public:
	LinearRegression(double *, double *, const int);
	virtual ~LinearRegression();
	double getM(const int);
	double getMLO(double *, double *, const int, const int [], const int);
	double getB(const int);
	double getBLO(double *, double *, const int, const int [], const int);
	double getSE(double *, double *, const int);
	double getSELO(double *, double *, const double, const int, const int [], const int);
	double getTStatistic(const double *, const double *, const int &, const int &, const double &, const double &, const double &, const double &);
private:
	double XXsum;
	double Xsum;
	double Ysum;
	double XYsum;
};

#endif /* LINEARREGRESSION_H_ */
