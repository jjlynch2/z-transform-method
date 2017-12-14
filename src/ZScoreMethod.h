/*
 * ZScoreMethod.h
 *
 *  Created on: Jan 16, 2017
 *      Author: Julia
 */

#ifndef ZSCOREMETHOD_H_
#define ZSCOREMETHOD_H_

class ZScoreMethod {
public:
	ZScoreMethod();
	virtual ~ZScoreMethod();
	void TwoTpValsToZ(double **,  double **, const int, const int);
	void OneTpValsToZ(double **, double **, const int, const int);
	void TwoTpValsToOneToZ(double **, double **, const int, const int);
	void combineZs(double **, double *, double *, double **, const int , const int);
	double combineZ(double [], double *, double **, const int);
	void combineZs2T(double **, double *, double *, double **, const int , const int);
	double combineZ2T(double [], double [], double **,  const int);
	double combineZ2TOrdered(double [], double [], double **, int [], const int, const int);
	double getZVal2T(const double);
	double getZVal2TP(const double);
	double getZVal1T(const double);
private:
	double getDenominator(double *, double **, const int);
	double getDenominatorOrdered(double *, double **, int [], const int, const int);
	double getNumerator(double **, const int, double *, const int);
	double getNumerator(double [],  double *, const int);
	double getNumeratorOrdered(double [], double *, int [], const int, const int);
};

#endif /* ZSCOREMETHOD_H_ */
