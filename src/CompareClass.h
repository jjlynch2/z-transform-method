/*
 * CompareClass.h
 *
 *  Created on: Feb 24, 2017
 *      Author: Julia
 */

#ifndef COMPARECLASS_H_
#define COMPARECLASS_H_

struct CompareClass {

	public:
	CompareClass(double * & Z): myArr(Z){;;;}
	bool operator()(const int & i1, const int & i2)
	{
		return (myArr[i1] < myArr[i2]);
	}

private:
	double * myArr;
};

#endif /* COMPARECLASS_H_ */
