/*
 * SEinventory.h
 *
 *  Created on: Jan 14, 2017
 *      Author: Julia
 */

#ifndef SEINVENTORY_H_
#define SEINVENTORY_H_

class SEinventory {
public:
	SEinventory();
	SEinventory(int, int, int, string &, string &);
	virtual ~SEinventory();
	void addMeasurements(string &, float * &);
	string getRecordId(int) const;
	float getRecordMeas(int, int) const;
	string getElementType(void) const;
	string getElementSide(void) const;
private:
	char * IDs;
	double ** Measurements;
	int numMeasPerRecord;
	int numRecords;
	int recordCapacity;
	int idLength;
	string elementType;
	string elementSide;
};

#endif /* SEINVENTORY_H_ */
