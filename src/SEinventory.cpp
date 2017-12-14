/*
 * SEinventory.cpp
 *
 *  Created on: Jan 14, 2017
 *      Author: Julia
 */


#include <iostream>
#include <string>
#include <cstring>
using namespace std;
#include "SEinventory.h"



SEinventory::SEinventory(): IDs(NULL), Measurements(NULL), numMeasPerRecord(0), numRecords(0), recordCapacity(0), idLength(0), elementType(NULL), elementSide(NULL) {
	// TODO Auto-generated constructor stub

}

SEinventory::SEinventory(int nm, int rc, int ilen, string & et, string & es ):  numMeasPerRecord(nm),  numRecords(0), recordCapacity(rc), idLength(ilen),
		elementType(et), elementSide(es){

		Measurements = new double * [numMeasPerRecord];

		for(int i = 0; i < numMeasPerRecord; i++)
		{
			Measurements[i] = new double [recordCapacity];
		}

		IDs = new char [recordCapacity * (idLength+1)];
}

SEinventory::~SEinventory() {

	for(int i = 0; i < numMeasPerRecord; i++)
		delete [] Measurements[i];

	delete [] Measurements;
	delete [] IDs;
}

void SEinventory::addMeasurements(string & id, float * & measures) {

	strncpy(IDs+(numRecords*(idLength+1)), id.c_str(), idLength);

	for(int i = 0; i < numMeasPerRecord; i++)
		Measurements[i][numRecords] = measures[i];

	numRecords++;
}

string SEinventory::getRecordId(int record) const{

	int rPos = record * (idLength+1);
	return string (IDs+rPos);
}

float SEinventory::getRecordMeas(int record, int meas) const{
	return Measurements[meas][record];


}

string SEinventory::getElementType(void) const{
	return elementType;
}

string SEinventory::getElementSide(void) const{
	return elementSide;
}

