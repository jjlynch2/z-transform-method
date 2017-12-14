/*
 * CSVIterator.cpp
 *
 *  Created on: Jan 15, 2017
 *      Author: Julia
 */
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
using namespace std;
#include "CSVRow.h"
#include "CSVIterator.h"

CSVIterator::CSVIterator(istream& str) :m_str(str.good()?&str:NULL) {
	++(*this);
}

CSVIterator::CSVIterator() :m_str(NULL) {}

CSVIterator::~CSVIterator() {
	// TODO Auto-generated destructor stub
}

CSVIterator& CSVIterator::operator++(){
	if (m_str)
	{
		if (!((*m_str) >> m_row))
		{
			m_str = NULL;
		}
	}
	return *this;
}

CSVIterator CSVIterator::operator++(int){

	CSVIterator  tmp(*this);++(*this);
	return tmp;
}

CSVRow const& CSVIterator::operator*() const
{
	return m_row;
}

CSVRow const* CSVIterator::operator->() const
{
	return &m_row;
}

bool CSVIterator::operator==(CSVIterator const& rhs) {
	return ((this == &rhs) || ((this->m_str == NULL) && (rhs.m_str == NULL)));
}

bool CSVIterator::operator!=(CSVIterator const& rhs) {
	return !((*this) == rhs);
}
