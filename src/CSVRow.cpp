/*
 * CSVRow.cpp
 *
 *  Created on: Jan 14, 2017
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

CSVRow::CSVRow() {
	// TODO Auto-generated constructor stub

}

CSVRow::~CSVRow() {
	// TODO Auto-generated destructor stub
}

string const & CSVRow::operator[](size_t index) const {
            return m_data[index];
}

size_t CSVRow::size() const
{
    return m_data.size();
}

void CSVRow::readNextRow(std::istream& str) {
            std::string         line;
            std::getline(str, line);

            std::stringstream   lineStream(line);
            std::string         cell;

            m_data.clear();
            while(std::getline(lineStream, cell, '\t'))
            {
                m_data.push_back(cell);
            }
            // This checks for a trailing comma with no data after it.
            if (!lineStream && cell.empty())
            {
                // If there was a trailing comma then add an empty element.
                m_data.push_back(string(""));
            }
}

istream& operator>>(istream& str, CSVRow& data)
{
    data.readNextRow(str);
    return str;
}
