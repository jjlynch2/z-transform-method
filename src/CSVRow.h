/*
 * CSVRow.h
 *
 *  Created on: Jan 14, 2017
 *      Author: Julia
 */

#ifndef CSVROW_H_
#define CSVROW_H_

class CSVRow {
public:
	CSVRow();
	virtual ~CSVRow();
	string const & operator[] (size_t) const;
	size_t size() const;
	void readNextRow(istream &);
	friend istream & operator>>(istream&, CSVRow&);
private:
        vector<string>    m_data;
};

#endif /* CSVROW_H_ */
