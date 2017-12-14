/*
 * CSVIterator.h
 *
 *  Created on: Jan 14, 2017
 *      Author: Julia
 */

#ifndef CSVITERATOR_H_
#define CSVITERATOR_H_

class CSVIterator {
public:
	typedef std::input_iterator_tag     iterator_category;
	typedef CSVRow                      value_type;
	typedef size_t                 difference_type;
	typedef CSVRow*                     pointer;
	typedef CSVRow&                     reference;
	CSVIterator(istream&);
	CSVIterator();
	virtual ~CSVIterator();

	CSVIterator& operator++();
	CSVIterator operator++(int);
	CSVRow const& operator*() const;
	CSVRow const* operator->() const;
	bool operator==(CSVIterator const&);
	bool operator!=(CSVIterator const& rhs);
private:
	istream*       m_str;
	CSVRow         m_row;
};

#endif /* CSVITERATOR_H_ */
