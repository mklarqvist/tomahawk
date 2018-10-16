/*
Copyright (C) 2017-current Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
==============================================================================*/
#ifndef TOMAHAWK_GENERIC_ITERATOR_H_
#define TOMAHAWK_GENERIC_ITERATOR_H_

#include <cstddef>
#include <iterator>

namespace tomahawk{

/**<
 * Raw iterator with random access.
 */
template<typename DataType>
class yonRawIterator : public std::iterator<std::random_access_iterator_tag,
                                            DataType, ptrdiff_t, DataType*, DataType&>
{
public:
    typedef yonRawIterator       self_type;
    typedef DataType             value_type;
    typedef value_type&          reference;
    typedef const value_type&    const_reference;
    typedef value_type*          pointer;
    typedef const value_type*    const_pointer;
    typedef std::ptrdiff_t       difference_type;
    typedef std::size_t          size_type;

public:
	yonRawIterator(DataType* ptr = nullptr){m_ptr = ptr;}
	yonRawIterator(const self_type& rawIterator) = default;
    ~yonRawIterator(){}

    self_type& operator=(const self_type& rawIterator) = default;
    self_type& operator=(pointer ptr){ m_ptr = ptr; return (*this); }

    operator bool() const {
        if(m_ptr) return true;
        else return false;
    }

    bool operator==(const self_type& rawIterator)const{ return (m_ptr == rawIterator.getConstPtr()); }
    bool operator!=(const self_type& rawIterator)const{ return (m_ptr != rawIterator.getConstPtr()); }

    self_type& operator+=(const ptrdiff_t& movement){ m_ptr += movement;return (*this); }
    self_type& operator-=(const ptrdiff_t& movement){ m_ptr -= movement;return (*this); }
    self_type& operator++(){ ++m_ptr;return (*this); }
    self_type& operator--(){ --m_ptr;return (*this); }
    self_type  operator++(int){ auto temp(*this);++m_ptr;return temp; }
    self_type  operator--(int){ auto temp(*this);--m_ptr;return temp; }
    self_type  operator+(const ptrdiff_t& movement){ auto oldPtr = m_ptr;m_ptr+=movement; auto temp(*this);m_ptr = oldPtr; return temp; }
    self_type  operator-(const ptrdiff_t& movement){ auto oldPtr = m_ptr;m_ptr-=movement; auto temp(*this);m_ptr = oldPtr; return temp; }
    ptrdiff_t operator-(const self_type& rawIterator){ return std::distance(rawIterator.getPtr(),this->getPtr()); }
    reference operator*(){ return *m_ptr; }
    const_reference operator*()const{ return *m_ptr; }
    pointer operator->(){ return m_ptr; }
    pointer getPtr()const{ return m_ptr; }
    const_pointer getConstPtr()const{ return m_ptr; }

protected:
    pointer m_ptr;
};

}

#endif
