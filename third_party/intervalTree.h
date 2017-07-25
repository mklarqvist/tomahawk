/*

Copyright (c) 2011 Erik Garrison

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
#ifndef __INTERVAL_TREE_H
#define __INTERVAL_TREE_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <memory>

namespace Tomahawk{
namespace Algorithm{

template <class T, typename K = std::size_t>
class Interval {
public:
    K start;
    K stop;
    T value;
    Interval(K s, K e, const T& v)
        : start(s)
        , stop(e)
        , value(v)
    { }
};

class ContigInterval : public Interval<ContigInterval*, S32> {
public:
	ContigInterval() : Interval<ContigInterval*, S32>(-1, -1, nullptr), contigID(-1){}

    void operator()(S32 s, S32 e, S32 id){
    	this->start = s;
    	this->stop  = e;
    	this->contigID = id;
    }

    friend std::ostream& operator<<(std::ostream& os, const ContigInterval& i){
    	os << i.contigID << ":" << i.start << "-" << i.stop << ';' << i.value;
    	return(os);
    }

public:
    S32 contigID;
};

template <class T, typename K>
K intervalStart(const Interval<T,K>& i) {
    return i.start;
}

template <class T, typename K>
K intervalStop(const Interval<T,K>& i) {
    return i.stop;
}

template <class T, typename K>
std::ostream& operator<<(std::ostream& out, const Interval<T,K>& i) {
	out << "Interval(" << i.start << ", " << i.stop << "): " << i.value;
	return out;
}

template <class T, typename K = std::size_t>
class IntervalStartSorter {
public:
    bool operator() (const Interval<T,K>& a, const Interval<T,K>& b) {
        return a.start < b.start;
    }
};

template <class T, typename K = std::size_t>
class IntervalTree {

public:
    typedef T interval;
    typedef std::vector<interval> intervalVector;
    typedef IntervalTree<T,K> intervalTree;

    intervalVector intervals;
    std::unique_ptr<intervalTree> left;
    std::unique_ptr<intervalTree> right;
    K center;

    IntervalTree<T,K>(void)
        : left(nullptr)
        , right(nullptr)
        , center(0)
    { }

private:
    std::unique_ptr<intervalTree> copyTree(const intervalTree& orig){
        return std::unique_ptr<intervalTree>(new intervalTree(orig));
    }

public:

    IntervalTree<T,K>(const intervalTree& other)
    :   intervals(other.intervals),
        left(other.left ? copyTree(*other.left) : nullptr),
        right(other.right ? copyTree(*other.right) : nullptr),
        center(other.center)
    {
    }

public:
    IntervalTree<T,K>& operator=(const intervalTree& other) {
        center = other.center;
        intervals = other.intervals;
        left = other.left ? copyTree(*other.left) : nullptr;
        right = other.right ? copyTree(*other.right) : nullptr;
        return *this;
    }

    // Note: changes the order of ivals
    IntervalTree<T,K>(
            intervalVector& ivals,
            std::size_t depth = 16,
            std::size_t minbucket = 64,
            K leftextent = 0,
            K rightextent = 0,
            std::size_t maxbucket = 512
            )
        : left(nullptr)
        , right(nullptr)
    {

        --depth;
        IntervalStartSorter<T,K> intervalStartSorter;
        if (depth == 0 || (ivals.size() < minbucket && ivals.size() < maxbucket)) {
            std::sort(ivals.begin(), ivals.end(), intervalStartSorter);
            intervals = ivals;
        } else {
            if (leftextent == 0 && rightextent == 0) {
                // sort intervals by start
              std::sort(ivals.begin(), ivals.end(), intervalStartSorter);
            }

            K leftp = 0;
            K rightp = 0;
            K centerp = 0;

            if (leftextent || rightextent) {
                leftp = leftextent;
                rightp = rightextent;
            } else {
                leftp = ivals.front().start;
                std::vector<K> stops;
                stops.resize(ivals.size());
                transform(ivals.begin(), ivals.end(), stops.begin(), intervalStop<T,K>);
                rightp = *max_element(stops.begin(), stops.end());
            }

            //centerp = ( leftp + rightp ) / 2;
            centerp = ivals.at(ivals.size() / 2).start;
            center = centerp;

            intervalVector lefts;
            intervalVector rights;

            for (typename intervalVector::const_iterator i = ivals.begin(); i != ivals.end(); ++i) {
                const interval& interval = *i;
                if (interval.stop < center) {
                    lefts.push_back(interval);
                } else if (interval.start > center) {
                    rights.push_back(interval);
                } else {
                    intervals.push_back(interval);
                }
            }

            if (!lefts.empty()) {
                left = std::unique_ptr<intervalTree>(new intervalTree(lefts, depth, minbucket, leftp, centerp));
            }
            if (!rights.empty()) {
                right = std::unique_ptr<intervalTree>(new intervalTree(rights, depth, minbucket, centerp, rightp));
            }
        }
    }

    intervalVector findOverlapping(K start, K stop) const {
		intervalVector ov;
		this->findOverlapping(start, stop, ov);
		return ov;
    }

    bool findAny(K start, K stop) const{
    	if (!intervals.empty() && ! (stop < intervals.front().start)) {
			for (typename intervalVector::const_iterator i = intervals.begin(); i != intervals.end(); ++i) {
				const interval& interval = *i;
				if (interval.stop >= start && interval.start <= stop)
					return true;
			}
		}

    	if (left && start <= center) {
			left->findOverlapping(start, stop);
		}

		if (right && stop >= center) {
			right->findOverlapping(start, stop);
		}

		return false;
    }

    void findOverlapping(K start, K stop, intervalVector& overlapping) const {
        if (!intervals.empty() && ! (stop < intervals.front().start)) {
            for (typename intervalVector::const_iterator i = intervals.begin(); i != intervals.end(); ++i) {
                const interval& interval = *i;
                if (interval.stop >= start && interval.start <= stop) {
                    overlapping.push_back(interval);
                }
            }
        }

        if (left && start <= center) {
            left->findOverlapping(start, stop, overlapping);
        }

        if (right && stop >= center) {
            right->findOverlapping(start, stop, overlapping);
        }

    }

    intervalVector findContained(K start, K stop) const {
	intervalVector contained;
	this->findContained(start, stop, contained);
	return contained;
    }

    void findContained(K start, K stop, intervalVector& contained) const {
        if (!intervals.empty() && ! (stop < intervals.front().start)) {
            for (typename intervalVector::const_iterator i = intervals.begin(); i != intervals.end(); ++i) {
                const interval& interval = *i;
                if (interval.start >= start && interval.stop <= stop) {
                    contained.push_back(interval);
                }
            }
        }

        if (left && start <= center) {
            left->findContained(start, stop, contained);
        }

        if (right && stop >= center) {
            right->findContained(start, stop, contained);
        }

    }

    ~IntervalTree(void) = default;

};

}
}

#endif
