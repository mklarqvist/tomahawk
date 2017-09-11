#ifndef TOMAHAWKOUTPUTSORTMERGEQUEUECONTAINER_H_
#define TOMAHAWKOUTPUTSORTMERGEQUEUECONTAINER_H_

namespace Tomahawk{
namespace Algorithm{
namespace Output{

template <class T>
struct TomahawkOutputSortMergeQueueContainer {
	typedef T entry_type;
	typedef TomahawkOutputSortMergeQueueContainer<entry_type> self_type;

public:
	TomahawkOutputSortMergeQueueContainer(const entry_type& data,
										  U32 streamID,
										  bool (*compFunc)(const entry_type& a, const entry_type& b) = T::operator<)
	: streamID(streamID)
	, data(data)
	, compFunc(compFunc)
    {}

    bool operator<(const self_type& a) const{
        return!(this->compFunc(this->data, a.data));
    }

public:
    U32 streamID;
    entry_type data;
    bool (*compFunc)(const entry_type& a, const entry_type& b);
};

}
}
}

#endif /* TOMAHAWKOUTPUTSORTMERGEQUEUECONTAINER_H_ */
