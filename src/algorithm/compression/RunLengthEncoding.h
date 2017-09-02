#ifndef ALGORITHM_COMPRESSION_RUNLENGTHENCODING_H_
#define ALGORITHM_COMPRESSION_RUNLENGTHENCODING_H_

namespace Tomahawk{
namespace Algorithm{

template <class T>
inline int PACK3(const BYTE& ref, char* target, T length){
	if(length <= 7){
		*target++ = 0 | ((ref & 15) << 3) | (length & 7); // highest bit is 0
		return 1;
	}

	char* target0 = target;
	*target++ = 128 | ((ref & 15) << 3) | (length & 7);
	length >>= 3;

	while(true){
		if(length <= 7){
			*target++ = 0 | (length & 127); // highest bit is 0
			length >>= 7;
			break;
		}

		*target++ = 128 | (length & 127);
		length >>= 7;
	}
	return(target - target0);
}

template <class T>
inline int UNPACK3(char* target, T length, BYTE& ref){
	if((*target & 128) == 0){
		length = *target & 7;
		ref = *target >> 3;
		return 1;
	}

	char* target0 = target;
	length = *target & 7;
	ref = (*target >> 3) & 15;
	++target;
	U32 offset = 3;

	while(true){
		length |= (*target & 127) << offset;
		offset += 7;

		if((*target & 128) == 0) break;
		++target;
	}

	return(target - target0);
}

}
}

#endif /* ALGORITHM_COMPRESSION_RUNLENGTHENCODING_H_ */
