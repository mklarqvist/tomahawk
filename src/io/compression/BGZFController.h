#ifndef BGZFCONTROLLER_H_
#define BGZFCONTROLLER_H_

#include "../../support/helpers.h"
#include "../../third_party/zlib/zconf.h"
#include "../../third_party/zlib/zlib.h"
#include "../BasicBuffer.h"
#include "GZFHeader.h"

namespace Tomahawk {
namespace IO {

class BGZFController {
	typedef BGZFController self_type;
		typedef IO::BasicBuffer buffer_type;
		typedef BGZFHeader header_type;

	public:
		BGZFController();
		BGZFController(const char* data, const U32 length);
		~BGZFController();

		void Clear();
		bool Inflate(buffer_type& input, buffer_type& output, const header_type& header) const;
		bool Inflate(buffer_type& input, buffer_type& output) const;
		U32 InflateSize(buffer_type& input) const;
		bool InflateBlock(std::ifstream& stream, buffer_type& input);

		friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
			stream.write(entry.buffer.data, entry.buffer.pointer);
			return stream;
		}

	private:
		bool __Inflate(buffer_type& input, buffer_type& output, const header_type& header) const;

	public:
		buffer_type buffer;
	};

} /* namespace IO */
} /* namespace Tomahawk */

#endif /* BGZFCONTROLLER_H_ */
