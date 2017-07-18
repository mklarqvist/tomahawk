/*
POLYGON - Next generation sequencing quality control and preprocessing
Copyright (C) 2015-2016, Marcus D. R. Klarqvist.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/
#ifndef GZCONTROLLER_H_
#define GZCONTROLLER_H_

#include <fstream>

#include "../helpers.h"
#include "../third_party/zlib/zlib.h"
#include "BasicBuffer.h"

namespace Tomahawk{
namespace IO{

class GZController{
	typedef GZController self_type;
	typedef IO::BasicBuffer buffer_type;

public:
	GZController();
	GZController(const char* data, const U32 length);
	GZController(const U32 largest_block_size);
	~GZController();

	void Clear();
	bool Inflate(buffer_type& input, buffer_type& output) const;
	bool Deflate(buffer_type& meta, buffer_type& rle);
	U32 InflateSize(buffer_type& input) const;

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream.write(entry.buffer_.data, entry.buffer_.pointer);
		return stream;
	}

public:
	buffer_type buffer_;
};

}
}

#endif /* GZCONTROLLER_H_ */
