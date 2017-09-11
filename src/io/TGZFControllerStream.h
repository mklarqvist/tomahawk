#ifndef IO_TGZFCONTROLLERSTREAM_H_
#define IO_TGZFCONTROLLERSTREAM_H_

#include "TGZFController.h"

namespace Tomahawk {
namespace IO {

class TGZFControllerStream : public TGZFController{
protected:
	enum TGZF_STATE {TGZF_OK, TGZF_HEADER, TGZF_END, TGZF_INIT, TGZF_ERROR};

public:
	TGZFControllerStream();
	~TGZFControllerStream();

	bool Inflate(std::ifstream& stream, const BYTE* output, const U32& avail_out, U32& return_size);

	void reset(void){
		this->total_out = 0;
		this->bytes_read = 0;
		this->BSIZE = 0;
		this->STATE = TGZF_STATE::TGZF_INIT;
	}

protected:
	bool InflateOpen(std::ifstream& stream);
	bool __Inflate(std::ifstream& stream, const BYTE* output, const U32 avail_out, U32& return_size);

protected:
	TGZF_STATE STATE;
	U32 chunk_size;
	U32 total_out;
	U32 bytes_read;
	U32 BSIZE;
	z_stream d_stream;
};

}
}



#endif /* IO_TGZFCONTROLLERSTREAM_H_ */
