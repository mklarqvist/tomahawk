#ifndef TOMAHAWK_TOMAHAWKOUTPUTMANAGER_H_
#define TOMAHAWK_TOMAHAWKOUTPUTMANAGER_H_

#include "../io/BasicWriters.h"
#include "TomahawkBlockManager.h"
#include "../io/GZController.h"

#define SLAVE_FLUSH_LIMIT	10000000	// 10 MB default flush limit
#define SLAVE_FLUSH_LIMIT_NATURAL 65536

namespace Tomahawk{
namespace IO {

template <class T>
struct TomahawkOutputManager{
	typedef IO::GenericWriterInterace writer_type;
	typedef TomahawkBlock<const T> controller_type;
	typedef TomahawkOutputManager<T> self_type;
	typedef Support::TomahawkOutputLD helper_type;
	typedef IO::BasicBuffer buffer_type;
	typedef GZController tgzf_controller;

	// Function pointer to write class function
	typedef void (self_type::*outFunction)(const controller_type& a, const controller_type& b, const helper_type& helper);
	typedef void (self_type::*flushFunction)(void);

public:
	TomahawkOutputManager(writer_type& writer,
						  const writer_type::compression type) :
		outCount(0),
		progressCount(0),
		function(type == writer_type::compression::natural ? &self_type::AddNatural : &self_type::AddBinary),
		flush(type == writer_type::compression::natural ? &self_type::FinaliseNatural : &self_type::FinaliseBinary),
		writer(writer),
		buffer(2*SLAVE_FLUSH_LIMIT),
		sprintf_buffer(new char[255])
	{

	}

	~TomahawkOutputManager(){
		this->Finalise();
		this->buffer.deleteAll();
		delete [] this->sprintf_buffer;
	}

	inline void Add(const controller_type& a, const controller_type& b, const helper_type& helper){ (this->*function)(a, b, helper); }
	inline void Finalise(void){ (this->*flush)(); }
	inline const U64& GetCounts(void) const{ return this->outCount; }
	inline void ResetProgress(void){ this->progressCount = 0; }
	inline const U32& GetProgressCounts(void) const{ return this->progressCount; }

private:
	void FinaliseNatural(void){
		this->writer << this->buffer;
		this->buffer.reset();
	}

	void FinaliseBinary(void){
		if(this->buffer.size() > 0){
			if(!this->compressor.Deflate(this->buffer)){
				std::cerr << "failed deflate" << std::endl;
				exit(1);
			}
			this->writer << compressor.buffer_;
			this->buffer.reset();
			this->compressor.Clear();
		}
	}

	void AddNatural(const controller_type& a, const controller_type& b, const helper_type& helper){
		this->buffer += std::to_string(helper.controller);
		this->buffer += '\t';
		this->buffer += std::to_string(-10*(log10(a.meta[a.metaPointer].MAF) + log10(b.meta[b.metaPointer].MAF)));
		this->buffer += '\t';

		/*
		this->buffer += std::to_string(a.meta[a.metaPointer].MAF);
		this->buffer += '\t';
		this->buffer += std::to_string(b.meta[b.metaPointer].MAF);
		this->buffer += '\t';
*/

		//  Todo: need Totempole for mapping contigID to name
		this->buffer += std::to_string(a.support->contigID);
		this->buffer += '\t';
		this->buffer += std::to_string(a.meta[a.metaPointer].position);
		this->buffer += '\t';
		this->buffer += std::to_string(b.support->contigID);
		this->buffer += '\t';
		this->buffer += std::to_string(b.meta[b.metaPointer].position);
		this->buffer += '\t';

		// If data is phased output in integer form
		if((helper.controller & 1) == 1){
			this->buffer += std::to_string((U32)helper.alleleCounts[0]);
			this->buffer += '\t';
			this->buffer += std::to_string((U32)helper.alleleCounts[1]);
			this->buffer += '\t';
			this->buffer += std::to_string((U32)helper.alleleCounts[4]);
			this->buffer += '\t';
			this->buffer += std::to_string((U32)helper.alleleCounts[5]);
		} else {
			this->buffer += std::to_string(helper.alleleCounts[0]);
			this->buffer += '\t';
			this->buffer += std::to_string(helper.alleleCounts[1]);
			this->buffer += '\t';
			this->buffer += std::to_string(helper.alleleCounts[4]);
			this->buffer += '\t';
			this->buffer += std::to_string(helper.alleleCounts[5]);
		}

		this->buffer += '\t';
		this->buffer += std::to_string(helper.D);
		this->buffer += '\t';
		this->buffer += std::to_string(helper.Dprime);
		this->buffer += '\t';
		this->buffer += std::to_string(helper.R2);
		this->buffer += '\t';
		U32 l = std::sprintf(this->sprintf_buffer, "%E", helper.P);
		this->buffer.Add(this->sprintf_buffer, l);
		this->buffer += '\t';
		this->buffer += std::to_string(helper.chiSqFisher);
		this->buffer += '\t';
		this->buffer += std::to_string(helper.chiSqModel);
		this->buffer += '\n';
		++this->outCount;
		++this->progressCount;

		if(this->buffer.size() > SLAVE_FLUSH_LIMIT_NATURAL){
			this->writer << this->buffer;
			this->buffer.reset();
		}
	}

	void AddBinary(const controller_type& a, const controller_type& b, const helper_type& helper){
		const U32 writePosA = a.meta[a.metaPointer].position << 2 | a.meta[a.metaPointer].phased << 1 | a.meta[a.metaPointer].missing;
		const U32 writePosB = b.meta[b.metaPointer].position << 2 | b.meta[b.metaPointer].phased << 1 | b.meta[b.metaPointer].missing;
		this->buffer += helper.controller;
		this->buffer += (double)(-10*(log10(a.meta[a.metaPointer].MAF) + log10(b.meta[b.metaPointer].MAF)));
		this->buffer += a.support->contigID;
		this->buffer += writePosA;
		this->buffer += b.support->contigID;
		this->buffer += writePosB;
		this->buffer << helper;
		++this->outCount;
		++this->progressCount;

		if(this->buffer.size() > SLAVE_FLUSH_LIMIT){
			if(!this->compressor.Deflate(this->buffer)){
				std::cerr << Helpers::timestamp("ERROR", "TGZF") << "Failed deflate..." << std::endl;
				exit(1);
			}
			this->writer << compressor.buffer_;
			this->buffer.reset();
			this->compressor.Clear();
		}
	}

	U64 outCount;			// lines written
	U32 progressCount;		// lines added since last flush
	outFunction function;	// add function pointer
	flushFunction flush;	// flush function pointer
	writer_type& writer;	// writer interface
	buffer_type buffer;		// internal buffer
	tgzf_controller compressor;// compressor
	char* sprintf_buffer;	// special buffer used for sprintf writing scientific output in natural mode
};

}
}

#endif /* TOMAHAWK_TOMAHAWKOUTPUTMANAGER_H_ */
