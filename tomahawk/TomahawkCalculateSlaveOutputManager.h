#ifndef TOMAHAWK_TOMAHAWKCALCULATESLAVEOUTPUTMANAGER_H_
#define TOMAHAWK_TOMAHAWKCALCULATESLAVEOUTPUTMANAGER_H_

#include "TomahawkBlockManager.h"

//#define SLAVE_FLUSH_LIMIT	65536
#define SLAVE_FLUSH_LIMIT	10000000
//#define SLAVE_FLUSH_LIMIT	5000

namespace Tomahawk{

template <class T>
struct TomahawkCalculateSlaveOutputManager{
	typedef IO::GenericWriterInterace writer_type;
	typedef TomahawkBlock<const T> controller_type;
	typedef TomahawkCalculateSlaveOutputManager<T> self_type;
	typedef Support::TomahawkOutputLD helper_type;
	typedef IO::BasicBuffer buffer_type;

	// Function pointer to write class function
	typedef void (self_type::*outFunction)(const controller_type& a, const controller_type& b, const helper_type& helper);

public:
	TomahawkCalculateSlaveOutputManager(writer_type& writer,
										const writer_type::compression type) :
		writer_output_type(type),
		outCount(0),
		progressCount(0),
		function(type == writer_type::compression::natural ? &self_type::AddNatural : &self_type::AddBinary),
		writer(writer),
		buffer(2*SLAVE_FLUSH_LIMIT),
		sprintf_buffer(new char[255])
	{}

	~TomahawkCalculateSlaveOutputManager(){
		this->Finalise();
		this->buffer.deleteAll();
		delete [] this->sprintf_buffer;
	}

	inline void Add(const controller_type& a, const controller_type& b, const helper_type& helper){ (this->*function)(a, b, helper); }
	void Finalise(void){
		this->writer << this->buffer;
		this->writer.flush();
		this->buffer.reset();
	}

	bool writeHeaders(void){
		if(this->writer_output_type == writer_type::compression::natural){
			// Write tab-delimited header
		} else if(this->writer_output_type == writer_type::compression::binary) {
			// Write binary header
			// Write index as we progress
		} else{
			std::cerr << Helpers::timestamp("ERROR", "WRITER") << "Illegal compression method: " << (int)this->writer_output_type << "..." << std::endl;
			return false;
		}

		return true;
	}

	inline const U64& GetCounts(void) const{ return this->outCount; }
	inline void ResetProgress(void){ this->progressCount = 0; }
	inline const U32& GetProgressCounts(void) const{ return this->progressCount; }

	template <class K>
	self_type& operator<<(const K& data){
		this->buffer += data;
		return(*this);
	}

	void flushIfFull(void){
		if(this->buffer.size() > SLAVE_FLUSH_LIMIT){
			// Move to writer
			//std::cout.write(&this->outstreamBuffer.data[0], this->outstreamBuffer.size());
			this->writer << this->buffer;
			this->buffer.reset();
		}
	}

private:
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

		this->buffer += std::to_string(a.support->contigID);
		this->buffer += '\t';
		this->buffer += std::to_string(a.meta[a.metaPointer].position);
		this->buffer += '\t';
		this->buffer += std::to_string(b.support->contigID);
		this->buffer += '\t';
		this->buffer += std::to_string(b.meta[b.metaPointer].position);
		this->buffer += '\t';
		this->buffer += std::to_string(helper.alleleCounts[0]);
		this->buffer += '\t';
		this->buffer += std::to_string(helper.alleleCounts[1]);
		this->buffer += '\t';
		this->buffer += std::to_string(helper.alleleCounts[4]);
		this->buffer += '\t';
		this->buffer += std::to_string(helper.alleleCounts[5]);
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

		if(this->buffer.size() > SLAVE_FLUSH_LIMIT){
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
			this->writer << this->buffer;
			this->buffer.reset();
		}
	}

	const writer_type::compression writer_output_type;
	U64 outCount;
	U32 progressCount;
	outFunction function;
	writer_type& writer;
	buffer_type buffer;
	char* sprintf_buffer;
};


}

#endif /* TOMAHAWK_TOMAHAWKCALCULATESLAVEOUTPUTMANAGER_H_ */
