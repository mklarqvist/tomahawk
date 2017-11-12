#ifndef TOTEMPOLEENTRY_H_
#define TOTEMPOLEENTRY_H_

namespace Tomahawk{
namespace Totempole{

struct TotempoleEntry{
	typedef TotempoleEntry self_type;

public:
	TotempoleEntry() :
		byte_offset(0),
		byte_offset_end(0),
		minPosition(0),
		maxPosition(0),
		contigID(0),
		n_variants(0),
		l_meta(0),
		l_meta_complex(0),
		l_gt_rle(0),
		l_gt_simple(0),
		n_info_patterns(0),
		n_format_patterns(0),
		n_filter_patterns(0),
		n_info_streams(0),
		n_format_streams(0),
		n_filter_streams(0)
	{}
	~TotempoleEntry(){}

	inline bool isValid(void) const{ return(this->byte_offset != 0); }

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << entry.byte_offset << '\t' << entry.byte_offset_end << '\t' << entry.contigID << '\t' <<
				  entry.minPosition << '-' << entry.maxPosition << '\t' << entry.n_variants << '\t' <<
				  entry.l_meta << '\t' <<
				  entry.l_gt_rle << '\t' << entry.l_gt_simple << '\t' << entry.l_meta_complex;
		return(stream);
	}

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.byte_offset_end), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.minPosition), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.maxPosition), sizeof(U64));
		stream.write(reinterpret_cast<const char*>(&entry.contigID),    sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.n_variants),  sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&entry.l_meta), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.l_meta_complex), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.l_gt_rle), sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.l_gt_simple), sizeof(U32));

		stream.write(reinterpret_cast<const char*>(&entry.n_info_patterns), sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&entry.n_format_patterns), sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&entry.n_filter_patterns), sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&entry.n_info_streams), sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&entry.n_format_streams), sizeof(U16));
		stream.write(reinterpret_cast<const char*>(&entry.n_filter_streams), sizeof(U16));


		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.byte_offset), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.byte_offset_end), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.minPosition), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.maxPosition), sizeof(U64));
		stream.read(reinterpret_cast<char*>(&entry.contigID),    sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.n_variants),       sizeof(U16));
		stream.read(reinterpret_cast<char*>(&entry.l_meta), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.l_meta_complex), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.l_gt_rle), sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.l_gt_simple), sizeof(U32));

		stream.read(reinterpret_cast<char*>(&entry.n_info_patterns), sizeof(U16));
		stream.read(reinterpret_cast<char*>(&entry.n_format_patterns), sizeof(U16));
		stream.read(reinterpret_cast<char*>(&entry.n_filter_patterns), sizeof(U16));
		stream.read(reinterpret_cast<char*>(&entry.n_info_streams), sizeof(U16));
		stream.read(reinterpret_cast<char*>(&entry.n_format_streams), sizeof(U16));
		stream.read(reinterpret_cast<char*>(&entry.n_filter_streams), sizeof(U16));

		return(stream);
	}

	void reset(void){
		this->byte_offset = 0;
		this->byte_offset_end = 0;
		this->contigID = 0;
		this->minPosition = 0;
		this->maxPosition = 0;
		this->n_variants = 0;
		this->l_meta = 0;
		this->l_gt_rle = 0;
		this->l_gt_simple = 0;
		this->l_meta_complex = 0;
		this->n_info_patterns = 0;
		this->n_format_patterns = 0;
		this->n_filter_patterns = 0;
		this->n_info_streams = 0;
		this->n_format_streams = 0;
		this->n_filter_streams = 0;
	}

public:
	U64 byte_offset;	// tellg() position in stream for start of record in Tomahawk file
	U64 byte_offset_end;// tellg() position in stream for start of record in Tomahawk file
	U64 minPosition;	// smallest bp position
	U64 maxPosition;	// largest bp position
	S32 contigID;		// contig identifier
	U16 n_variants;     // number of variants in this block
	U32 l_meta;
	U32 l_meta_complex;
	U32 l_gt_rle;
	U32 l_gt_simple;

	// TODO:
	U16 n_info_patterns;
	U16 n_format_patterns;
	U16 n_filter_patterns;
	U16 n_info_streams;
	U16 n_format_streams;
	U16 n_filter_streams;
	// INFO patterns
	// FORMAT patterns
	// FILTER patterns
	// Each INFO/FORMAT/FILTER stream_offset, type-width, controller byte

	// PPA length and offset are implicit
	//T* PPA;
};

}
}

#endif /* TOTEMPOLEENTRY_H_ */
