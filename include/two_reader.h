#ifndef TWO_READER_H_
#define TWO_READER_H_

#include <fstream>
#include <thread>

#include "buffer.h"
#include "core.h"
#include "index.h"
#include "zstd_codec.h"
#include "vcf_utils.h"
#include "intervals.h"

namespace tomahawk {

class twk_two_filter {
public:

	typedef bool(twk_two_filter::*filter_func)(const twk1_two_t* rec) const;

	twk_two_filter() :
		filter_vec(0),
		minR2(0), maxR2(100), minR(-100), maxR(100), minD(-100), maxD(100),
		minDprime(0), maxDprime(100), minP(0), maxP(1),
		hA_min(0), hA_max(999999999), hB_min(0), hB_max(999999999),
		hC_min(0), hC_max(999999999), hD_min(0), hD_max(999999999),
		mhc_min(0), mhc_max(999999999),
		minChi(0), maxChi(std::numeric_limits<double>::max()),
		minChiModel(0), maxChiModel(std::numeric_limits<double>::max()),
		flag_include(std::numeric_limits<uint32_t>::max()), flag_exclude(0)
	{}

	inline twk_two_filter& SetR2(const double low, const double high){ minR2 = low; maxR2 = high; filter_vec |= 1 << 0; return(*this); }
	inline twk_two_filter& SetR2Low(const double low){ minR2 = low; filter_vec |= 1 << 0; return(*this); }
	inline twk_two_filter& SetR2High(const double high){ maxR2 = high; filter_vec |= 1 << 0; return(*this); }

	inline twk_two_filter& SetD(const double low, const double high){ minD = low; maxD = high; filter_vec |= 1 << 1; return(*this); }
	inline twk_two_filter& SetDLow(const double low){ minD = low; filter_vec |= 1 << 1; return(*this); }
	inline twk_two_filter& SetDHigh(const double high){ maxD = high; filter_vec |= 1 << 1; return(*this); }

	inline twk_two_filter& SetDprime(const double low, const double high){ minDprime = low; maxDprime = high; filter_vec |= 1 << 2; return(*this); }
	inline twk_two_filter& SetDprimeLow(const double low){ minDprime = low; filter_vec |= 1 << 2; return(*this); }
	inline twk_two_filter& SetDprimeHigh(const double high){ maxDprime = high; filter_vec |= 1 << 2; return(*this); }

	inline twk_two_filter& SetP(const double low, const double high){ minP = low; maxP = high; filter_vec |= 1 << 3; return(*this); }
	inline twk_two_filter& SetPLow(const double low){ minP = low; filter_vec |= 1 << 3; return(*this); }
	inline twk_two_filter& SetPHigh(const double high){ maxP = high; filter_vec |= 1 << 3; return(*this); }

	inline twk_two_filter& SetHapA(const double low, const double high){ hA_min = low; hA_max = high; filter_vec |= 1 << 4; return(*this); }
	inline twk_two_filter& SetHapALow(const double low){ hA_min = low; filter_vec |= 1 << 4; return(*this); }
	inline twk_two_filter& SetHapAHigh(const double high){ hA_max = high; filter_vec |= 1 << 4; return(*this); }

	inline twk_two_filter& SetHapB(const double low, const double high){ hB_min = low; hB_max = high; filter_vec |= 1 << 5; return(*this); }
	inline twk_two_filter& SetHapBLow(const double low){ hB_min = low; filter_vec |= 1 << 5; return(*this); }
	inline twk_two_filter& SetHapBHigh(const double high){ hB_max = high; filter_vec |= 1 << 5; return(*this); }

	inline twk_two_filter& SetHapC(const double low, const double high){ hC_min = low; hC_max = high; filter_vec |= 1 << 6; return(*this); }
	inline twk_two_filter& SetHapCLow(const double low){ hC_min = low; filter_vec |= 1 << 6; return(*this); }
	inline twk_two_filter& SetHapCHigh(const double high){ hC_max = high; filter_vec |= 1 << 6; return(*this); }

	inline twk_two_filter& SetHapD(const double low, const double high){ hD_min = low; hD_max = high; filter_vec |= 1 << 7; return(*this); }
	inline twk_two_filter& SetHapDLow(const double low){ hD_min = low; filter_vec |= 1 << 7; return(*this); }
	inline twk_two_filter& SetHapDHigh(const double high){ hD_max = high; filter_vec |= 1 << 7; return(*this); }

	inline twk_two_filter& SetR(const double low, const double high){ minR = low; maxR = high; filter_vec |= 1 << 8; return(*this); }
	inline twk_two_filter& SetRLow(const double low){ minR = low; filter_vec |= 1 << 8; return(*this); }
	inline twk_two_filter& SetRHigh(const double high){ maxR = high; filter_vec |= 1 << 8; return(*this); }

	inline twk_two_filter& SetUpperTrig(){ filter_vec |= 1 << 9; return(*this); }
	inline twk_two_filter& SetLowerTrig(){ filter_vec |= 1 << 10; return(*this); }

	inline twk_two_filter& SetMHC(const double low, const double high){ mhc_min = low; mhc_max = high; filter_vec |= 1 << 11; return(*this); }
	inline twk_two_filter& SetMHCLow(const double low){ mhc_min = low; filter_vec |= 1 << 11; return(*this); }
	inline twk_two_filter& SetMHCHigh(const double high){ mhc_max = high; filter_vec |= 1 << 11; return(*this); }

	inline twk_two_filter& SetFlag(const double low, const double high){ flag_include = low; flag_exclude = high; filter_vec |= 1 << 12; return(*this); }
	inline twk_two_filter& SetFlagInclude(const double low){ flag_include = low; filter_vec |= 1 << 12; return(*this); }
	inline twk_two_filter& SetFlagExclude(const double high){ flag_exclude = high; filter_vec |= 1 << 12; return(*this); }

	inline twk_two_filter& SetChiSq(const double low, const double high){ minChi = low; maxChi = high; filter_vec |= 1 << 13; return(*this); }
	inline twk_two_filter& SetChiSqLow(const double low){ minChi = low; filter_vec |= 1 << 13; return(*this); }
	inline twk_two_filter& SetChiSqHigh(const double high){ maxChi = high; filter_vec |= 1 << 13; return(*this); }

	inline twk_two_filter& SetChiSqModel(const double low, const double high){ minChiModel = low; maxChiModel = high; filter_vec |= 1 << 14; return(*this); }
	inline twk_two_filter& SetChiSqModelLow(const double low){ minChiModel = low; filter_vec |= 1 << 14; return(*this); }
	inline twk_two_filter& SetChiSqModelHigh(const double high){ maxChiModel = high; filter_vec |= 1 << 14; return(*this); }

	/**<
	 * Constructs the list of functional pointers to the actual functons to
	 * use for filtering records. It is mandatory to invoke this function
	 * prior invoking `Filter`.
	 * @return Returns TRUE.
	 */
	bool Build(){
		funcs.clear();
		if((filter_vec >> 0)  & 1) funcs.push_back(&twk_two_filter::FilterR2);
		if((filter_vec >> 1)  & 1) funcs.push_back(&twk_two_filter::FilterD);
		if((filter_vec >> 2)  & 1) funcs.push_back(&twk_two_filter::FilterDprime);
		if((filter_vec >> 3)  & 1) funcs.push_back(&twk_two_filter::FilterP);
		if((filter_vec >> 4)  & 1) funcs.push_back(&twk_two_filter::FilterHapA);
		if((filter_vec >> 5)  & 1) funcs.push_back(&twk_two_filter::FilterHapB);
		if((filter_vec >> 6)  & 1) funcs.push_back(&twk_two_filter::FilterHapC);
		if((filter_vec >> 7)  & 1) funcs.push_back(&twk_two_filter::FilterHapD);
		if((filter_vec >> 8)  & 1) funcs.push_back(&twk_two_filter::FilterR);
		if((filter_vec >> 9)  & 1) funcs.push_back(&twk_two_filter::FilterUpperTrig);
		if((filter_vec >> 10) & 1) funcs.push_back(&twk_two_filter::FilterLowerTrig);
		if((filter_vec >> 11) & 1) funcs.push_back(&twk_two_filter::FilterMHC);
		if((filter_vec >> 12) & 1) funcs.push_back(&twk_two_filter::FilterFlags);
		if((filter_vec >> 13) & 1) funcs.push_back(&twk_two_filter::FilterChiSq);
		if((filter_vec >> 14) & 1) funcs.push_back(&twk_two_filter::FilterChiSqModel);
		return true;
	}

	bool Filter(const twk1_two_t* rec) const{
		if(funcs.size() == 0) return true;
		for(int i = 0; i < funcs.size(); ++i){
			if((this->*funcs[i])(rec) == false)
				return false;
		}
		return true;
	}

	inline bool FilterR2(const twk1_two_t* rec) const{ return(rec->R2 >= minR2 && rec->R2 <= maxR2); }
	inline bool FilterD(const twk1_two_t* rec) const{ return(rec->D >= minD && rec->D <= maxD); }
	inline bool FilterDprime(const twk1_two_t* rec) const{ return(rec->Dprime >= minDprime && rec->Dprime <= maxDprime); }
	inline bool FilterP(const twk1_two_t* rec) const{ return(rec->P >= minP && rec->P <= maxP); }
	inline bool FilterHapA(const twk1_two_t* rec) const{ return(rec->cnt[0] >= hA_min && rec->cnt[0] <= hA_max); }
	inline bool FilterHapB(const twk1_two_t* rec) const{ return(rec->cnt[1] >= hB_min && rec->cnt[1] <= hB_max); }
	inline bool FilterHapC(const twk1_two_t* rec) const{ return(rec->cnt[2] >= hC_min && rec->cnt[2] <= hC_max); }
	inline bool FilterHapD(const twk1_two_t* rec) const{ return(rec->cnt[3] >= hD_min && rec->cnt[3] <= hD_max); }
	inline bool FilterR(const twk1_two_t* rec) const{ return(rec->R >= minR2 && rec->R <= maxR2); }
	inline bool FilterUpperTrig(const twk1_two_t* rec) const{
		return(rec->ridA <= rec->ridB && (rec->ridA == rec->ridB && rec->Apos < rec->Bpos));
	}
	inline bool FilterLowerTrig(const twk1_two_t* rec) const{
		return(rec->ridB <= rec->ridA && (rec->ridA == rec->ridB && rec->Bpos < rec->Apos));
	}

	bool FilterMHC(const twk1_two_t* rec) const{
		const double* m;
		m = rec->cnt[0] > rec->cnt[1] ? &rec->cnt[0] : &rec->cnt[1];
		m = rec->cnt[2] > *m ? &rec->cnt[2] : m;
		m = rec->cnt[3] > *m ? &rec->cnt[3] : m;
		double c = 0;
		for(int i = 0; i < 4; ++i) c += (m == &rec->cnt[i] ? 0 : rec->cnt[i]);
		return(c >= mhc_min && c<= mhc_max);
	}

	inline bool FilterFlags(const twk1_two_t* rec) const{
		return((rec->controller & flag_include) && ((rec->controller & flag_exclude) == 0));
	}

	inline bool FilterChiSq(const twk1_two_t* rec) const{ return(rec->ChiSqFisher >= minChi && rec->ChiSqFisher <= maxChi); }
	inline bool FilterChiSqModel(const twk1_two_t* rec) const{ return(rec->ChiSqModel >= minChiModel && rec->ChiSqModel <= maxChiModel); }

public:
	uint32_t filter_vec;
	double minR2, maxR2, minR, maxR, minD, maxD, minDprime, maxDprime, minP, maxP;
	double hA_min, hA_max, hB_min, hB_max, hC_min, hC_max, hD_min, hD_max;
	double mhc_min, mhc_max;
	double minChi, maxChi, minChiModel, maxChiModel;
	uint32_t flag_include, flag_exclude;
private:
	std::vector<filter_func> funcs;
};

class twk1_two_iterator {
public:
	twk1_two_iterator() : offset(0), stream(nullptr), rcd(nullptr){}
	~twk1_two_iterator(){ }

	bool NextBlockRaw();
	bool NextBlock();
	bool NextRecord();
	inline const twk1_two_block_t& GetBlock(void) const{ return(this->blk); }

public:
	uint64_t offset;
	ZSTDCodec zcodec; // support codec
	twk_buffer_t buf; // support buffer
	twk_oblock_two_t oblk;// block wrapper
	twk1_two_block_t blk; // block
	std::istream* stream; // stream pointer
	twk1_two_t* rcd;
};

struct twk_two_settings {
	twk_two_settings() : out_type('b'), n_threads(std::thread::hardware_concurrency()){}

	char out_type;
	int32_t n_threads;
	std::string in, out;
	std::vector<std::string> ivals; // unparsed interval strings
	twk_intervals_two intervals;
	twk_two_filter filter;
};

/**<
 * Reader of twk output files.
 */
class two_reader {
public:
	two_reader() : buf(nullptr), stream(nullptr){}
	~two_reader(){ delete stream; }

	/**<
	 * Open a target two file. File header, index, and footer will be read
	 * and parsed as part of the opening procedure. If these pass without errors
	 * then return TRUE. Returns FALSE otherwise.
	 * @param file Input target two file.
	 * @return     Returns TRUE upon success or FALSE otherwise.
	 */
	bool Open(std::string file);

	inline bool NextBlock(){ return(it.NextBlock()); }
	inline bool NextBlockRaw(){ return(it.NextBlockRaw()); }
	inline bool NextRecord(){ return(it.NextRecord()); }

public:
	std::streambuf* buf;
	std::istream* stream;
	std::ifstream fstream;
	io::VcfHeader hdr;
	IndexOutput index;

	twk_two_filter filter;
	twk1_two_iterator it;
	ZSTDCodec zcodec;
};

}

#endif /* TWO_READER_H_ */
