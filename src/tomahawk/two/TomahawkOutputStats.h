#ifndef TOMAHAWKOUTPUT_TOMAHAWKOUTPUTSTATS_H_
#define TOMAHAWKOUTPUT_TOMAHAWKOUTPUTSTATS_H_

namespace Tomahawk{
namespace TWO{

struct TomahawkOutputStatsData{
	typedef TomahawkOutputStatsData self_type;

	TomahawkOutputStatsData() : n_bins(0), bins(nullptr){}
	TomahawkOutputStatsData(const U32 bins) : n_bins(bins), bins(new U64[bins]){
		memset(this->bins, 0, sizeof(U64)*bins);
	}
	~TomahawkOutputStatsData(){ delete [] this->bins; }

	bool setBins(const U32& bins){
		if(bins == 0)
			return false;

		if(this->bins != nullptr){
			delete [] this->bins;
			this->bins = nullptr;
		}

		this->bins = new U64[bins];
		memset(this->bins, 0, sizeof(U64)*bins);

		return true;
	}

	U64& operator[](const U32& p){ return(this->bins[p]); }
	const U64& operator[](const U32& p) const{ return(this->bins[p]); }

	const U64 getTotal(void) const{
		U64 tot = 0;
		for(U32 i = 0; i < this->n_bins; ++i)
			tot += this->bins[i];

		return(tot);
	}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		for(U32 i = 0; i < entry.n_bins; ++i)
			stream << i << '\t' << entry.bins[i] << '\n';

		return(stream);
	}

	U32 n_bins;
	U64* bins;
};

struct TomahawkOutputStats{
	typedef TomahawkOutputStats self_type;
	typedef TomahawkOutputStatsData data_type;

	TomahawkOutputStats() : n_bins(0){}
	TomahawkOutputStats(const U32 bins) : n_bins(bins), within(bins), across(bins), global(bins){}

	bool setBins(const U32& bins){
		if(bins == 0)
			return false;

		if(!this->within.setBins(bins)) return false;
		if(!this->across.setBins(bins)) return false;
		if(!this->global.setBins(bins)) return false;

		return true;
	}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		for(U32 i = 0; i < entry.n_bins; ++i)
			stream << i << '\t' << entry.within[i] << '\t' << entry.across[i] << '\t' << entry.global[i] << '\n';

		return(stream);
	}

	U32 n_bins;
	data_type within;
	data_type across;
	data_type global;
};

struct TomahawkOutputStatsContainer{
	typedef TomahawkOutputStats stats_type;
	typedef IO::OutputEntry entry_type;

	TomahawkOutputStatsContainer(const U32& bins) : n_bins(bins), R2(this->n_bins), D(this->n_bins), Dprime(this->n_bins){}

	void operator+=(const entry_type& entry){
		float R2 = entry.R2;
		if(R2 > 1) R2 = 1;

		float D = entry.D;
		if(entry.D > 1) D = 1;
		if(entry.D < -1) D = -1;

		float DP = entry.Dprime;
		if(DP > 1) DP = 1;
		if(DP < -1) DP = -1;

		++this->R2.global[(U32)(R2 * this->n_bins)];
		++this->Dprime.global[abs((S32)(DP * this->n_bins))];
		++this->D.global[abs((S32)(D * this->n_bins))];
		if(entry.AcontigID == entry.BcontigID){
			++this->R2.within[(U32)(R2 * this->n_bins)];
			++this->Dprime.within[abs((S32)(DP * this->n_bins))];
			++this->D.within[abs((S32)(D * this->n_bins))];
		}
		else{
			++this->R2.across[(U32)(R2 * this->n_bins)];
			++this->Dprime.across[abs((S32)(DP * this->n_bins))];
			++this->D.across[abs((S32)(D * this->n_bins))];
		}
	}

	U32 n_bins;
	stats_type R2, D, Dprime;
};

}
}

#endif /* TOMAHAWKOUTPUT_TOMAHAWKOUTPUTSTATS_H_ */
