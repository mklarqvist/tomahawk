#include <thread>

#include "core.h"

namespace tomahawk {

/****************************
*  Core genotype
****************************/
twk1_t::twk1_t() :
	gt_ptype(0), gt_flipped(0), gt_phase(0), gt_missing(0),
	alleles(0), pos(0), ac(0), an(0), rid(0),
	n_het(0), n_hom(0), gt(nullptr){}
twk1_t::~twk1_t(){ delete gt; }

twk1_t& twk1_t::operator=(const twk1_t& other){
	gt_ptype = other.gt_ptype;
	gt_phase = other.gt_phase;
	gt_missing = other.gt_missing;
	alleles = other.alleles;
	pos = other.pos;
	ac = other.ac;
	an = other.an;
	rid = other.rid;
	n_het = other.n_het;
	n_hom = other.n_hom;
	delete[] gt; gt = nullptr;
	gt = other.gt->Clone();
	return(*this);
}

twk1_t& twk1_t::operator=(twk1_t&& other) noexcept{
	gt_ptype = other.gt_ptype;
	gt_phase = other.gt_phase;
	gt_missing = other.gt_missing;
	alleles = other.alleles;
	pos = other.pos;
	ac = other.ac;
	an = other.an;
	rid = other.rid;
	n_het = other.n_het;
	n_hom = other.n_hom;
	delete[] gt; gt = nullptr;
	other.gt->Move(gt);
	other.gt = nullptr;
	return(*this);
}

void twk1_t::clear(){
	gt_ptype = 0, gt_phase = 0, gt_missing = 0; alleles = 0; pos = 0; ac = 0; an = 0;
	n_het = 0, n_hom = 0;
	if(gt != nullptr) gt->clear();
}

twk_buffer_t& operator<<(twk_buffer_t& buffer, const twk1_t& self){
	assert(self.gt != nullptr);
	uint8_t pack = self.gt_ptype << 3 | self.gt_flipped << 2 | self.gt_phase << 1 | self.gt_missing;
	SerializePrimitive(pack, buffer);
	SerializePrimitive(self.alleles, buffer);
	SerializePrimitive(self.pos, buffer);
	SerializePrimitive(self.ac, buffer);
	SerializePrimitive(self.an, buffer);
	SerializePrimitive(self.rid, buffer);
	SerializePrimitive(self.n_het, buffer);
	SerializePrimitive(self.n_hom, buffer);
	buffer << *self.gt;
	return(buffer);
}

twk_buffer_t& operator>>(twk_buffer_t& buffer, twk1_t& self){
	delete self.gt; self.gt = nullptr;
	uint8_t pack = 0;
	DeserializePrimitive(pack, buffer);
	self.gt_ptype   = pack >> 3;
	self.gt_flipped = (pack >> 2 & 1);
	self.gt_phase   = (pack >> 1) & 1;
	self.gt_missing = pack & 1;
	DeserializePrimitive(self.alleles, buffer);
	DeserializePrimitive(self.pos, buffer);
	DeserializePrimitive(self.ac, buffer);
	DeserializePrimitive(self.an, buffer);
	DeserializePrimitive(self.rid, buffer);
	DeserializePrimitive(self.n_het, buffer);
	DeserializePrimitive(self.n_hom, buffer);
	switch(self.gt_ptype){
	case(1): self.gt = new twk1_igt_t<uint8_t>; break;
	case(2): self.gt = new twk1_igt_t<uint16_t>; break;
	case(4): self.gt = new twk1_igt_t<uint32_t>; break;
	default: std::cerr << "illegal gt primtiive type" << std::endl; exit(1);
	}

	buffer >> *self.gt;

	return(buffer);
}

/****************************
*  Core containers
****************************/
twk1_block_t::twk1_block_t() : n(0), m(0), rid(0), minpos(0), maxpos(0), rcds(nullptr){}
twk1_block_t::twk1_block_t(const uint32_t p): n(0), m(p), rid(0), minpos(0), maxpos(0), rcds(new twk1_t[p]){}
twk1_block_t::~twk1_block_t(){ delete[] rcds; }

twk1_block_t& twk1_block_t::operator=(twk1_block_t&& other){
	//std::cerr << "invoking move ctor" << std::endl;
	delete[] rcds; rcds = nullptr;
	std::swap(rcds, other.rcds);
	n = other.n; m = other.m; rid = other.rid;
	minpos = other.minpos; maxpos = other.maxpos;
	return(*this);
}

void twk1_block_t::Add(const twk1_t& rec){
	if(this->n == this->m) this->resize(); // will also trigger when n=0 and m=0
	this->maxpos = rec.pos + 1; // right non-inclusive
	this->rcds[this->n++] = rec;
}

void twk1_block_t::resize(void){
	if(this->rcds == nullptr){
		this->rcds = new twk1_t[500];
		this->n = 0;
		this->m = 500;
		return;
	}
	twk1_t* temp = rcds;
	rcds = new twk1_t[m*2];
	for(int i = 0; i < n; ++i) rcds[i] = std::move(temp[i]); // move records over
	delete[] temp;
	m*=2;
}

void twk1_block_t::clear(){
	delete[] rcds; rcds = nullptr;
	n = 0; m = 0; rid = 0; minpos = 0; maxpos = 0;
}

twk_buffer_t& operator<<(twk_buffer_t& buffer, const twk1_block_t& self){
	SerializePrimitive(self.n, buffer);
	SerializePrimitive(self.m, buffer);
	SerializePrimitive(self.rid, buffer);
	for(int i = 0; i < self.n; ++i) buffer << self.rcds[i];
	return(buffer);
}

twk_buffer_t& operator>>(twk_buffer_t& buffer, twk1_block_t& self){
	delete[] self.rcds; self.rcds = nullptr;
	DeserializePrimitive(self.n, buffer);
	DeserializePrimitive(self.m, buffer);
	DeserializePrimitive(self.rid, buffer);
	self.rcds = new twk1_t[self.m];
	for(int i = 0; i < self.n; ++i) buffer >> self.rcds[i];
	return(buffer);
}


twk_oblock_t::twk_oblock_t() : n(0), nc(0){}

void twk_oblock_t::Write(std::ostream& stream, const uint32_t n, const uint32_t nc, const twk_buffer_t& buffer){
	uint8_t marker = 1;
	SerializePrimitive(marker, stream);
	SerializePrimitive(n, stream);
	SerializePrimitive(nc, stream);
	stream.write(buffer.data(), buffer.size());
}

std::ostream& operator<<(std::ostream& stream, const twk_oblock_t& self){
	uint8_t marker = 1;
	SerializePrimitive(marker, stream);
	SerializePrimitive(self.n, stream);
	SerializePrimitive(self.nc, stream);
	assert(self.nc == self.bytes.size());
	stream.write(self.bytes.data(), self.bytes.size());
	return(stream);
}

std::istream& operator>>(std::istream& stream, twk_oblock_t& self){
	// marker has to be read outside
	self.bytes.reset();
	DeserializePrimitive(self.n,  stream);
	DeserializePrimitive(self.nc, stream);
	self.bytes.reset(); self.bytes.resize(self.nc);
	stream.read(self.bytes.data(), self.nc);
	self.bytes.n_chars_ = self.nc;
	return(stream);
}

//

twk_ld_settings::twk_ld_settings() :
	square(true), window(false), low_memory(false), bitmaps(false),
	force_phased(false), forced_unphased(false), force_cross_intervals(false),
	c_level(1), bl_size(500), b_size(10000), l_window(1000000),
	n_threads(std::thread::hardware_concurrency()), cycle_threshold(0),
	ldd_load_type(TWK_LDD_ALL), out("-"),
	minP(1), minR2(0.1), maxR2(100), minDprime(0), maxDprime(100),
	n_chunks(1), c_chunk(0)
{}

std::string twk_ld_settings::GetString() const{
	std::string s =  "square=" + std::string((square ? "TRUE" : "FALSE"))
				  + ",window=" + std::string((window ? "TRUE" : "FALSE"))
				  + ",low_memory=" + std::string((low_memory ? "TRUE" : "FALSE"))
				  + ",bitmaps=" + std::string((bitmaps ? "TRUE" : "FALSE"))
				  + ",force_phased=" + std::string((force_phased ? "TRUE" : "FALSE"))
				  + ",force_unphased=" + std::string((forced_unphased ? "TRUE" : "FALSE"))
				  + ",compression_level=" + std::to_string(c_level)
				  + ",block_size=" + std::to_string(bl_size)
				  + ",output_block_size=" + std::to_string(b_size)
				  + (window ? std::string(",window_size=") + std::to_string(l_window) : "")
				  + ",minP=" + std::to_string(minP)
				  + ",minR2=" + std::to_string(minR2)
				  + ",maxR2=" + std::to_string(maxR2)
				  + ",minDprime=" + std::to_string(minDprime)
				  + ",maxDprime=" + std::to_string(maxDprime)
				  + ",n_chunks=" + std::to_string(n_chunks)
				  + ",c_chunk=" + std::to_string(c_chunk)
				  + ",n_threads=" + std::to_string(n_threads)
				  + ",ldd_type=" + std::to_string((int)ldd_load_type)
				  + ",cycle_threshold=" + std::to_string(cycle_threshold);
	return(s);
}

//
twk_igt_vec::twk_igt_vec():
	n(0),
	front_zero(0),
	tail_zero(0),
	data(nullptr),
	mask(nullptr)
{
}

twk_igt_vec::~twk_igt_vec(){
	aligned_free(data);
	aligned_free(mask);
}

bool twk_igt_vec::Build(const twk1_t& rec,
                        const uint32_t& n_samples)
{
	if(n == 0){
		n = ceil((double)(n_samples*2)/64);
		n += (n*64) % 128; // must be divisible by 128-bit register
		data = reinterpret_cast<uint64_t*>(aligned_malloc(n*sizeof(uint64_t), SIMD_ALIGNMENT));
		if(rec.gt_missing){
			mask = reinterpret_cast<uint64_t*>(aligned_malloc(n*sizeof(uint64_t), SIMD_ALIGNMENT));
		}
	}

	memset(data, 0, n*sizeof(uint64_t));
	if(rec.gt_missing)
		memset(mask, 0, n*sizeof(uint64_t));

	uint32_t cumpos = 0;
	for(int i = 0; i < rec.gt->n; ++i){
		const uint32_t len  = rec.gt->GetLength(i);
		const uint8_t  refA = rec.gt->GetRefA(i);
		const uint8_t  refB = rec.gt->GetRefB(i);

		if(refA == 0 && refB == 0){
			cumpos += 2*len;
			continue;
		}

		for(int j = 0; j < 2*len; j+=2){
			if(refA == 1){ this->SetData(cumpos + j + 0); }
			if(refB == 1){ this->SetData(cumpos + j + 1); }
			if(refA == 2){ this->SetMask(cumpos + j + 0); this->SetMask(cumpos + j + 1); }
			if(refB == 2){ this->SetMask(cumpos + j + 0); this->SetMask(cumpos + j + 1); }
		}
		cumpos += 2*len;
	}

	if(rec.gt_missing){
		// Invert mask
		/*for(int i = 0; i < n; ++i){
			mask[i] = ~mask[i];
		}*/
	}
	assert(cumpos == n_samples*2);


	const uint32_t byteAlignedEnd  = (2*n_samples/SIMD_WIDTH) * (SIMD_WIDTH / 64);

	int j = 0;
	if(rec.gt_missing){
		// Search from left->right
		for(; j < byteAlignedEnd; ++j){
			if(this->data[j] != 0 || this->mask[j] != 0)
				break;
		}

		// Front of zeroes
		this->front_zero = ((j - 1 < 0 ? 0 : j - 1)*4)/GENOTYPE_TRIP_COUNT;
		if(j != byteAlignedEnd){
			j = byteAlignedEnd - 1;
			for(; j > 0; --j){
				if(this->data[j] != 0 || this->mask[j] != 0)
					break;
			}
		}
	}
	// If no data is missing then we have no mask. Run this subroutine instead
	// to avoid memory problems.
	else {
		// Search from left->right
		for(; j < byteAlignedEnd; ++j){
			if(this->data[j] != 0)
				break;
		}

		// Front of zeroes
		this->front_zero = ((j - 1 < 0 ? 0 : j - 1)*4)/GENOTYPE_TRIP_COUNT;
		if(j != byteAlignedEnd){
			j = byteAlignedEnd - 1;
			for(; j > 0; --j){
				if(this->data[j] != 0)
					break;
			}
		}
	}

	// Tail of zeroes
	this->tail_zero = byteAlignedEnd == j ? 0 : ((byteAlignedEnd - (j+1))*4)/GENOTYPE_TRIP_COUNT;

	return(true);
}

// twk1_two
twk1_two_t::twk1_two_t() :
	controller(0), ridA(0), ridB(0), Amiss(0), Aphased(0), Apos(0), Bmiss(0),
	Bphased(0), Bpos(0), R(0), R2(0), D(0), Dprime(0), P(0),
	ChiSqModel(0), ChiSqFisher(0)
{
	memset(cnt, 0, sizeof(double)*4);
}

void twk1_two_t::clear(){
	controller = 0; ridA = 0; ridB = 0;
	Amiss = 0; Aphased = 0; Apos = 0;
	Bmiss = 0; Bphased = 0; Bpos = 0;
	R = 0; R2 = 0; D = 0; Dprime = 0; P = 0;
	ChiSqModel = 0; ChiSqFisher = 0;
	memset(cnt, 0, sizeof(double)*4);
}

bool twk1_two_t::operator<(const twk1_two_t& other) const{
	if(ridA < other.ridA) return true;
	if(other.ridA < ridA) return false;
	if(ridB < other.ridB) return true;
	if(other.ridB < ridB) return false;
	if(Apos < other.Apos) return true;
	if(other.Apos < Apos) return false;
	if(Bpos < other.Bpos) return true;
	if(other.Bpos < Bpos) return false;
	return false;
}

twk_buffer_t& operator<<(twk_buffer_t& os, const twk1_two_t& entry){
	SerializePrimitive(entry.controller, os);
	SerializePrimitive(entry.ridA, os);
	SerializePrimitive(entry.ridB, os);
	uint32_t packA = entry.Apos << 2 | entry.Aphased << 1 | entry.Amiss;
	uint32_t packB = entry.Bpos << 2 | entry.Bphased << 1 | entry.Bmiss;
	SerializePrimitive(packA, os);
	SerializePrimitive(packB, os);
	SerializePrimitive(entry.cnt[0], os);
	SerializePrimitive(entry.cnt[1], os);
	SerializePrimitive(entry.cnt[2], os);
	SerializePrimitive(entry.cnt[3], os);
	SerializePrimitive(entry.D, os);
	SerializePrimitive(entry.Dprime, os);
	SerializePrimitive(entry.R, os);
	SerializePrimitive(entry.R2, os);
	SerializePrimitive(entry.P, os);
	SerializePrimitive(entry.ChiSqFisher, os);
	SerializePrimitive(entry.ChiSqModel, os);
	return os;
}

twk_buffer_t& operator>>(twk_buffer_t& os, twk1_two_t& entry){
	DeserializePrimitive(entry.controller, os);
	DeserializePrimitive(entry.ridA, os);
	DeserializePrimitive(entry.ridB, os);
	uint32_t packA = 0;
	uint32_t packB = 0;
	DeserializePrimitive(packA, os);
	DeserializePrimitive(packB, os);
	entry.Apos = packA >> 2;
	entry.Bpos = packB >> 2;
	entry.Aphased = (packA >> 1) & 1;
	entry.Bphased = (packB >> 1) & 1;
	entry.Amiss = packA & 1;
	entry.Bmiss = packB & 1;
	DeserializePrimitive(entry.cnt[0], os);
	DeserializePrimitive(entry.cnt[1], os);
	DeserializePrimitive(entry.cnt[2], os);
	DeserializePrimitive(entry.cnt[3], os);
	DeserializePrimitive(entry.D, os);
	DeserializePrimitive(entry.Dprime, os);
	DeserializePrimitive(entry.R, os);
	DeserializePrimitive(entry.R2, os);
	DeserializePrimitive(entry.P, os);
	DeserializePrimitive(entry.ChiSqFisher, os);
	DeserializePrimitive(entry.ChiSqModel, os);
	return(os);
}

std::ostream& twk1_two_t::PrintLD(std::ostream& os) const {
	os << controller << '\t' << ridA << '\t' << Apos << '\t' << ridB << '\t' << Bpos << '\t'
	   << cnt[0] << '\t' << cnt[1] << '\t' << cnt[2] << '\t' << cnt[3] << '\t' << D << '\t'
	   << Dprime << '\t' << R << '\t' << R2 << '\t' << P << '\t' << ChiSqFisher << '\t' << ChiSqModel << '\n';
	return(os);
}

std::ostream& twk1_two_t::PrintLDJson(std::ostream& os) const {
	os << '[' << controller << ',' << ridA << ',' << Apos << ',' << ridB << ',' << Bpos << ','
	   << cnt[0] << ',' << cnt[1] << ',' << cnt[2] << ',' << cnt[3] << ',' << D << ','
	   << Dprime << ',' << R << ',' << R2 << ',' << P << ',' << ChiSqFisher << ',' << ChiSqModel << ']' << '\n';
	return(os);
}

std::ostream& operator<<(std::ostream& os, const twk1_two_t& entry){
	os << entry.controller << '\t' << entry.ridA << '\t' << entry.Apos << '\t' << entry.ridB << '\t' << entry.Bpos << '\t'
	   << entry.cnt[0] << '\t' << entry.cnt[1] << '\t' << entry.cnt[2] << '\t' << entry.cnt[3] << '\t'
	   << entry.D << '\t' << entry.Dprime << '\t' << entry.R << '\t' << entry.R2 << '\t' << entry.P << '\t' << entry.ChiSqFisher << '\t' << entry.ChiSqModel;
	return os;
}

//
twk_oblock_two_t::twk_oblock_two_t() : n(0), nc(0){}

void twk_oblock_two_t::Write(std::ostream& stream, const uint32_t n, const uint32_t nc, const twk_buffer_t& buffer){
	uint8_t marker = 1;
	SerializePrimitive(marker, stream);
	SerializePrimitive(n, stream);
	SerializePrimitive(nc, stream);
	stream.write(buffer.data(), buffer.size());
}

std::ostream& operator<<(std::ostream& stream, const twk_oblock_two_t& self){
	uint8_t marker = 1;
	SerializePrimitive(marker, stream);
	SerializePrimitive(self.n, stream);
	SerializePrimitive(self.nc, stream);
	assert(self.nc == self.bytes.size());
	stream.write(self.bytes.data(), self.bytes.size());
	return(stream);
}

std::istream& operator>>(std::istream& stream, twk_oblock_two_t& self){
	// marker has to be read outside
	self.bytes.reset();
	DeserializePrimitive(self.n,  stream);
	DeserializePrimitive(self.nc, stream);
	self.bytes.reset(); self.bytes.resize(self.nc);
	stream.read(self.bytes.data(), self.nc);
	self.bytes.n_chars_ = self.nc;
	return(stream);
}

//
twk1_two_block_t::twk1_two_block_t() : n(0), m(0), rcds(nullptr){}
twk1_two_block_t::twk1_two_block_t(const uint32_t p): n(0), m(p), rcds(new twk1_two_t[p]){}
twk1_two_block_t::~twk1_two_block_t(){ delete[] rcds; }

twk1_two_block_t& twk1_two_block_t::Add(const twk1_two_t& rec){
	if(n == m) this->resize(); // will also trigger when n=0 and m=0
	rcds[n++] = rec;
	return(*this);
}

void twk1_two_block_t::resize(const uint32_t p){
	if(p < n) return;

	//std::cerr << "resizing to " << p << " @ " << n << "/" << m << std::endl;

	twk1_two_t* temp = rcds;
	rcds = new twk1_two_t[p];
	for(int i = 0; i < n; ++i) rcds[i] = std::move(temp[i]); // move records over
	delete[] temp;
	m = p;
}

void twk1_two_block_t::resize(void){
	if(this->rcds == nullptr){
		this->rcds = new twk1_two_t[500];
		this->n = 0;
		this->m = 500;
		return;
	}
	//std::cerr << "resizing=" << n << "/" << m << std::endl;

	twk1_two_t* temp = rcds;
	rcds = new twk1_two_t[m*2];
	for(int i = 0; i < n; ++i) rcds[i] = std::move(temp[i]); // move records over
	delete[] temp;
	m*=2;
}

void twk1_two_block_t::reset(){
	n = 0;
}

void twk1_two_block_t::clear(){
	delete[] rcds; rcds = nullptr;
	n = 0; m = 0;
}

bool twk1_two_block_t::Sort(){
	//std::cerr << "sorting=" << n << std::endl;
	std::sort(start(), end());
	return(true);
}

twk_buffer_t& operator<<(twk_buffer_t& buffer, const twk1_two_block_t& self){
	SerializePrimitive(self.n, buffer);
	SerializePrimitive(self.m, buffer);
	for(int i = 0; i < self.n; ++i) buffer << self.rcds[i];
	return(buffer);
}

twk_buffer_t& operator>>(twk_buffer_t& buffer, twk1_two_block_t& self){
	delete[] self.rcds; self.rcds = nullptr;
	DeserializePrimitive(self.n, buffer);
	DeserializePrimitive(self.m, buffer);
	self.rcds = new twk1_two_t[self.m];
	for(int i = 0; i < self.n; ++i) buffer >> self.rcds[i];
	return(buffer);
}

}
