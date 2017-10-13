#ifndef ALGORITHM_GENOTYPESAMPLER_H_
#define ALGORITHM_GENOTYPESAMPLER_H_

#include <random>

namespace Tomahawk{
namespace Algorithm{

const char TOMAHAWK_SAMPLER_LOOKUP[4] = {1, 4, 5};

template <class T>
class GenotypeSampler{
	typedef GenotypeSampler self_type;

public:
	explicit GenotypeSampler(const T size) :
		size(size),
		prev_draw_size(0),
		mt(this->rd()),
		genotype_dist(0, 2),
		indices(new T[size]),
		output(new T[size]),
		RLE(nullptr),
		bitvector(nullptr),
		draw_dists(new std::uniform_int_distribution<T>[size])
	{
		// Initiate to 0:n-1
		for(T i = 0; i < this->size; ++i){
			this->indices[i] = i;
			this->draw_dists[i] = std::uniform_int_distribution<T>(0, this->size - i);
		}
	}

	~GenotypeSampler(){
		delete [] this->indices;
		delete [] this->output;
		delete [] this->draw_dists;
	}

	// Fisher-Yates permutation
	bool permute(void){
		for(T i = 0; i < this->size; ++i)
			std::swap(this->indices[this->draw_dists[i](this->mt)], this->indices[i]);

		return true;
	}

	inline void reset(void){ memset(this->output, 0, sizeof(T)*this->size); }
	inline void reset(const T size){
		for(T i = 0; i < size; ++i)
			this->output[this->indices[i]] = 0;
	}

	T* draw(const T size, bool permute = true){
		// Permute data
		if(permute) this->permute();

		// Reset output data vector
		this->reset();

		// Place genotypes at drawn positions
		for(U32 i = 0; i < size; ++i)
			this->output[this->indices[i]] = TOMAHAWK_SAMPLER_LOOKUP[this->genotype_dist(this->mt)];

		// Maintain history
		this->prev_draw_size = size;

		return(this->output);
	}

	friend std::ostream& operator<<(std::ostream& stream, const self_type& self){
		for(U32 i = 0; i < self.size - 1; ++i)
			stream << self.output[i] << ' ';
		stream << self.output[self.size - 1];
		return(stream);
	}

private:
	void encode(void);
	void encodeRLE(void);
	void encodeBitVector(void);

private:
	const T size;
	T prev_draw_size;
	std::random_device rd;
	std::mt19937 mt; // marseinne twister prng
	std::uniform_int_distribution<T> genotype_dist;
	T* indices; //indices
	T* output; //indices
	/*
	 * RLE representation. In the worst case the
	 * length is O(n)
	 */
	char* RLE;
	/*
	 * Bit-vector representation. This has a fixed
	 * memory cost of O(n/register width in b)
	 */
	char* bitvector;
	std::uniform_int_distribution<T>* draw_dists;
};

}
}

#endif /* ALGORITHM_GENOTYPESAMPLER_H_ */
