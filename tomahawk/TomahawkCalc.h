#ifndef TOMAHAWK_TOMAHAWKCALC_H_
#define TOMAHAWK_TOMAHAWKCALC_H_

#include "TomahawkReader.h"

namespace Tomahawk {

class TomahawkCalc{
	typedef TomahawkCalc self_type;
	typedef TomahawkCalcParameters parameter_type;
	typedef std::pair<U32,U32> pair_type;
	typedef std::vector<pair_type> pair_vector;
	typedef IO::GenericWriterInterace writer_type;
	typedef Tomahawk::Balancer balancer_type;
	typedef TotempoleReader totempole_reader;
	typedef Interface::ProgressBar progress_type;
	typedef TomahawkReader reader_type;

	// Used to keep track of char pointer offsets in buffer
	// and what totempole entry is associated with that position
	struct DataOffsetPair{
		DataOffsetPair(const char* data, const TotempoleEntry& entry) : entry(entry), data(data){}
		~DataOffsetPair(){}

		const TotempoleEntry& entry;
		const char* data;
	};

public:
	TomahawkCalc();
	~TomahawkCalc();

	bool Open(const std::string input);
	bool Calculate(pair_vector& blocks);
	bool Calculate(std::vector<U32>& blocks);
	bool Calculate();
	bool SelectWriterOutputType(const writer_type::type writer_type);
	void SetOutputType(writer_type::compression type){ this->parameters.compression_type = type; }
	bool OpenWriter(void);
	bool OpenWriter(const std::string destination);
	bool ValidateParameters(void);

private:
	bool CalculateWrapper();
	template <class T> bool Calculate();
	bool WriteTwoHeader(void);
	bool WriteTwoHeaderNatural(void);
	bool WriteTwoHeaderBinary(void);

private:
	U32 threads;
	progress_type progress;
	balancer_type balancer;
	parameter_type parameters;
	reader_type reader;
	writer_type* writer;
};

} /* namespace Tomahawk */

#endif /* TOMAHAWK_TOMAHAWKCALC_H_ */
