#ifndef TOMAHAWK_TOMAHAWKSUPPORT_H_
#define TOMAHAWK_TOMAHAWKSUPPORT_H_

namespace Tomahawk{
namespace Support{

#pragma pack(1)
template <class T>
struct TomahawkRun{
public:
	TomahawkRun();	// Disallowed ctor
	~TomahawkRun(); // Disallowed dtor

	T alleleA: Constants::TOMAHAWK_ALLELE_PACK_WIDTH,
	  alleleB: Constants::TOMAHAWK_ALLELE_PACK_WIDTH,
	  runs:    sizeof(T)*8 - Constants::TOMAHAWK_SNP_PACK_WIDTH;
};

#pragma pack(1)
template <class T>
struct TomahawkRunPacked{
public:
	TomahawkRunPacked();	// Disallowed ctor
	~TomahawkRunPacked();	// Disallowed dtor

	T alleles: Constants::TOMAHAWK_SNP_PACK_WIDTH,
	  runs:    sizeof(T)*8 - Constants::TOMAHAWK_SNP_PACK_WIDTH;
};

} // end support

namespace Constants{

const std::string PASS("PASS");
const std::string GT("GT");
const std::string QUAL("100");

}

}



#endif /* TOMAHAWK_TOMAHAWKSUPPORT_H_ */
