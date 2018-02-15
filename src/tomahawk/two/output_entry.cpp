#include "output_entry.h"

namespace Tomahawk{
namespace IO{

// These memcpy works because the struct is aligned without padding
OutputEntry::OutputEntry(const char* const data_buffer){
	memcpy(this, data_buffer, sizeof(self_type));
}

// Copy data from stream
OutputEntry::OutputEntry(const IO::BasicBuffer& data_buffer){
	memcpy(this, data_buffer.data(), sizeof(self_type));
}

OutputEntry::OutputEntry(const self_type* const other){
	memcpy(this, other, sizeof(self_type));
}

// Comparator function
// Called from sort helper only
bool OutputEntry::operator<(const self_type& other) const{
	if (this->AcontigID < other.AcontigID) return true;
	if (other.AcontigID < this->AcontigID) return false;

	if (this->Aposition < other.Aposition) return true;
	if (other.Aposition < this->Aposition) return false;

	if (this->BcontigID < other.BcontigID) return true;
	if (other.BcontigID < this->BcontigID) return false;

	if (this->Bposition < other.Bposition) return true;
	if (other.Bposition < this->Bposition) return false;

	return false;
}

bool OutputEntry::operator<=(const self_type& other) const{
	if (this->AcontigID <= other.AcontigID) return true;
	if (other.AcontigID <= this->AcontigID) return false;

	if (this->Aposition <= other.Aposition) return true;
	if (other.Aposition <= this->Aposition) return false;

	if (this->BcontigID <= other.BcontigID) return true;
	if (other.BcontigID <= this->BcontigID) return false;

	if (this->Bposition <= other.Bposition) return true;
	if (other.Bposition <= this->Bposition) return false;

	return true;
}

bool OutputEntry::operator==(const self_type& other) const{
	if (this->AcontigID != other.AcontigID) return false;
	if (this->Aposition != other.Aposition) return false;
	if (this->BcontigID != other.BcontigID) return false;
	if (this->Bposition != other.Bposition) return false;
	return true;
}

void OutputEntry::swapDirection(void){
	U32 Ac = this->AcontigID;
	U32 Bc = this->BcontigID;
	this->AcontigID = Bc;
	this->BcontigID = Ac;
	U32& A = *reinterpret_cast<U32*>(((char*)this + sizeof(U16) +   sizeof(U32)));
	U32& B = *reinterpret_cast<U32*>(((char*)this + sizeof(U16) + 3*sizeof(U32)));
	std::swap(A,B);
}

}
}
