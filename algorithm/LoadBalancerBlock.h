#ifndef ALGORITHM_LOADBALANCERBLOCK_H_
#define ALGORITHM_LOADBALANCERBLOCK_H_

#include <vector>
#include "../support/TypeDefinitions.h"

namespace Tomahawk{

struct LoadBalancerBlock{
	LoadBalancerBlock() : fromRow(0), toRow(0), fromColumn(0), toColumn(0), staggered(false), fromRowAbsolute(0), toRowAbsolute(0), fromColumnAbsolute(0), toColumnAbsolute(0){}
	LoadBalancerBlock(const U32 fromRow, const U32 toRow, const U32 fromColumn, const U32 toColumn) :
		fromRow(fromRow),
		toRow(toRow),
		fromColumn(fromColumn),
		toColumn(toColumn),
		staggered(false),
		fromRowAbsolute(0),
		toRowAbsolute(0),
		fromColumnAbsolute(0),
		toColumnAbsolute(0)
	{}

	LoadBalancerBlock(const U32 fromRow, const U32 toRow, const U32 fromColumn, const U32 toColumn, const U32 fromRowAbsolute, const U32 toRowAbsolute, const U32 fromColumnAbsolute, const U32 toColumnAbsolute, bool stagger = false) :
			fromRow(fromRow),
			toRow(toRow),
			fromColumn(fromColumn),
			toColumn(toColumn),
			staggered(stagger),
			fromRowAbsolute(fromRowAbsolute),
			toRowAbsolute(toRowAbsolute),
			fromColumnAbsolute(fromColumnAbsolute),
			toColumnAbsolute(toColumnAbsolute)
	{}

	~LoadBalancerBlock(){}

	inline U32 getRows(void) const{ return(this->toRow - this->fromRow); }
	inline U32 getColumns(void) const{ return(this->toColumn - this->fromColumn); }
	inline U32 getSize(void) const{
		if(this->isDiagonal())
			return( ((this->getRows() * this->getColumns()) - this->getRows()) / 2 + this->getRows() );
		else
			return(this->getRows() * this->getColumns());
	}

	inline LoadBalancerBlock& operator()(const U32 fromRow, const U32 toRow, const U32 fromColumn, const U32 toColumn, const bool diagonal){
		this->fromRow = fromRow;
		this->toRow = toRow;
		this->fromColumn = fromColumn;
		this->toColumn= toColumn;
		return(*this);
	}

	inline bool isDiagonal(void) const{
		return(this->fromRow == this->fromColumn && this->toRow == this->toColumn);
	}

	inline bool isDiagonalAbsolute(void) const{
		return(this->fromRowAbsolute == this->fromColumnAbsolute && this->toRowAbsolute == this->toColumnAbsolute);
	}

	friend std::ostream& operator<<(std::ostream& os, const LoadBalancerBlock& block){
		os << block.fromRow << '-' << block.toRow << '\t' << block.fromColumn << '-' << block.toColumn << '\t' << block.fromRowAbsolute << '-' << block.toRowAbsolute << '\t' << block.fromColumnAbsolute << '-' << block.toColumnAbsolute << '\t' << (int)block.staggered;
		return(os);
	}

	// Relative order
	U32 fromRow;
	U32 toRow;
	U32 fromColumn;
	U32 toColumn;
	bool staggered;

	// Absolute order
	U32 fromRowAbsolute;
	U32 toRowAbsolute;
	U32 fromColumnAbsolute;
	U32 toColumnAbsolute;
};

}

#endif /* ALGORITHM_LOADBALANCERBLOCK_H_ */
