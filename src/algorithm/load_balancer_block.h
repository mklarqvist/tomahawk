#ifndef ALGORITHM_LOAD_BALANCER_BLOCK_H_
#define ALGORITHM_LOAD_BALANCER_BLOCK_H_

#include <vector>
#include "support/type_definitions.h"

namespace tomahawk{

struct LoadBalancerBlock{
	LoadBalancerBlock() :
		fromRow(0),
		toRow(0),
		fromColumn(0),
		toColumn(0),
		staggered(false)
	{}

	LoadBalancerBlock(const U32 fromRow, const U32 toRow, const U32 fromColumn, const U32 toColumn) :
		fromRow(fromRow),
		toRow(toRow),
		fromColumn(fromColumn),
		toColumn(toColumn),
		staggered(false)
	{}

	LoadBalancerBlock(const U32 fromRow, const U32 toRow, const U32 fromColumn, const U32 toColumn, bool staggered) :
		fromRow(fromRow),
		toRow(toRow),
		fromColumn(fromColumn),
		toColumn(toColumn),
		staggered(staggered)
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
		this->fromRow    = fromRow;
		this->toRow      = toRow;
		this->fromColumn = fromColumn;
		this->toColumn   = toColumn;
		return(*this);
	}

	inline bool isDiagonal(void) const{
		return(this->fromRow == this->fromColumn && this->toRow == this->toColumn);
	}

	friend std::ostream& operator<<(std::ostream& os, const LoadBalancerBlock& block){
		os << "[" << block.fromRow << '-' << block.toRow << ", " << block.fromColumn << '-' << block.toColumn << "] staggered: " << (int)block.staggered;
		return(os);
	}

	// Relative order
	U32 fromRow, toRow;
	U32 fromColumn, toColumn;
	bool staggered;
};

struct LoadBalancerThread{
	LoadBalancerThread() :
		row(0),
		fromColumn(0),
		toColumn(0)
	{}

	LoadBalancerThread(const U32 row, const U32 fromColumn, const U32 toColumn) :
		row(row),
		fromColumn(fromColumn),
		toColumn(toColumn)
	{}

	~LoadBalancerThread(){}

	friend std::ostream& operator<<(std::ostream& os, const LoadBalancerThread& block){
		os << block.row << '\t' << block.fromColumn << '-' << block.toColumn;
		return(os);
	}

public:
	// Relative order
	U32 row;
	U32 fromColumn;
	U32 toColumn;
};

}

#endif /* ALGORITHM_LOAD_BALANCER_BLOCK_H_ */
