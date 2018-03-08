#pragma once

#define BEGIN_SFEM_NS namespace sfem {
#define END_SFEM_NS }

#include <cassert>

#include <array>
#include <vector>
#include <string>

////////////////////////////////////////////////////////////////////////////////
BEGIN_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Useful TYPEDEFs

// array type
template<typename T>
using Array = std::vector<T>;



////////////////////////////////////////////////////////////////////////////////
// Global setup

void SetSpaceDim(int nd);
int SpaceDim();

enum CoordType {
	Cartesian,
	RZ,
};

void SetCoordSys(int coord);
int CoordSys();
inline bool IsCartCoord() { return CoordSys() == CoordType::Cartesian; }
inline bool IsRZCoord() { return CoordSys() == CoordType::RZ; }


////////////////////////////////////////////////////////////////////////////////
// 
struct FEEnvironment
{

};




////////////////////////////////////////////////////////////////////////////////
END_SFEM_NS;
////////////////////////////////////////////////////////////////////////////////

