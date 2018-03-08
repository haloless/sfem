#pragma once

#include "sfem.h"

BEGIN_SFEM_NS;



enum BasicGeomType
{
	Undefined = -1,

	Point = 0,
	Line,
	Triangle,
	Quadrilateral,

	NumGeom
};


int GeometryDimension(int geom);



END_SFEM_NS;




