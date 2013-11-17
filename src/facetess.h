#ifndef FACETESS3D_H_
#define FACETESS3D_H_

#include "cgal.h"

namespace OpenSCAD {
namespace facetess {

typedef enum tesstype_e {
	NONE,
	CGAL_NEF_STANDARD, // built-in to CGAL's NefPolyhedron->Polyhedron converter.
	CONSTRAINED_DELAUNAY_TRIANGULATION,
	STRAIGHT_SKELETON,
	EARCLIP
} tesstype;

bool is_triangulation( tesstype_e t );

typedef enum tessellater_status_e {
        TESSELLATER_OK,
        BODY_NOT_SIMPLE,
        BODY_NOT_COPLANAR,
        BODY_ONLY_COLLINEAR,
        BODY_LACKS_POINTS,
        HOLE_ONLY_COLLINEAR,
        HOLE_NOT_SIMPLE,
        HOLE_OUTSIDE_BODY,
        HOLE_INSIDE_HOLE,
        HOLE_LACKS_POINTS,
        HOLE_NOT_COPLANAR,
        PROJECTION_2D3D_FAILED,
	CGAL_ERROR,
	CANNOT_EARCLIP_HOLE,
	NONSIMPLE_EAR,
	POLYGON_WITHOUT_EAR,
	UNKNOWN_TESSELLATER,
} tessellater_status;

/* Tessellate the input 3d polygon (with holes) into one or more output
   3d polygons (without holes).

   The first Polygon_3 in the input represents the 'outline' or 'body' contour.
   If there are subsequent Polygon_3s in the input, they represent 'holes'.
   The Body and Holes must all be simple, ordinary non-intersecting polygons.

   The output is a sequence of Polygon_3, none of which will have holes.
*/
tessellater_status tessellate(
	std::vector<CGAL_Polygon_3> &input_pgon3d,
	std::vector<CGAL_Polygon_3> &output_pgons3d,
	tesstype tess );

} // namespace facetess
} // namespace OpenSCAD

#endif
