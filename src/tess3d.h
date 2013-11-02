#ifndef TESS3D_H_
#define TESS3D_H_

#include "cgal.h"

namespace OpenSCAD {

typedef enum tessellation_e {
	TESS_NONE,
	TESS_CGAL_NEF_STANDARD, // internal to CGAL's Nef->Polyhedron converter
	TESS_CONSTRAINED_DELAUNAY_TRIANGULATION,
	TESS_STRAIGHT_SKELETON
} tessellation;

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
} tessellater_status;

/* Tessellate the input polygon (with holes) into one or more output polygons
   (without holes).

   The first Polygon_3 in the input represents the 'outline' or 'body' contour.
   If there are subsequent Polygon_3s in the input, they represent 'holes'.
   The Body and Holes must all be simple, ordinary non-intersecting polygons.

   The output is a sequence of Polygon_3, none of which will have holes.
*/
tessellater_status tessellate(
        std::vector<CGAL_Polygon_3> &input_pgon3d,
        std::vector<CGAL_Polygon_3> &output_pgons3d,
        tessellation tesstype );

} // namespace OpenSCAD

#endif
