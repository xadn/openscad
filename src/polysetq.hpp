#ifndef POLYSETQ_H_
#define POLYSETQ_H_

#include "polyset.h"
#include "cgal.h"
#include "linalg.h"

/* PolySetQ is an auxiliary to PolySet for dealing with CGAL data
structures, especially during import/export. It is a bit similar to
PolySet, but it does not form the backbone of OpenSCAD like PolySet does.

Veretxes are x,y,z coordinates in space. A polygon is a sequence of vertexes.
A volume is a sequence of Polygons.

Polygons dont store vertexes, but instead indexes to the Vertex list.
Volumes dont store polygons, but instead indexes to the Polygon list.

Polygons are assumed to be 'flat' aka 'planar' and without holes.
Volumes are assumed to represent collections of polygons that form
valid, closed 3d shape such that every 'inside' the shape is separated
from 'empty space' by a polygon. This is also called a '2-manifold' or
'Polyhedron'. Volumes are assumed not to intersect. Polygons are assumed
to be 'counterclockwise' in point order, if viewed from outside the
shape.

Problems like self intersections, invalid shapes, badly oriented
polygons, etc, are the responsibility of the caller.

Example: a triangular thingy
 Code:
 PolySetQ q;
 q.append_volume();
 q.append_polygon();
 q.append_vertex(Point(0,0,0));
 q.append_vertex(Point(1,0,0));
 q.append_vertex(Point(0,1,0));
 q.append_polygon();
 q.append_vertex(Point(0,0,0));
 q.append_vertex(Point(1,0,0));
 q.append_vertex(Point(0,0,1));
 q.append_polygon();
 q.append_vertex(Point(1,0,0));
 q.append_vertex(Point(0,1,0));
 q.append_vertex(Point(0,0,1));
 q.append_polygon();
 q.append_vertex(Point(0,1,0));
 q.append_vertex(Point(0,0,0));
 q.append_vertex(Point(0,0,1));

 Resulting data structure:
 Vertex list: [0,0,0][1,0,0][0,1,0][0,0,1] (4 items, each with x,y,z coordinate)
 Polygon list: [0,1,2][0,1,3][1,2,3][2,0,3] (4 items, each a list of point indexes)
 Volume list: [0,1,2,3] (1 item, a list of polygon indexes)

*/
class VolumeQ : public std::vector<size_t>
{
public:
	Color4f color;
	VolumeQ( Color4f c ) : color(c) { }
	VolumeQ() { this->color=Color4f(128,200,30); }
};
class PolySetQ
{
public:
	std::map<CGAL_HDS::Vertex::Point,size_t> vertmap;
	std::vector<CGAL_HDS::Vertex::Point> vertlist;
	size_t vertcount;
	typedef std::vector<size_t> Polygon;
	std::vector<Polygon> polygons;
	std::vector< VolumeQ > volumes;
	void append_volume() {
		VolumeQ v;
		volumes.push_back(v);
	}
	void append_poly() {
		Polygon p;
		polygons.push_back(p);
		size_t current_volume = volumes.size()-1;
		size_t current_polygon = polygons.size()-1;
		volumes[current_volume].push_back( current_polygon );
	}
	void append_vertex(CGAL_HDS::Vertex::Point &p) {
		if (vertmap.count(p)==0) {
			vertmap[p] = vertcount;
			vertlist.push_back(p);
			vertcount++;
		}
		size_t current_polygon = polygons.size()-1;
		polygons[current_polygon].push_back( vertmap[p] );
	}
	void merge_polygons( PolySetQ &Q ) {
		// merge polygons from Q into current volume, if any
		for (size_t i=0;i<Q.polygons.size();i++) {
			this->append_poly();
			for (size_t j=0;j<Q.polygons[i].size();j++) {
				CGAL_Point_3 p = Q.vertlist[Q.polygons[i][j]];
				this->append_vertex( p );
			}
		}
	}
	void clear() {
		vertmap.clear();
		vertlist.clear();
		polygons.clear();
		volumes.clear();
		vertcount = 0;
	}
	std::string dump() {
		std::stringstream s;
		for (size_t i=0;i<vertlist.size();i++){
			s << "vertex " << i << ": " << vertlist[i] << "\n";
		}
		for (size_t i=0;i<polygons.size();i++){
			s << "pgon " << i << ":";
			for (size_t j=0;j<polygons[i].size();j++){
				s << " " << polygons[i][j];
			} s << "\n";
		}
		for (size_t i=0;i<volumes.size();i++){
			s << "volume " << i << ":";
			for (size_t j=0;j<volumes[i].size();j++){
				s << " " << volumes[i][j];
			} s << "\n";
		}
		return s.str();
	}
};

#endif
