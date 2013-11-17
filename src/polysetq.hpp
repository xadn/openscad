#ifndef POLYSETQ_H_
#define POLYSETQ_H_

#include "polyset.h"
#include "cgal.h"

/* Auxiliary polyset for dealing with CGAL data structures.
 The vertexes are in a list. A Polygon is a set of indexes into that
 list. A 'volume' is a set of indexes into the polygon list.

 Example: a pyramid
 Vertex list: [0,0,0][1,0,0][0,1,0][0,0,1] (4 items, each with x,y,z coordinate)
 Polygon list: [0,1,2][0,1,3][1,2,3][2,0,3] (4 items, each a list of points)
 Volume list: [0,1,2,3] (1 item, a list of polygons)

 Issues undealt with: invalid shapes, self intersections, normals, etc. Those
 are the responsibility of the caller.
*/
class PolySetQ
{
public:
	std::map<CGAL_HDS::Vertex::Point,size_t> vertmap;
	std::vector<CGAL_HDS::Vertex::Point> vertlist;
	size_t vertcount;
	typedef std::vector<size_t> Polygon;
	std::vector<Polygon> polygons;
	typedef std::vector<size_t> Volume;
	std::vector<Volume> volumes;
	void append_volume() {
		Volume v;
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
