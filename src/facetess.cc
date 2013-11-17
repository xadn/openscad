/*
tessellation of planar faces of 3d CGAL polyhedrons.
tessellation can be with triangles, or other shapes.

Terminology:

Vertex --> a point in space, 2d space or 3d space
Polygon --> a sequence of vertexes
Face --> A list of polygons, the first is a body, the others are holes

Process:

 1. Project input 3d face down into 2d
 2. Convert kernel being used, if needed
 3. Tessellate 2d face into polygons without holes
 4. Convert kernel back to original
 5. Project back up into 3d




Projection Details:

Polygons in 3d space can be 'projected' down into 2d space. A problem is
when they are 'too thin' in one direction and the projection would just
be a line or tiny sliver. To solve this, we project into one of three 2d
planes: xy, yz, or xz. This gives a 2d polygon.

Then we 'rotate' this 2d polygon down so everything is x,y. The rotation
is done by flipping coordinates.

If polygon is 'thin' in z direction, every point is transformed like this:
 3d->2d->3d
  x->x->x
  y->y->y
  z->0->z

If polygon is 'thin' in y direction, every point is transformed like this:
 3d->2d->3d
  x->x->x
  y->0->y
  z->y->z

If polygon is 'thin' in x direction, every point is transformed like this:
 3d->2d->3d
  x->0->x
  y->y->y
  z->x->z


Another feature, to avoid precision problems, is to store a map between
the points in the original kernel and the points in the less-precise kernel
needed to do tessellation.

*/


#include "facetess.h"
#include "printutils.h"
using namespace OpenSCAD::facetess;

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesher_no_edge_refinement_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Dimension.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/connect_holes.h>
#include <CGAL/enum.h>
#include <CGAL/intersections.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Object.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Segment_2_Segment_2_intersection.h>
#include <iostream>

#include <CGAL/Delaunay_mesh_criteria_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesh_area_criteria_2.h>
#include <CGAL/Delaunay_mesh_local_size_criteria_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel CDTKernel;

#include <CGAL/Cartesian.h>
typedef CGAL::Cartesian<CGAL::Gmpq> TessKernel;

typedef TessKernel::FT TessKernel_Number;
typedef TessKernel::Line_2 TessKernel_Line_2;
typedef TessKernel::Point_2 TessKernel_Point_2;
typedef TessKernel::Segment_2 TessKernel_Segment_2;
typedef CGAL::Polygon_2<TessKernel> TessKernel_Polygon_2;
typedef std::vector<TessKernel_Polygon_2> TessKernel_Face_2;
typedef CGAL::Iso_rectangle_2<TessKernel> TessKernel_Iso_rectangle_2;
typedef CGAL::Direction_2<TessKernel> TessKernel_Direction_2;
typedef CGAL::Ray_2<TessKernel> TessKernel_Ray_2;

typedef TessKernel::Point_3 TessKernel_Point_3;
typedef TessKernel::Vector_3 TessKernel_Vector_3;
typedef std::vector<TessKernel_Point_3> TessKernel_Polygon_3;
typedef std::vector<TessKernel_Polygon_3> TessKernel_Face_3;
typedef CGAL::Plane_3<TessKernel> TessKernel_Plane_3;
typedef CGAL::Line_3<TessKernel> TessKernel_Line_3;
typedef CGAL::Direction_3<TessKernel> TessKernel_Direction_3;

typedef CDTKernel::FT CDTKernel_Number;
typedef CDTKernel::Point_2 CDTKernel_Point_2;
typedef CDTKernel::Vector_2 CDTKernel_Vector_2;
typedef CDTKernel::Segment_2 CDTKernel_Segment_2;
typedef CGAL::Polygon_2<CDTKernel> CDTKernel_Polygon_2;
typedef std::vector<CDTKernel_Polygon_2> CDTKernel_Face_2;
typedef CGAL::Iso_rectangle_2<CDTKernel> CDTKernel_Iso_rectangle_2;
typedef CGAL::Triangulation_vertex_base_2<CDTKernel> Vertbase;
typedef CGAL::Delaunay_mesh_face_base_2<CDTKernel> Facebase;
typedef CGAL::Triangulation_data_structure_2<Vertbase, Facebase> Tridatastruct;
typedef CGAL::Constrained_Delaunay_triangulation_2<CDTKernel, Tridatastruct> CDTriangulation;
typedef CDTriangulation::Vertex_handle CDT_Vertex_handle;

typedef CGAL::Delaunay_mesh_criteria_2<CDTriangulation> MeshCriteria;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDTriangulation> sMeshCriteria;
typedef CGAL::Delaunay_mesh_local_size_criteria_2<CDTriangulation> lsMeshCriteria;
typedef CGAL::Delaunay_mesh_area_criteria_2<CDTriangulation> aMeshCriteria;

bool OpenSCAD::facetess::is_triangulation( tesstype_e t ){
        switch (t) {
                case NONE: return true;
                case CGAL_NEF_STANDARD: return true;
                case CONSTRAINED_DELAUNAY_TRIANGULATION: return true;
                case STRAIGHT_SKELETON: return false;
                case EARCLIP: return true;
        }
        return false;
}

template <class T> class DummyMeshCriteria {
public:
        typedef double Quality;
        class Is_bad {
        public:
                CGAL::Mesh_2::Face_badness operator()(const Quality) const {
                        return CGAL::Mesh_2::NOT_BAD;
                }
                CGAL::Mesh_2::Face_badness operator()(const typename T ::Face_handle&, Quality&q) const {
                        q = 1;
                        return CGAL::Mesh_2::NOT_BAD;
                }
        };
        Is_bad is_bad_object() const { return Is_bad(); }
};

CDTKernel_Number sqr( const CDTKernel_Number & x) { return x * x; }
TessKernel_Number sqr( const TessKernel_Number & x) { return x * x; }
CDTKernel_Number blueq( const CDTKernel_Point_2 &a, const CDTKernel_Point_2 &b ) {
	return sqr(a.x()-b.x()) + sqr(a.y()-b.y());
}
TessKernel_Number blueq( const TessKernel_Point_2 &a, const TessKernel_Point_2 &b ) {
	return sqr(a.x()-b.x()) + sqr(a.y()-b.y());
}
CDTKernel_Number blueq( const CDTKernel_Vector_2 &v ) {
	return sqr(v.x())+sqr(v.y());
}
CDTKernel_Number spread( const CDTKernel_Vector_2 &a, const CDTKernel_Vector_2 &b ) {
	CDTKernel_Number dotproduct = a.x()*b.x()+a.y()*b.y();
	return sqr(dotproduct) / ( blueq(a) * blueq( b ) );
}

/* 
from cgal docs

The termination of the algorithm is guaranteed when criteria are 
shape criteria corresponding to a bound on smallest angles not less than 
\( 20.7\) degrees (this corresponds to a radius-edge ratio bound not 
less than \( \sqrt{2}\)). Any size criteria that are satisfied by small 
enough triangle can be added to the set of criteria without compromising 
the termination.

------> note. 20.7 degrees.. a spread of 1/10 will be OK there.

Note that, in the presence of input angles smaller than \( 60\) degrees, 
some bad shaped triangles can appear in the final mesh in the 
neighboring of small angles. To achieve termination and the respect of 
size criteria everywhere, the Is_bad predicate has to return 
CGAL::Mesh_2::IMPERATIVELY_BAD when size criteria are not satisfied, and 
CGAL::Mesh_2::BAD when shape criteria are not satisfied.

MeshingCriteria_2 also provides a type Quality designed to code a 
quality measure for triangles. The type Quality must be less-than 
comparable as the meshing algorithm will order bad triangles by quality, 
to split those with smallest quality first. The predicate Is_bad 
computes the quality of the triangle as a by-product. */

template <class T> class TestMeshCriteria {
public:
        typedef double Quality;
        class Is_bad {
        public:
                CGAL::Mesh_2::Face_badness operator()(const Quality q) const {
			PRINTDB("operator() q: %d", q);
			double epsilon = 0.1;
			if (q>epsilon) return CGAL::Mesh_2::IMPERATIVELY_BAD;
			if (q>0.04) return CGAL::Mesh_2::BAD;
			return CGAL::Mesh_2::BAD;
                }
                CGAL::Mesh_2::Face_badness operator()(const typename T ::Face_handle& fh, Quality&q) const {
			double epsilon = 0.1;
			const CDTriangulation::Point &a = fh->vertex(0)->point();
			const CDTriangulation::Point &b = fh->vertex(1)->point();
			const CDTriangulation::Point &c = fh->vertex(2)->point();
			PRINTDB("operator() criteria %s %s %s",a%b%c);
			CDTKernel_Vector_2 Va( b, c );
			CDTKernel_Vector_2 Vb( c, a );
			CDTKernel_Vector_2 Vc( a, b );
			CDTKernel_Number Qa = blueq( Va );
			CDTKernel_Number Qb = blueq( Vb );
			CDTKernel_Number Qc = blueq( Vc );
			q = Qa*Qb*Qc;
			if (Qa*Qb*Qc>epsilon) return CGAL::Mesh_2::IMPERATIVELY_BAD;
			CDTKernel_Number s1 = spread( Va, Vb );
			CDTKernel_Number s2 = spread( Vb, Vc );
			CDTKernel_Number s3 = spread( Vc, Va );
			CDTKernel_Number minspread = std::min(s1,std::min(s2,s3));
			//CDTKernel_Number arch = archimedes( Qa, Qb, Qc );
			if (minspread<10) return CGAL::Mesh_2::BAD;
//			if ((Qa+Qb+Qc)<0.005) return CGAL::Mesh_2::IMPERATIVELY_BAD;
			//if (arch<0.005) return CGAL::Mesh_2::IMPERATIVELY_BAD;
/*			find min spread
			if minspread < 10 return BAD
			find size using trip quad
			if size < x return IMPERATIVELY BAD
			q = minspread
*/
//			notbadq = q;
			return CGAL::Mesh_2::NOT_BAD;
                }
        };
        Is_bad is_bad_object() const { return Is_bad(); }
};


template <typename SomePoint_2, typename SomePolygon_2>
bool inside( SomePoint_2 &pt, SomePolygon_2 &P )
{
	CGAL::Bounded_side b = CGAL::bounded_side_2( P.vertices_begin(), P.vertices_end(), pt );
	return ( b == CGAL::ON_BOUNDED_SIDE );
}

template <typename SomePoint_2, typename SomePolygon_2>
bool inside( SomePolygon_2 &tocheck, SomePolygon_2 &outside )
{
	bool result = true;
	for (size_t i=0;i<tocheck.size();i++){
		SomePoint_2 pt = tocheck[i];
		if (! inside<SomePoint_2,SomePolygon_2>( pt, outside ) )
			result = false;
	}
	return result;
}

tessellater_status triangulate_cdt( CDTKernel_Face_2 &input, std::vector<CDTKernel_Polygon_2> &output)
{
	PRINTD( "triangulation" );
	// assumes polygons are simple and holes are inside body.
	tessellater_status status = TESSELLATER_OK;
	CGAL::Failure_behaviour old_behaviour = CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
	try {
	if (input.size()==1 && input[0].size()==3) {
		output.push_back(input[0]);
		PRINTD("Already a triangle. Not triangulating");
		return status;
	}
	CDTriangulation cdt;
	PRINTDB( "Polygons in input (including holes): %i",input.size());
	std::vector<CDT_Vertex_handle> vhandles;
	for ( size_t i=0;i<input.size();i++ ) {
		PRINTDB( "inserting cdt, creating vhandle: %i ", i );
		PRINTDB( "points in pgon: %i", input[i].size() );
		vhandles.clear();
		for ( size_t j=0;j<input[i].size();j++ ) {
			PRINTDB( " %s", input[i][j] );;
			vhandles.push_back( cdt.insert( input[i][j] ) );
		}
		for ( size_t k=0;k<vhandles.size();k++ ) {
			int vidx1 = k;
			int vidx2 = (k+1)%vhandles.size();
			cdt.insert_constraint( vhandles[vidx1], vhandles[vidx2] );
		}
	}

	std::vector<CDTriangulation::Point> list_of_seeds;
	for ( size_t   i = 1  ;i<input.size();i++ ) {
		PRINTD( "generating seed for hole");
		CDTKernel_Iso_rectangle_2 bbox = bounding_box( input[i].vertices_begin(), input[i].vertices_end() );
		PRINTDB( "finding bounding box: %s", bbox );
		CGAL::Random_points_in_iso_rectangle_2<CDTKernel_Point_2> generator( bbox.min(), bbox.max() );
		CDTriangulation::Point seedpt = *generator;
		while( !inside<CDTriangulation::Point,CDTKernel_Polygon_2>( seedpt, input[i] ) ) {
			generator++;
			seedpt = *generator;
			PRINTDB( "checking inside: %s", *generator );
		}
		list_of_seeds.push_back( seedpt );
	}

	PRINTDB( "# of seeds: %i", list_of_seeds.size());
	for (size_t i=0;i<list_of_seeds.size();i++) {
		PRINTDB( "seed: %s", list_of_seeds[i] );
	}
	PRINTD( "Meshing..." );
//	CGAL::refine_Delaunay_mesh_2_without_edge_refinement(
	CGAL::refine_Delaunay_mesh_2(
//	  cdt, list_of_seeds.begin(), list_of_seeds.end(), MeshCriteria(0.125,0.5)); // freezes
//	  cdt, list_of_seeds.begin(), list_of_seeds.end(), DummyMeshCriteria<CDTriangulation>()); // freezes
	  cdt, list_of_seeds.begin(), list_of_seeds.end(), sMeshCriteria(0.125,0.5)); // freezes
//	  cdt, list_of_seeds.begin(), list_of_seeds.end(), lsMeshCriteria());
//	  cdt, list_of_seeds.begin(), list_of_seeds.end(), aMeshCriteria()); // freezes
//	  cdt, list_of_seeds.begin(), list_of_seeds.end(), TestMeshCriteria<CDTriangulation>());

	PRINTDB( "Meshed. Number of vertices: %i", cdt.number_of_vertices());
	PRINTDB( "Number of finite faces: %i", cdt.number_of_faces());

	CDTriangulation::Finite_faces_iterator fit;
	for( fit=cdt.finite_faces_begin(); fit!=cdt.finite_faces_end(); ++fit )
	{
		if(fit->is_in_domain()) {
			CDTKernel_Polygon_2 tmp;
			tmp.push_back( cdt.triangle( fit )[0] );
			tmp.push_back( cdt.triangle( fit )[1] );
			tmp.push_back( cdt.triangle( fit )[2] );
			output.push_back( tmp );
		}
	}
	} // try
	catch (const CGAL::Assertion_exception &e) {
		PRINTB("CGAL error in facetess::triangulate %s",e.what());
		status = CGAL_ERROR;
	}
	CGAL::set_error_behaviour(old_behaviour);
	return status;
}

static int fibindex = 0; // fake random earclipping using Fibonacci's sequence
tessellater_status clip_ear( TessKernel_Polygon_2 &pgon, TessKernel_Polygon_2 &newear )
{
	newear.clear();
	fibindex = (1 + fibindex);
	bool earok = false;
	size_t index0, index1, index2 = 0;
	TessKernel_Number maxarea(-1);
	size_t maxarea_index = 0;
	for (size_t i=0;i<pgon.size();i++) {
		newear.clear();
		index0 = (i+fibindex+0)%pgon.size();
		index1 = (i+fibindex+1)%pgon.size();
		index2 = (i+fibindex+2)%pgon.size();
		TessKernel_Point_2 p1 = pgon[index0];
		TessKernel_Point_2 p2 = pgon[index1];
		TessKernel_Point_2 p3 = pgon[index2];
		/*PRINTDB("earsearch. idxs: %i %i %i",index0%index1%index2);
		double ex1 = CGAL::to_double(p1.x());
		double ey1 = CGAL::to_double(p1.y());
		double ex2 = CGAL::to_double(p2.x());
		double ey2 = CGAL::to_double(p2.y());
		double ex3 = CGAL::to_double(p3.x());
		double ey3 = CGAL::to_double(p3.y());
		PRINTB("edata=[[%d,%d],[%d,%d],[%d,%d]]",ex1%ey1%ex2%ey2%ex3%ey3);
		PRINTB("edataRat=[%s %s %s]",p1%p2%p3);
		*/
		newear.push_back( p1 );
		newear.push_back( p2 );
		newear.push_back( p3 );

		if (newear.is_simple()) {
			PRINTD("Ear simple");
			if (newear.orientation()==CGAL::COUNTERCLOCKWISE) {
				PRINTD("Ear CCW");
				bool clean=true;
				for (size_t j=0;j<pgon.size();j++) {
					TessKernel_Point_2 testp = pgon[j];
					if (newear.has_on_bounded_side( testp )) {
						PRINTDB("Ear search fail. pt inside. idx %i, data %s",j%testp);
						clean = false;
						break;
					}
					if (newear.has_on_boundary( testp )) {
						if ( testp != p1 && testp != p2 && testp != p3 ) {
							PRINTDB("Ear search fail, pt on boundary: idx %i data %s",j%testp);
							clean = false;
							break;
						}
					}
				}
				if (clean) earok=true;
			} else {
				PRINTD("Earlcip fail. not CCW");
				PRINTDB("Orientation was: %i",newear.orientation());
			}
		} else {
			PRINTD("Earclip fail. Ear not simple");
		}
		if (earok) {
			if (newear.area()>maxarea) {
				maxarea_index = index1;
			}
			break;
		}
	}
	if (!earok) {
		PRINTD("- pgon had no ear - thats not true! -thats impossible! -");
		for (size_t i=0;i<pgon.size();i++) {
			PRINTDB("noear.pgon.pt %i %s",i%pgon[i]);
		}
		return POLYGON_WITHOUT_EAR;
	} else {
		PRINTDB("chopping ear pt idx %i: %s",index1%pgon[index1]);
		TessKernel_Polygon_2::Vertex_iterator vi = pgon.vertices_begin();
		for (size_t j=0;j<maxarea_index;j++) vi++;
		pgon.erase( vi );
	}
/*	std::stringstream s; s<<"pdata=[";
	for (size_t j=0;j<pgon.size();j++) {
		double x = CGAL::to_double(pgon[j].x());
		double y = CGAL::to_double(pgon[j].y());
		s << " [ " << x << "," << y << "],";
	} s<<"]";*/
	//PRINTB("%s",s.str());
	return TESSELLATER_OK;
}

std::pair<size_t,size_t> find_top_hole_pt( TessKernel_Face_2 &input )
{
	std::pair<size_t,size_t> result(1,0); // indexes: hole, point within hole
	if (input.size()<=1) return result;
	TessKernel_Number highesty = input[result.first][result.second].y();
	for (size_t i=1;i<input.size();i++){
		for (size_t j=0;j<input[i].size();j++){
			const TessKernel_Point_2 p = input[i][j];
			if (p.y()>highesty) {
				highesty = p.y();
				result.first = i;
				result.second = j;
			}
		}
	}
	return result;
}

// assumes body and holes are all simple and non-intersecting
tessellater_status remove_top_hole( TessKernel_Face_2 &input )
{
	if (input.size()<=1) return TESSELLATER_OK;
	std::pair<size_t,size_t> toremove = find_top_hole_pt( input );
	size_t holeidx = toremove.first;
	size_t holept_index = toremove.second;
	TessKernel_Polygon_2 hole_to_rm = input[ holeidx ];
	TessKernel_Point_2 toremove_pt = hole_to_rm[ holept_index ];
	TessKernel_Ray_2 l = TessKernel_Ray_2( toremove_pt, TessKernel_Direction_2(0,1) );
	TessKernel_Polygon_2 body = input[0];
	const TessKernel_Point_2 *ptest;
	TessKernel_Polygon_2 newbody;
	PRINTDB("Remove top hole. pgon idx: %i pt idx: %i ",holeidx%holept_index );
	for (size_t i=0;i<input.size();i++)
		for (size_t j=0;j<input[i].size();j++)
			PRINTDB(" pgon %i pt %i data %s",i%j%input[i][j]);
/*	if (hole_to_rm.orientation()==CGAL::CLOCKWISE) {
		hole_to_rm.reverse_orientation();
		holept_index = hole_to_rm.size() - holept_index;
	}*/
	/// find new pt
	TessKernel_Point_2 newpt;
	size_t cutindex = -1;
	TessKernel_Number lowq(-1);
	for (size_t i=0;i<body.size();i++) {
		TessKernel_Segment_2 s( body[i],body[(i+1)%body.size()] );
		CGAL::Object obj = CGAL::intersection( l, s );
		ptest = CGAL::object_cast<TessKernel_Point_2>(&obj);
		if (ptest) {
			PRINTDB("ray from holept %s, intersects pgon @ %s",toremove_pt%*ptest);
			TessKernel_Number bq = blueq(toremove_pt,*ptest);
			if (bq < lowq || lowq == -1) {
				newpt = *ptest;
				lowq = bq;
				cutindex = i;
			}
		}
	}


	for (size_t i=0;i<body.size();i++) {
		newbody.push_back( body[i] );
		if (i==cutindex) {
			newbody.push_back( newpt );
			for (size_t j=0;j<hole_to_rm.size();j++) {
				size_t index = ( j+holept_index ) % hole_to_rm.size();
				newbody.push_back( hole_to_rm[index] );
			}
			newbody.push_back( hole_to_rm[ holept_index ] );
			newbody.push_back( newpt );
		}
	}

	std::vector<TessKernel_Polygon_2>::iterator vi=input.begin();
	for (size_t i=0;i<holeidx;i++) vi++;
	input.erase( vi );
	input[0] = newbody;
	PRINTD("Removed hole");
	PRINTDB("New Face body. simple: %i #pts: %i",input[0].is_simple()%input[0].size());
		/*std::stringstream s; s<<"posthole=[";
		for (size_t j=0;j<input[0].size();j++) {
			double x = CGAL::to_double(input[0][j].x());
			double y = CGAL::to_double(input[0][j].y());
			s << " [ " << x << "," << y << "],";
		} s<<"]";
		PRINTDB("%s",s.str());*/
	for (size_t i=0;i<input[0].size();i++) {
		PRINTDB(" pt %i: %s",i%input[0][i]);
	}
	PRINTDB("New Face # contours: %i",input.size());
	return TESSELLATER_OK;
}

tessellater_status remove_holes( TessKernel_Face_2 &input )
{
	if (input.size()==0) return TESSELLATER_OK;
	if (input.size()==1) return TESSELLATER_OK;
	size_t numholes = input.size()-1;
	PRINTDB("Num holes: %i",numholes );
	for (size_t h=0;h<numholes;h++) {
		remove_top_hole( input );
		PRINTDB("Num contours: %i",input.size() );
	}
	PRINTDB("Holes removed. #pts in new pgon: %i",input[0].size());
	return TESSELLATER_OK;
}

tessellater_status remove_adjacent_samepoint( TessKernel_Polygon_2 &pgon )
{
	TessKernel_Polygon_2::Vertex_iterator vi = pgon.vertices_begin();
	size_t index0, index1 = 0;
	for (size_t i=0;i<pgon.size();i++) {
		index0 = (i+0)%pgon.size();
		index1 = (i+1)%pgon.size();
		TessKernel_Point_2 p1 = pgon[index0];
		TessKernel_Point_2 p2 = pgon[index1];
		if (p1==p2) {
			pgon.erase(vi);
			break;
		}
		vi++;
	}
	return TESSELLATER_OK;
}

tessellater_status triangulate_earclip(
	TessKernel_Face_2 &input,
	std::vector<TessKernel_Polygon_2> &output)
{
	tessellater_status status = TESSELLATER_OK;
	if (input.size()==1 && input[0].size()==3) {
		output.push_back(input[0]);
		PRINTD("Already a triangle. Not triangulating");
		return status;
	}
	PRINTDB( "triangulate_earclip. input pts: %i. simple: %i", input[0].size()%input[0].is_simple());
	if (input[0].is_simple()) PRINTDB( "triangulate_earclip. orient:%i", input[0].orientation());
	if (input.size()>1) {
		remove_holes( input );
		PRINTDB( "triangulate_earclip after hole removal. input pts: %i. simple: %i", input[0].size()%input[0].is_simple());
		if (input[0].is_simple()) PRINTDB( "triangulate_earclip after hole removal. orient:%i", input[0].orientation());
	}
	size_t size_before=-1;
	while(size_before!=input[0].size()) {
		PRINTDB("Remove adjacent samepoints %i",input[0].size());
		size_before = input[0].size();
		remove_adjacent_samepoint( input[0] );
	}
	TessKernel_Polygon_2 pgon = input[0];
	for (size_t i=0;i<pgon.size();i++) PRINTDB("  pt %i %s ",i%pgon[i]);
	while (pgon.size()>2) {
		TessKernel_Polygon_2 ear;
		status = clip_ear( pgon, ear );
		if (status!=TESSELLATER_OK) return status;
		PRINTDB( " new ear: pts: %i simple: %i",ear.size()%ear.is_simple());
		if (ear.is_simple()) PRINTDB( " new ear orient: %i",ear.orientation());
		PRINTDB( " new pgon: pts: %i simple: %i",pgon.size()%pgon.is_simple());
		if (pgon.is_simple()) PRINTDB( " new pgon orient: %i",pgon.orientation());
		else {
			PRINTD("yikes. new pgon not simple");
		}
		//for (size_t i=0;i<pgon.size();i++) PRINTDB("  pt %i %s ",i%pgon[i]);
		output.push_back( ear );
	}
	//PRINTD("ears all clipped. last pgon: %s",pgon);
	//output.push_back( pgon );
	PRINTDB( "ears in output: %i",output.size());
	return status;
}


/*
void skeletonize( std::vector<TessKernel_Polygon_2> &pgons, std::vector<TessKernel_Polygon_2> &output_pgons2d)
{
	SsBuilder ssb ;
	for (size_t i=0;i<pgons.size();i++) {
		ssb.enter_contour(pgons[i].vertices_begin(),pgons[i].vertices_end());
	}

	boost::shared_ptr<Ss> ss = ssb.construct_skeleton();

	if ( ss ) {
		Ss::Face_const_iterator fi;
		for (fi = ss->faces_begin(); fi!=ss->faces_end(); fi++ ) {
			TessKernel_Polygon_2 tmp;
			Ss::Halfedge_const_handle v1(fi->halfedge());
			Ss::Halfedge_const_handle vi(v1);
			do {
				tmp.push_back( vi->vertex()->point() );
				vi = vi->next();
			} while ( vi != v1 );
			output_pgons2d.push_back( tmp );
		}
	}
	else {
		PRINTD( "construction of skeleton failed" );
	}
}

*/

typedef enum projection_type_e {
	XY_PROJECTION,
	XZ_PROJECTION,
	YZ_PROJECTION
} projection_type;

// find three noncollinear points within the polygon. set p,q,r.
bool find_noncollinear_pts( const TessKernel_Polygon_3 &pgon, TessKernel_Point_3 &p, TessKernel_Point_3 &q, TessKernel_Point_3 &r )
{
	p=pgon[0]; q=p; r=p;
	for ( size_t i=1;i<pgon.size();i++ ) {
		q = pgon[i];
		if ( p!=q ) break;
	}
	for ( size_t i=1;i<pgon.size();i++ ) {
		r = pgon[i];
		if (!CGAL::collinear( p, q, r )) break;
	}
	if (CGAL::collinear( p, q, r )) return false;
	return true;
}


/*


the problem.

1. do CDT in high precision gmpq.

2. CDT freezes on certain inputs

3. the code is 'dropping' some faces from example models

4. are some input triangles being sub-tessellated without needing to be?

*/

bool check_intersect( TessKernel_Segment_2 &s1, TessKernel_Segment_2 &s2 )
{
	TessKernel_Point_2 s1p0 = s1.point(0);
	TessKernel_Point_2 s1p1 = s1.point(1);
	TessKernel_Point_2 s2p0 = s2.point(0);
	TessKernel_Point_2 s2p1 = s2.point(1);
	std::stringstream s;
	s << s1p0;
	std::string cs1p0 = s.str();
	s.str("");
	s << s1p1;
	std::string cs1p1 = s.str();
	s.str("");
	s << s2p0;
	std::string cs2p0 = s.str();
	s.str("");
	s << s2p1;
	std::string cs2p1 = s.str();
	bool cs1 = (cs1p0==cs1p1);
	bool cs2 = (cs2p0==cs2p1);
	//PRINTDB("check intersect. s1[%s,%s] s2[%s,%s]", s1p0%s1p1%s2p0%s2p1);
	//PRINTDB("s1 str comp %s %s ",cs1%cs2);
	if (cs1) return true;
	if (cs2) return true;
	if (CGAL::do_intersect( s1, s2 )) {
		PRINTDB("intersection detected. s1[%s,%s] s2[%s,%s]", s1p0%s1p1%s2p0%s2p1);
		if (cs1p0==cs2p0||cs1p1==cs2p0||cs1p0==cs2p1||cs1p1==cs2p1) {
			if ( (cs1p0==cs2p0&&cs1p1==cs2p1) || (cs1p0==cs2p1&&cs1p1==cs2p0) )
			//if ( (s1p0==s2p0&&s1p1==s2p1) || (s1p0==s2p1&&s1p1==s2p0) )
				return true;
			else
				return false;
		}
		return true;
	}
	return false;
}

// CGAL's is_simple() is buggy. this is very slow but it actually works.
bool simple_check( TessKernel_Polygon_2 &pgon )
{
	if (pgon.size()<3) return false;
	std::vector<TessKernel_Segment_2> segs;
	for ( size_t i=0;i<pgon.size();i++ ) {
		TessKernel_Point_2 p1 = pgon[i];
		TessKernel_Point_2 p2 = pgon[(i+1)%pgon.size()];
		TessKernel_Segment_2 s = TessKernel_Segment_2( p1, p2 );
		segs.push_back( s );
	}
	for ( size_t i=0;i<segs.size();i++ ) {
		TessKernel_Segment_2 s1 = segs[i];
		for ( size_t j=i+2;j<segs.size();j++ ) {
			TessKernel_Segment_2 s2 = segs[j%segs.size()];
			if ( j!=((j+1)%segs.size()) && check_intersect( s1, s2 ) ) {
				return false;
			}
		}
	}
	return true;
}

tessellater_status project_3d_to_2d(
	TessKernel_Face_3 &input, TessKernel_Face_2 &output,
	std::map<TessKernel_Point_2,TessKernel_Point_3> &vertmap,
	TessKernel_Plane_3 &polygon_plane3d,
	projection_type &projection)
{
	tessellater_status status = TESSELLATER_OK;

	// ---------------- project 3d down to 2d

	// imagine the 'normal' vector to the polygon. now imagine shining
	// a light on it from up, then from the back, and then from the side.
	// imagine the 'shadows' of that normal hit the xy plane, yz plane,
	// or xz plane. the "smallest" shadow tells us which 2d plane we want
	// to project our polygon onto, xy, yz, or xz, so that it will be 'fattest'

	PRINTD( "finding three non-collinear points, to find normal");
	if (input[0].size()<3) return BODY_LACKS_POINTS;
	TessKernel_Vector_3 normal;
	TessKernel_Point_3 p, q, r;
	if (!find_noncollinear_pts( input[0], p, q, r ))
		return BODY_ONLY_COLLINEAR;
	normal = CGAL::normal(p,q,r);
	polygon_plane3d = TessKernel_Plane_3( p, q, r );

	TessKernel_Number xyshadow = normal.x()*normal.x()+normal.y()*normal.y();
	TessKernel_Number yzshadow = normal.y()*normal.y()+normal.z()*normal.z();
	TessKernel_Number xzshadow = normal.x()*normal.x()+normal.z()*normal.z();
	TessKernel_Number minshadow = std::min( xyshadow, std::min( yzshadow, xzshadow ) );
	if (minshadow==xyshadow) projection = XY_PROJECTION;
	else if (minshadow==yzshadow) projection = YZ_PROJECTION;
	else projection = XZ_PROJECTION;
	PRINTDB( "shadowssize: xy, yz, xz: %i %i %i", xyshadow % yzshadow % xzshadow );
	PRINTDB( "min shadow: %i", minshadow );

	for (size_t i=0;i<input.size();++i) {
		TessKernel_Polygon_3 tmppoly3d = input[i];
		TessKernel_Polygon_2 tmppoly2d;
		for (size_t j=0;j<tmppoly3d.size();++j) {
			TessKernel_Point_3 tmp3d = tmppoly3d[j];
			TessKernel_Point_2 tmp2d;
			TessKernel_Number x = tmp3d.x();
			TessKernel_Number y = tmp3d.y();
			TessKernel_Number z = tmp3d.z();
			if ( projection == XY_PROJECTION )
				tmp2d = TessKernel_Point_2( x, y);
			else if ( projection == XZ_PROJECTION )
				tmp2d = TessKernel_Point_2( x, z );
			else if ( projection == YZ_PROJECTION )
				tmp2d = TessKernel_Point_2( y, z );
			PRINTDB( "pgon %i pt %i: %s --> %s", i%j%tmp3d % tmp2d );
			vertmap[ tmp2d ] = tmp3d;
			tmppoly2d.push_back( tmp2d );
		}
		output.push_back( tmppoly2d );
	}
	return status;
}

template <typename SomePoint_2,typename SomePolygon_2>
tessellater_status quality_check( std::vector<SomePolygon_2> &pgons )
{
	PRINTD( "quality check");
	if (!pgons[0].is_simple()) return BODY_NOT_SIMPLE;
	//if (!simple_check(pgons[0])) return BODY_NOT_SIMPLE;
	PRINTDB( "body pgon. points: %i", pgons[0].size() );
	PRINTDB( "body is convex: %i", pgons[0].is_convex() );
	PRINTDB( "body is simple: %i", pgons[0].is_simple() );
	int numholes = pgons.size()-1;
	PRINTDB( "holes: %i", numholes );
	for (size_t i=1;i<pgons.size();i++) {
		PRINTDB( "hole pgon: %i points: %i", i%pgons[i].size() );
		PRINTDB( "hole is convex: %i", pgons[i].is_convex() );
		PRINTDB( "hole is simple: %i", pgons[i].is_simple() );
		if (pgons[i].size()<3) return HOLE_LACKS_POINTS;
		if (!pgons[i].is_simple()) return HOLE_NOT_SIMPLE;
		if (!inside<SomePoint_2,SomePolygon_2>(pgons[i], pgons[0])) return HOLE_OUTSIDE_BODY;
		for (size_t j=1;j<pgons.size();j++) {
			if (inside<SomePoint_2,SomePolygon_2>(pgons[i], pgons[j])) return HOLE_INSIDE_HOLE;
		}
	}
	return TESSELLATER_OK;
}

template <typename SomePolygon_2>
tessellater_status fix_orientation( std::vector<SomePolygon_2> &pgons,
	CGAL::Orientation &original_body_orientation )
{
	// must come after simplicity check - it will crash on some non-simple pgons
	PRINTD( "checking orientation, fixing if necessary");
	original_body_orientation = pgons[0].orientation();
	if (pgons[0].orientation() != CGAL::COUNTERCLOCKWISE)
		pgons[0].reverse_orientation();
	for (size_t i=1;i<pgons.size();i++) {
		if (pgons[i].orientation() != CGAL::CLOCKWISE)
			pgons[i].reverse_orientation();
	}
	return TESSELLATER_OK;
}

CGAL::Cartesian_converter<TessKernel,CDTKernel> converter1;
CGAL::Cartesian_converter<CDTKernel,TessKernel> converter2;
void convert_kernel( TessKernel_Face_2 &input,
	CDTKernel_Face_2 &output,
	std::map<CDTKernel_Point_2,TessKernel_Point_2> &vertmap )
{
	for (size_t i=0;i<input.size();++i) {
		CDTKernel_Polygon_2 cdtk_pgon;
		for (size_t j=0;j<input[i].size();++j) {
			TessKernel_Point_2 tk_pt = input[i][j];
			CDTKernel_Point_2 cdtk_pt = converter1( tk_pt );
			vertmap[ cdtk_pt ] = tk_pt;
			cdtk_pgon.push_back( cdtk_pt );
		}
		output.push_back( cdtk_pgon  );
	}
}

void convert_kernel( std::vector<CDTKernel_Polygon_2> &input,
	std::vector<TessKernel_Polygon_2> &output,
	std::map<CDTKernel_Point_2,TessKernel_Point_2> &vertmap )
{
	for (size_t i=0;i<input.size();++i) {
		TessKernel_Polygon_2 tk_pgon;
		for (size_t j=0;j<input[i].size();++j) {
			CDTKernel_Point_2 cdtk_pt = input[i][j];
			TessKernel_Point_2 tk_pt;
			if (vertmap.count(cdtk_pt)) tk_pt = vertmap[ cdtk_pt ];
			else tk_pt = converter2( cdtk_pt );
			tk_pgon.push_back( tk_pt );
		}
		output.push_back( tk_pgon  );
	}
}

tessellater_status cdt_tess( TessKernel_Face_2 &tk_input2dface,
	std::vector<TessKernel_Polygon_2> &tk_output_pgons2d,
	CGAL::Orientation &original_body_orientation )
{
	std::map<CDTKernel_Point_2,TessKernel_Point_2> vertmap_k;
	CDTKernel_Face_2 cdtface;
	std::vector<CDTKernel_Polygon_2> cdtpgons;
 	convert_kernel( tk_input2dface, cdtface, vertmap_k );
	tessellater_status status = quality_check<CDTKernel_Point_2,CDTKernel_Polygon_2>( cdtface );
	if (status!=TESSELLATER_OK) {
		PRINTD("Tessellater - CDT Quality check failed");
		return status;
	}
	fix_orientation<CDTKernel_Polygon_2>( cdtface, original_body_orientation );
	status = triangulate_cdt( cdtface, cdtpgons );
	if (status==TESSELLATER_OK) {
		// fix orentation for output.
		// must do orientation fix w same kernel of output
		for (size_t i=0;i<cdtpgons.size();i++) {
			if (cdtpgons[i].orientation() == original_body_orientation)
				cdtpgons[i].reverse_orientation();
		}
		convert_kernel( cdtpgons, tk_output_pgons2d, vertmap_k );
	}
	return status;
}

tessellater_status earclip_tess( TessKernel_Face_2 &tk_input2dface,
	std::vector<TessKernel_Polygon_2> &tk_output_pgons2d,
	CGAL::Orientation &original_body_orientation )
{
	tessellater_status status = quality_check<TessKernel_Point_2,TessKernel_Polygon_2>( tk_input2dface );
	if (status!=TESSELLATER_OK) {
		PRINTD("Tessellater - TK Quality check failed");
		return status;
	}
	fix_orientation<TessKernel_Polygon_2>( tk_input2dface, original_body_orientation );
	triangulate_earclip( tk_input2dface, tk_output_pgons2d );
	PRINTD( "restore orientation so the 3d normals will be OK");
	// must do orientation w same kernel of output
	for (size_t i=0;i<tk_output_pgons2d.size();i++) {
		if (!tk_output_pgons2d[i].is_simple()) {
			PRINTDB("Error: Non simple ear!: %s",tk_output_pgons2d[i]);
			status=NONSIMPLE_EAR;
			break;
		}
		if (tk_output_pgons2d[i].orientation() == original_body_orientation)
			tk_output_pgons2d[i].reverse_orientation();
	}
	return status;
}

tessellater_status do_tessellation(
	TessKernel_Face_3 &input3dface,
	std::vector<TessKernel_Polygon_3> &output_pgons3d,
	tesstype tess )
{
	tessellater_status status = TESSELLATER_OK;
	TessKernel_Face_2 tk_input2dface;
	std::map<TessKernel_Point_2,TessKernel_Point_3> vertmap;
	TessKernel_Plane_3 polygon_plane3d;
	projection_type projection;
	project_3d_to_2d( input3dface, tk_input2dface, vertmap, polygon_plane3d, projection );

	CGAL::Orientation original_body_orientation = CGAL::CLOCKWISE;

	PRINTD( "tessellate into simple polygons-without-holes");
	std::vector<TessKernel_Polygon_2> tk_output_pgons2d;
	if (tess==CONSTRAINED_DELAUNAY_TRIANGULATION) {
		status = cdt_tess( tk_input2dface, tk_output_pgons2d, original_body_orientation );
	} else if (tess==EARCLIP) {
		status = earclip_tess( tk_input2dface, tk_output_pgons2d, original_body_orientation );
	} else {
		status = UNKNOWN_TESSELLATER;
	}
	if (status!= TESSELLATER_OK) {
		PRINTD("Tessellation failed. no output polygons");
		return status;
	}
	PRINTD( "project from 2d back into 3d");
	for (size_t i=0;i<tk_output_pgons2d.size();i++) {
		TessKernel_Polygon_2 pgon2d = tk_output_pgons2d[i];
		TessKernel_Polygon_3 pgon3d;
		PRINTDB(" Polygon %i. Pts: %i Orient: %s",i%pgon2d.size()%pgon2d.orientation());
		for (size_t j=0;j<pgon2d.size();j++) {
			TessKernel_Point_2 pt2d = pgon2d[j];
			PRINTDB("  2dpt %s",pt2d);
			TessKernel_Point_3 pt3d;
			if (vertmap.count( pt2d )) {
				pt3d = vertmap[ pt2d ];
				PRINTDB("  cached. mapping back to 3d->%s",pt3d);
			} else {
				PRINTDB("  2dpt %s, not cached..calculating 3d...",pt2d);
				TessKernel_Line_3 line3d;
				TessKernel_Number x = pt2d.x();
				TessKernel_Number y = pt2d.y();
				if ( projection == XY_PROJECTION ) {
					TessKernel_Point_3 tmp3d = TessKernel_Point_3( x,y,0 );
					line3d = TessKernel_Line_3( tmp3d, TessKernel_Direction_3(0,0,1) );
				} else if ( projection == XZ_PROJECTION ) {
					TessKernel_Point_3 tmp3d = TessKernel_Point_3( x,0,y );
					line3d = TessKernel_Line_3( tmp3d, TessKernel_Direction_3(0,1,0) );
				} else if ( projection == YZ_PROJECTION ) {
					TessKernel_Point_3 tmp3d = TessKernel_Point_3( 0,x,y );
					line3d = TessKernel_Line_3( tmp3d, TessKernel_Direction_3(1,0,0) );
				}
				CGAL::Object obj = CGAL::intersection( line3d, polygon_plane3d );
				const TessKernel_Point_3 *point_test = CGAL::object_cast<TessKernel_Point_3>(&obj);
				if (point_test) {
					pt3d = *point_test;
					PRINTDB(" 3d point: %s",pt3d);
					vertmap[ pt2d ] = pt3d;
				} else {
					return PROJECTION_2D3D_FAILED;
				}
			} //notcached
			pgon3d.push_back( pt3d );
		}
		output_pgons3d.push_back( pgon3d );
	}

	std::map<int,int> toremove;
	for (size_t i=0;i<output_pgons3d.size();i++){
		std::map<TessKernel_Point_3,int> dups;
		for (size_t j=0;j<output_pgons3d[i].size();j++){
			TessKernel_Point_3 p = output_pgons3d[i][j];
			if (dups.count(p)==0) {
				dups[p]=1;
			} else {
				PRINTD("ERROR: dup 3d pt after tessellation. removing.");
				toremove[i]=1;
			}
		}
	}
	std::vector<TessKernel_Polygon_3> outputs;
	for (size_t i=0;i<output_pgons3d.size();i++){
		if (toremove.count(i)==0) outputs.push_back(output_pgons3d[i]);
	}
	output_pgons3d = outputs;
	return status;
}

// this function is the interface b/t openscad and the tessellater code.
tessellater_status OpenSCAD::facetess::tessellate(
        std::vector<CGAL_Polygon_3> &input_pgon3d,
        std::vector<CGAL_Polygon_3> &output_pgons3d,
        tesstype tess )
{
	tessellater_status status = TESSELLATER_OK;
	PRINTDB( "tessellate. tesstype: %i", tess );
	std::vector<TessKernel_Polygon_3> tess_pgon3d_out;
	status = do_tessellation( input_pgon3d, tess_pgon3d_out, tess );

	PRINTDB( "tessellated. status: %i", status );
	if (status!=TESSELLATER_OK) {
		PRINT("ERROR: Tessellator was not able to process facet. Points:");
		for (size_t i=0;i<input_pgon3d.size();i++) {
			for (size_t j=0;j<input_pgon3d[i].size();j++) {
				TessKernel_Point_3 tess_point = input_pgon3d[i][j];
				PRINTB("ERROR: pgon: %i pt %i coords %s", i%j%tess_point);
			}
		}
		PRINT("ERROR: (end of point list)");
		tess_pgon3d_out = input_pgon3d;
	}

	output_pgons3d.clear();
	PRINTDB("Output polygons: %i",tess_pgon3d_out.size());
	for (size_t i=0;i<tess_pgon3d_out.size();i++) {
		CGAL_Polygon_3 oscad_pgon;
		TessKernel_Polygon_3 tess_pgon = tess_pgon3d_out[i];
		PRINTDB( "Polygon %i, # pts: %i", i%tess_pgon.size() );
		for (size_t j=0;j<tess_pgon.size();j++) {
			TessKernel_Point_3 tess_point = tess_pgon[j];
			CGAL_Point_3 oscad_point;
			oscad_point = tess_point;
			oscad_pgon.push_back( oscad_point );
		}
		output_pgons3d.push_back( oscad_pgon );
	}
	return status;
}
