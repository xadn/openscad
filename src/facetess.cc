/*
 tessellation of planar faces of 3d CGAL polyhedrons.
 tessellation can be with triangles, or other shapes.

 the process is like so:

 1. Change OpenSCAD CGAL kernel points into the kernel needed for tessellation
 2. Project polygon down into 2d
 3. Tessellate in 2d
 4. Project back up into 3d
 5. Change tess kernel points back into OpenSCAD kernel points


notes on projections:

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

*/


#include "facetess.h"
#include "printutils.h"
using namespace OpenSCAD::facetess;

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Dimension.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/partition_2.h>
#include <CGAL/connect_holes.h>
#include <CGAL/enum.h>
//#include <CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/intersections.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Object.h>

#include <CGAL/Segment_2.h>
#include <CGAL/Segment_2_Segment_2_intersection.h>
#include <iostream>


//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//typedef CGAL::Exact_predicates_inexact_constructions_kernel TessKernel;

#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
typedef CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt TessKernel;

//#include <CGAL/Simple_cartesian.h>
//typedef CGAL::Simple_cartesian<CORE::BigRational> TessKernel;

// CGAL delaunay and Skeleton are picky about kernels

typedef TessKernel::Point_2 TK_Point_2;
typedef TessKernel::Segment_2 TK_Segment_2;
typedef TessKernel::Point_3 TK_Point_3;
typedef TessKernel::Vector_3 TK_Vector_3;
//typedef TessKernel::Vector_2 TK_Vector_2;
typedef CGAL::Polygon_2<TessKernel> TK_Polygon_2;
typedef std::vector<TK_Point_3> TK_Polygon_3;
typedef CGAL::Aff_transformation_2<TessKernel> Aff_transformation_2;
//typedef CGAL::Iso_cuboid_3<TessKernel> TK_Iso_cuboid_3;
//typedef CGAL::Polygon_2<TessKernel>::Vertex_iterator TK_Vertex_iterator;
typedef CGAL::Polygon_with_holes_2<TessKernel> Polygon_with_holes_2;
typedef CGAL::Partition_traits_2<TessKernel>::Polygon_2 Partition_Polygon_2;
//typedef Partition_Polygon_2::Vertex_iterator Partition_Vertex_iterator;
typedef CGAL::Plane_3<TessKernel> TK_Plane_3;
typedef CGAL::Line_3<TessKernel> TK_Line_3;

/*
typedef CGAL::Straight_skeleton_2<TessKernel> Ss;
typedef boost::shared_ptr<Ss> SsPtr;
typedef CGAL::Straight_skeleton_builder_traits_2<TessKernel> SsBuilderTraits;
typedef CGAL::Straight_skeleton_builder_2<SsBuilderTraits,Ss> SsBuilder;
*/

typedef CGAL::Direction_3<TessKernel> Direction_3;
typedef CGAL::Triangulation_vertex_base_2<TessKernel> Vertbase;
typedef CGAL::Delaunay_mesh_face_base_2<TessKernel> Facebase;
typedef CGAL::Triangulation_data_structure_2<Vertbase, Facebase> Tridatastruct;
typedef CGAL::Constrained_Delaunay_triangulation_2<TessKernel, Tridatastruct> CDTriangulation;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDTriangulation> MeshCriteria;
typedef CDTriangulation::Vertex_handle CDT_Vertex_handle;
typedef CDTriangulation::Point CDT_Point;
typedef CGAL::Iso_rectangle_2<TessKernel> Iso_rectangle_2;

template <class T> class DummyCriteria {
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

bool inside( TK_Point_2 &pt, TK_Polygon_2 &P )
{
	CGAL::Bounded_side b = CGAL::bounded_side_2( P.vertices_begin(), P.vertices_end(), pt );
	return ( b == CGAL::ON_BOUNDED_SIDE );
}

bool inside( TK_Polygon_2 &tocheck, TK_Polygon_2 &outside )
{
	bool result = true;
	for (size_t i=0;i<tocheck.size();i++){
		TK_Point_2 pt = tocheck[i];
		if (! inside( pt, outside ) ) result = false;
	}
	return result;
}

void triangulate( std::vector<TK_Polygon_2> &pgons, std::vector<TK_Polygon_2> &output_pgons2d)
{
	if (pgons.size()==1 && pgons[0].size()==3) {
		TK_Polygon_2 tmp = pgons[0];
		output_pgons2d.push_back(tmp);
		PRINTD("Already a triangle. Not triangulating");
		return;
	}
	// assumes polygons are simple and holes are inside body.
	CDTriangulation cdt;
	PRINTD( "triangulation" );
	PRINTDB( "polygons in input: %i",pgons.size());

	std::vector<CDT_Vertex_handle> vhandles;
	for ( size_t i=0;i<pgons.size();i++ ) {
		PRINTDB( "inserting cdt, creating vhandle: %i ", i );
		PRINTDB( "points in pgon: %i", pgons[i].size() );
		vhandles.clear();
		for ( size_t j=0;j<pgons[i].size();j++ ) {
			PRINTDB( " %s", pgons[i][j] );;
			vhandles.push_back( cdt.insert( pgons[i][j] ) );
		}
		for ( size_t k=0;k<vhandles.size();k++ ) {
			int vidx1 = k;
			int vidx2 = (k+1)%vhandles.size();
			cdt.insert_constraint( vhandles[vidx1], vhandles[vidx2] );
		}
	}

	std::vector<CDT_Point> list_of_seeds;
	for ( size_t i=1;i<pgons.size();i++ ) {
		PRINTD( "generating seed for hole");
		Iso_rectangle_2 bbox = bounding_box( pgons[i].vertices_begin(), pgons[i].vertices_end() );
		PRINTDB( "finding bounding box: %s", bbox );
		CGAL::Random_points_in_iso_rectangle_2<TK_Point_2> generator( bbox.min(), bbox.max() );
		CDT_Point seedpt = *generator;
		while( !inside( seedpt, pgons[i] ) ) {
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
	CGAL::refine_Delaunay_mesh_2(
//	  cdt, list_of_seeds.begin(), list_of_seeds.end(), MeshCriteria());
	  cdt, list_of_seeds.begin(), list_of_seeds.end(), DummyCriteria<CDTriangulation>());


	PRINTDB( "Meshed. Number of vertices: %i", cdt.number_of_vertices());
	PRINTDB( "Number of finite faces: %i", cdt.number_of_faces());

	//int mesh_faces_counter = 0;
	CDTriangulation::Finite_faces_iterator fit;
	for( fit=cdt.finite_faces_begin(); fit!=cdt.finite_faces_end(); ++fit )
	{
		if(fit->is_in_domain()) {
			TK_Polygon_2 tmp;
			tmp.push_back( cdt.triangle( fit )[0] );
			tmp.push_back( cdt.triangle( fit )[1] );
			tmp.push_back( cdt.triangle( fit )[2] );
			output_pgons2d.push_back( tmp );
		}
	}
}


/*
void skeletonize( std::vector<TK_Polygon_2> &pgons, std::vector<TK_Polygon_2> &output_pgons2d)
{
	SsBuilder ssb ;
	for (size_t i=0;i<pgons.size();i++) {
		ssb.enter_contour(pgons[i].vertices_begin(),pgons[i].vertices_end());
	}

	boost::shared_ptr<Ss> ss = ssb.construct_skeleton();

	if ( ss ) {
		Ss::Face_const_iterator fi;
		for (fi = ss->faces_begin(); fi!=ss->faces_end(); fi++ ) {
			TK_Polygon_2 tmp;
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
bool find_noncollinear_pts( TK_Polygon_3 &pgon, TK_Point_3 &p, TK_Point_3 &q, TK_Point_3 &r )
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

bool check_intersect( TK_Segment_2 &s1, TK_Segment_2 &s2 )
{
	TK_Point_2 s1p0 = s1.point(0);
	TK_Point_2 s1p1 = s1.point(1);
	TK_Point_2 s2p0 = s2.point(0);
	TK_Point_2 s2p1 = s2.point(1);
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
bool simple_check( TK_Polygon_2 &pgon )
{
	if (pgon.size()<3) return false;
	std::vector<TK_Segment_2> segs;
	for ( size_t i=0;i<pgon.size();i++ ) {
		TK_Point_2 p1 = pgon[i];
		TK_Point_2 p2 = pgon[(i+1)%pgon.size()];
		TK_Segment_2 s = TK_Segment_2( p1, p2 );
		segs.push_back( s );
	}
	for ( size_t i=0;i<segs.size();i++ ) {
		TK_Segment_2 s1 = segs[i];
		for ( size_t j=i+2;j<segs.size();j++ ) {
			TK_Segment_2 s2 = segs[j%segs.size()];
			if ( j!=((j+1)%segs.size()) && check_intersect( s1, s2 ) ) {
				return false;
			}
		}
	}
	return true;
}

tessellater_status do_tessellation(
	std::vector<TK_Polygon_3> &input_pgon3d,
	std::vector<TK_Polygon_3> &output_pgons3d,
	tesstype tess )
{
//	return TESSELLATER_OK;

	std::vector<TK_Polygon_3> &pgons = input_pgon3d;

	PRINTDB( "do_tessellation: input 3d pgons: %i", pgons.size());

	// ---------------- project 3d down to 2d

	// imagine the 'normal' vector to the polygon. now imagine shining
	// a light on it from up, then from the back, and then from the side.
	// imagine the 'shadows' of that normal hit the xy plane, yz plane,
	// or xz plane. the "smallest" shadow tells us which 2d plane we want
	// to project our polygon onto, xy, yz, or xz, so that it will be 'fattest'

	PRINTD( "finding three non-collinear points, to find normal");
	if (pgons[0].size()<3) return BODY_LACKS_POINTS;
	TK_Point_3 p, q, r;
	if (!find_noncollinear_pts( pgons[0], p, q, r ))
		return BODY_ONLY_COLLINEAR;
	TK_Plane_3 polygon_plane3d( p, q, r );
	TK_Vector_3 normal;
	projection_type projection;
	normal = CGAL::normal(p,q,r);

	TessKernel::FT xyshadow = normal.x()*normal.x()+normal.y()*normal.y();
	TessKernel::FT yzshadow = normal.y()*normal.y()+normal.z()*normal.z();
	TessKernel::FT xzshadow = normal.x()*normal.x()+normal.z()*normal.z();
	TessKernel::FT minshadow = std::min( xyshadow, std::min( yzshadow, xzshadow ) );
	if (minshadow==xyshadow) projection = XY_PROJECTION;
	else if (minshadow==yzshadow) projection = YZ_PROJECTION;
	else projection = XZ_PROJECTION;
	PRINTDB( "shadowssize: xy, yz, xz: %i %i %i", xyshadow % yzshadow % xzshadow );
	PRINTDB( "min shadow: %i", minshadow );

	std::vector<TK_Polygon_2> input_pgon2d;
	std::map<TK_Point_2,TK_Point_3> vertmap;

	PRINTD( "projecting into to 2d" );
	for (size_t i=0;i<pgons.size();++i) {
		TK_Polygon_3 tmppoly3d = pgons[i];
		TK_Polygon_2 tmppoly2d;
		for (size_t j=0;j<tmppoly3d.size();++j) {
			TK_Point_3 tmp3d = tmppoly3d[j];
			TK_Point_2 tmp2d;
			if ( projection == XY_PROJECTION )
				tmp2d = TK_Point_2( tmp3d.x(), tmp3d.y() );
			else if ( projection == XZ_PROJECTION )
				tmp2d = TK_Point_2( tmp3d.x(), tmp3d.z() );
			else if ( projection == YZ_PROJECTION )
				tmp2d = TK_Point_2( tmp3d.y(), tmp3d.z() );
			PRINTDB( "%s --> %s", tmp3d % tmp2d );
			vertmap[tmp2d] = tmp3d;
			tmppoly2d.push_back( tmp2d );
		}
		input_pgon2d.push_back( tmppoly2d );
	}

	PRINTD( "checking for simplicity, self-intersection, etc");
	PRINTDB( "is convex: %i", input_pgon2d[0].is_convex() );
	PRINTDB( "is simple: %i", input_pgon2d[0].is_simple() );
	PRINTDB( "is simple (oscad): %i", simple_check(input_pgon2d[0]));
	if (!input_pgon2d[0].is_simple()) return BODY_NOT_SIMPLE;
	//if (!simple_check(input_pgon2d[0])) return BODY_NOT_SIMPLE;
	for (size_t i=1;i<input_pgon2d.size();i++) {
		if (input_pgon2d[i].size()<3) return HOLE_LACKS_POINTS;
		if (!input_pgon2d[i].is_simple()) return HOLE_NOT_SIMPLE;
		if (!inside(input_pgon2d[i], input_pgon2d[0])) return HOLE_OUTSIDE_BODY;
		for (size_t j=1;j<input_pgon2d.size();j++) {
			if (inside(input_pgon2d[i], input_pgon2d[j])) return HOLE_INSIDE_HOLE;
		}
	}

	PRINTD( "checking orientation, fixing if necessary");
	// must come after simplicity check - orientation will crash
	CGAL::Orientation original_body_orientation = input_pgon2d[0].orientation();
	if (input_pgon2d[0].orientation() != CGAL::COUNTERCLOCKWISE)
		input_pgon2d[0].reverse_orientation();
	for (size_t i=1;i<input_pgon2d.size();i++) {
		if (input_pgon2d[i].orientation() != CGAL::CLOCKWISE)
			input_pgon2d[i].reverse_orientation();
	}

	PRINTD( "tessellate into simple polygons-without-holes");
	std::vector<TK_Polygon_2> output_pgons2d;
	if (tess==CONSTRAINED_DELAUNAY_TRIANGULATION)
		triangulate( input_pgon2d, output_pgons2d );
//	else
//		skeletonize( input_pgon2d, output_pgons2d );

	PRINTD( "fix orientation so the 3d normals will be OK");
	for (size_t i=0;i<output_pgons2d.size();i++) {
		if (output_pgons2d[i].orientation() == original_body_orientation)
			output_pgons2d[i].reverse_orientation();
	}

	PRINTD( "project from 2d back into 3d");
	for (size_t i=0;i<output_pgons2d.size();i++) {
		TK_Polygon_2 pgon2d = output_pgons2d[i];
		TK_Polygon_3 pgon3d;
		for (size_t j=0;j<pgon2d.size();j++) {
			TK_Point_2 pt2d = pgon2d[j];
			TK_Point_3 pt3d;
			if (vertmap.count( pt2d )) {
				pt3d = vertmap[ pt2d ];
			} else {
				TK_Line_3 line3d;
				if ( projection == XY_PROJECTION ) {
					TK_Point_3 tmp3d = TK_Point_3( pt2d.x(), pt2d.y(), 0 );
					line3d = TK_Line_3( tmp3d, Direction_3(0,0,1) );
				} else if ( projection == XZ_PROJECTION ) {
					TK_Point_3 tmp3d = TK_Point_3( pt2d.x(), 0, pt2d.y() );
					line3d = TK_Line_3( tmp3d, Direction_3(0,1,0) );
				} else if ( projection == YZ_PROJECTION ) {
					TK_Point_3 tmp3d = TK_Point_3( 0, pt2d.x(), pt2d.y() );
					line3d = TK_Line_3( tmp3d, Direction_3(1,0,0) );
				}
				CGAL::Object obj = CGAL::intersection( line3d, polygon_plane3d );
				const TK_Point_3 *point_test = CGAL::object_cast<TK_Point_3>(&obj);
				if (point_test)
					pt3d = *point_test;
				else
					return PROJECTION_2D3D_FAILED;
			}
			pgon3d.push_back( pt3d );
		}
		output_pgons3d.push_back( pgon3d );
	}
	return TESSELLATER_OK;
}

/*
Converters: Convert the number types for the coordinates of a single 3d
point. This is required because CGAL's Delaunay algorithms only can use
certain types of kernel number types, like double or CORE:EXPR. double
is usable, but its not very precise which can cause issues.

Note: If you change the Kernel and Numbertype in cgal.h, you need to
modify these converter functions too.

OpenSCAD has mostly historically used GMPQ as it's number type for CGAL.
GMPQ is a rational type, i.e. the ratio of two 'big integers' (gmpz).

CORE::Expr is a very cool symbolic number type that can do 'as good' as
GMPQ, but is also usable by Delaunay (perhaps because it has square roots)

*/

TK_Point_3 converter1( CGAL_Point_3 p3 )
{
	std::stringstream tmp;

	tmp << p3.x();
	CORE::BigRat xm( tmp.str() );
	tmp.str("");

	tmp << p3.y();
	CORE::BigRat ym( tmp.str() );
	tmp.str("");
	tmp << p3.z();
	CORE::BigRat zm( tmp.str() );
	tmp.str("");

	CORE::Expr xout(xm);
	CORE::Expr yout(ym);
	CORE::Expr zout(zm);

	TK_Point_3 p( xout, yout, zout );
	return p;

/*	TessKernel::FT x( CGAL::to_double(p3.x().numerator())/CGAL::to_double(p3.x().denominator()) );
	TessKernel::FT y( CGAL::to_double(p3.y().numerator())/CGAL::to_double(p3.y().denominator()) );
	TessKernel::FT z( CGAL::to_double(p3.z().numerator())/CGAL::to_double(p3.z().denominator()) );
	TK_Point_3 p( x,y,z );
	return p;
*/
}

CGAL_Point_3 converter2( TK_Point_3 p3 )
{

	std::stringstream tmp;

	CORE::BigRat xm( p3.x().doubleValue() );
	CORE::BigRat ym( p3.y().doubleValue() );
	CORE::BigRat zm( p3.z().doubleValue() );

	tmp << xm;
	CGAL::Gmpq xout(tmp.str());
	tmp.str("");
	tmp << ym;
	CGAL::Gmpq yout(tmp.str());
	tmp.str("");
	tmp << zm;
	CGAL::Gmpq zout(tmp.str());
	tmp.str("");

	CGAL_Point_3 p( xout, yout, zout );
	return p;
}

std::map<TK_Point_3,CGAL_Point_3> kernel_vertmap1;
std::map<CGAL_Point_3,TK_Point_3> kernel_vertmap2;

// this function is the interface b/t openscad and the tessellater code.
tessellater_status OpenSCAD::facetess::tessellate(
        std::vector<CGAL_Polygon_3> &input_pgon3d,
        std::vector<CGAL_Polygon_3> &output_pgons3d,
        tesstype tess )
{
	CGAL::Failure_behaviour old_behaviour = CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
	try {
	PRINTDB( "tessellate. tesstype: %i", tess );
	//CGAL::Cartesian_converter<CGAL_Kernel3,TessKernel> converter1;
	//CGAL::Cartesian_converter<TessKernel,CGAL_Kernel3> converter2;
	std::vector<TK_Polygon_3> tess_pgon3d;
	std::vector<TK_Polygon_3> tess_pgon3d_out;
	PRINTDB( "Input 3d points: %i", input_pgon3d.size());
	PRINTD( "Converting points from CGAL_Kernel3 to TessKernel" );
	for (size_t i=0;i<input_pgon3d.size();i++) {
		CGAL_Polygon_3 oscad_pgon = input_pgon3d[i];
		TK_Polygon_3 tess_pgon;
		for (size_t j=0;j<oscad_pgon.size();j++) {
			CGAL_Point_3 oscad_point = oscad_pgon[j];
			TK_Point_3 tess_point = converter1( oscad_point );
			tess_pgon.push_back( tess_point );
			kernel_vertmap1[tess_point] = oscad_point;
			PRINTDB( " %s -> %s", oscad_point % tess_point );
		}
		tess_pgon3d.push_back( tess_pgon );
	}

	tessellater_status s = do_tessellation( tess_pgon3d, tess_pgon3d_out, tess );
	PRINTDB( "tessellated. status: %i", s );
	if (s!=TESSELLATER_OK) {
		PRINT("WARNING: Tessellator was not able to process facet");
		for (size_t i=0;i<tess_pgon3d.size();i++) {
			for (size_t j=0;j<tess_pgon3d[i].size();j++) {
				TK_Point_3 tess_point = tess_pgon3d[i][j];
				PRINTB("WARNING: %s", tess_point);
			}
		}
		tess_pgon3d_out = tess_pgon3d;
	}

	output_pgons3d.clear();
	PRINTDB("Output polygons: %i",tess_pgon3d_out.size());
	for (size_t i=0;i<tess_pgon3d_out.size();i++) {
		CGAL_Polygon_3 oscad_pgon;
		TK_Polygon_3 tess_pgon = tess_pgon3d_out[i];
		PRINTDB( "pgon %i", i );
		for (size_t j=0;j<tess_pgon.size();j++) {
			TK_Point_3 tess_point = tess_pgon[j];
			CGAL_Point_3 oscad_point;
			PRINTDB( "%i : %s", j%tess_point );
			if (kernel_vertmap1.count(tess_point)) {
				oscad_point = kernel_vertmap1[tess_point];
				PRINTDB( " point was saved in kernel vertex map: %s", oscad_point );
			} else {
				PRINTD( " point not in kernel vertex map. converting" );
				oscad_point = converter2( tess_point );
				PRINTDB( " new point: %s", oscad_point );
			}
			oscad_pgon.push_back( oscad_point );
		}
		output_pgons3d.push_back( oscad_pgon );
	}
	} // try
	catch (const CGAL::Assertion_exception &e) {
		PRINTB("CGAL error in facetess::tessellate %s",e.what());
	}
	CGAL::set_error_behaviour(old_behaviour);
	return TESSELLATER_OK;
}
