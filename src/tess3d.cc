#include "tess3d.h"

using namespace OpenSCAD;

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
#include <CGAL/Straight_skeleton_builder_2.h>
#include <CGAL/intersections.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Object.h>
#include <iostream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel TessKernel;
//apparently delaunay triangle meshing requires a kernel with sqrt?

typedef TessKernel::Point_2 Point_2;
typedef TessKernel::Point_3 Point_3;
typedef TessKernel::Vector_3 Vector_3;
typedef TessKernel::Vector_2 Vector_2;
typedef CGAL::Aff_transformation_2<TessKernel> Aff_transformation_2;
typedef CGAL::Iso_cuboid_3<TessKernel> Iso_cuboid_3;
typedef CGAL::Polygon_2<TessKernel> Polygon_2;
typedef std::vector<Point_3> Polygon_3;
typedef CGAL::Polygon_2<TessKernel>::Vertex_iterator Vertex_iterator;
typedef CGAL::Polygon_with_holes_2<TessKernel> Polygon_with_holes_2;
typedef CGAL::Partition_traits_2<TessKernel>::Polygon_2 Partition_Polygon_2;
typedef Partition_Polygon_2::Vertex_iterator Partition_Vertex_iterator;
typedef CGAL::Plane_3<TessKernel> Plane_3;
typedef CGAL::Line_3<TessKernel> Line_3;
typedef CGAL::Straight_skeleton_2<TessKernel> Ss;
typedef boost::shared_ptr<Ss> SsPtr;
typedef CGAL::Straight_skeleton_builder_traits_2<TessKernel> SsBuilderTraits;
typedef CGAL::Straight_skeleton_builder_2<SsBuilderTraits,Ss> SsBuilder;
typedef CGAL::Direction_3<TessKernel> Direction_3;

typedef CGAL::Triangulation_vertex_base_2<TessKernel> Vertbase;
typedef CGAL::Delaunay_mesh_face_base_2<TessKernel> Facebase;
typedef CGAL::Triangulation_data_structure_2<Vertbase, Facebase> Tridatastruct;
typedef CGAL::Constrained_Delaunay_triangulation_2<TessKernel, Tridatastruct> CDTriangulation;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDTriangulation> MeshCriteria;
typedef CDTriangulation::Vertex_handle CDT_Vertex_handle;
typedef CDTriangulation::Point CDT_Point;
typedef CGAL::Iso_rectangle_2<TessKernel> Iso_rectangle_2;

/* small note on projection

1. avoid matrices to keep things simple
2. need to deal with xy, yz, and xz planar polygons
3. so we just flip. first we flip from 3d->2d, then 2d->3d.

  3d->2d->3d

like so:

If polygon is 'thin' in z direction:
  x->x->x
  y->y->y
  z->0->z

If polygon is 'thin' in y direction:
  x->x->x
  y->0->y
  z->y->z

If polygon is 'thin' in x direction:
  x->0->x
  y->y->y
  z->x->z

*/

bool inside( Point_2 &pt, Polygon_2 &P )
{
	CGAL::Bounded_side b = CGAL::bounded_side_2( P.vertices_begin(), P.vertices_end(), pt );
	return ( b == CGAL::ON_BOUNDED_SIDE );
}

bool inside( Polygon_2 &tocheck, Polygon_2 &outside )
{
	bool result = true;
	for (int i=0;i<tocheck.size();i++){
		Point_2 pt = tocheck[i];
		if (! inside( pt, outside ) ) result = false;
	}
	return result;
}

void triangulate( std::vector<Polygon_2> &pgons, std::vector<Polygon_2> &output_pgons2d)
{
	std::stringstream debug;
	// assumes polygons are simple and holes are inside body.
	CDTriangulation cdt;
	debug << "triangulation\n";

	std::vector<CDT_Vertex_handle> vhandles;
	for ( int i=0;i<pgons.size();i++ ) {
		debug << "inserting cdt, creating vhandle: " << i << "\n";
		vhandles.clear();
		for ( int j=0;j<pgons[i].size();j++ ) {
			vhandles.push_back( cdt.insert( pgons[i][j] ) );
		}
		for ( int k=0;k<vhandles.size();k++ ) {
			int vidx1 = k;
			int vidx2 = (k+1)%vhandles.size();
			cdt.insert_constraint( vhandles[vidx1], vhandles[vidx2] );
		}
	}

	std::vector<CDT_Point> list_of_seeds;
	for ( int i=1;i<pgons.size();i++ ) {
		debug << "generating seed for hole\n";
		Iso_rectangle_2 bbox = bounding_box( pgons[i].vertices_begin(), pgons[i].vertices_end() );
		CGAL::Random_points_in_iso_rectangle_2<Point_2> generator( bbox.min(), bbox.max() );
		CDT_Point seedpt = *generator;
		while( !inside( seedpt, pgons[i] ) ) {
			generator++;
			seedpt = *generator;
			debug << "checking inside:" << *generator << "\n";
		}
		list_of_seeds.push_back( seedpt );
	}

	debug << "# of seeds: " << list_of_seeds.size() << std::endl;

	debug << "Meshing..." << std::endl;
	CGAL::refine_Delaunay_mesh_2(
	  cdt, list_of_seeds.begin(), list_of_seeds.end(), MeshCriteria());


	debug << "Meshed. Number of vertices: " << cdt.number_of_vertices() << std::endl;
	debug << "Number of finite faces: " << cdt.number_of_faces() << std::endl;

	int mesh_faces_counter = 0;
	CDTriangulation::Finite_faces_iterator fit;
	for( fit=cdt.finite_faces_begin(); fit!=cdt.finite_faces_end(); ++fit )
	{
		if(fit->is_in_domain()) {
			Polygon_2 tmp;
			tmp.push_back( cdt.triangle( fit )[0] );
			tmp.push_back( cdt.triangle( fit )[1] );
			tmp.push_back( cdt.triangle( fit )[2] );
			output_pgons2d.push_back( tmp );
		}
	}
}


void skeletonize( std::vector<Polygon_2> &pgons, std::vector<Polygon_2> &output_pgons2d)
{
	std::stringstream debug;
	SsBuilder ssb ;
	for (int i=0;i<pgons.size();i++) {
		ssb.enter_contour(pgons[i].vertices_begin(),pgons[i].vertices_end());
	}

	boost::shared_ptr<Ss> ss = ssb.construct_skeleton();

	if ( ss ) {
		Ss::Face_const_iterator fi;
		for (fi = ss->faces_begin(); fi!=ss->faces_end(); fi++ ) {
			Polygon_2 tmp;
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
		debug << "construction of skeleton failed\n";
	}
}

typedef enum projection_type_e {
	XY_PROJECTION,
	XZ_PROJECTION,
	YZ_PROJECTION
} projection_type;

// find three noncollinear points within the polygon. set p,q,r.
bool find_noncollinear_pts( Polygon_3 &pgon, Point_3 &p, Point_3 &q, Point_3 &r )
{
	p=pgon[0]; q=p; r=p;
	for ( int i=1;i<pgon.size();i++ ) {
		q = pgon[i];
		if ( p!=q ) break;
	}
	for ( int i=1;i<pgon.size();i++ ) {
		r = pgon[i];
		if (!CGAL::collinear( p, q, r )) break;
	}
	if (CGAL::collinear( p, q, r )) return false;
	return true;
}

tessellater_status tessellate(
        std::vector<Polygon_3> &input_pgon3d,
        std::vector<Polygon_3> &output_pgons3d,
        Tessellation tesstype )
{
	std::stringstream debug;
//	output_pgons3d = input_pgon3d;
//	return TESSELLATER_OK;

	std::vector<Polygon_3> &pgons = input_pgon3d;

	// ---------------- project 3d down to 2d

	// imagine the 'normal' vector to the polygon. now imagine shining
	// a light on it from up, then from the back, and then from the side.
	// imagine the 'shadows' of that normal hit the xy plane, yz plane,
	// or xz plane. the "smallest" shadow tells us which 2d plane we want
	// to project our polygon onto, xy, yz, or xz.

	debug << "finding three non-collinear points, to find normal";
	if (pgons[0].size()<3) return BODY_LACKS_POINTS;
	Point_3 p, q, r;
	if (!find_noncollinear_pts( pgons[0], p, q, r ))
		return BODY_ONLY_COLLINEAR;
	Plane_3 polygon_plane3d( p, q, r );
	Vector_3 normal;
	projection_type projection;
	normal = CGAL::normal(p,q,r);

	TessKernel::FT xyshadow = normal.x()*normal.x()+normal.y()*normal.y();
	TessKernel::FT yzshadow = normal.y()*normal.y()+normal.z()*normal.z();
	TessKernel::FT xzshadow = normal.x()*normal.x()+normal.z()*normal.z();
	TessKernel::FT minshadow = std::min( xyshadow, std::min( yzshadow, xzshadow ) );
	if (minshadow==xyshadow) projection = XY_PROJECTION;
	else if (minshadow==yzshadow) projection = YZ_PROJECTION;
	else projection = XZ_PROJECTION;
	debug << "shadowssize: xy, yz, xz:" << xyshadow << "," << yzshadow << "," << xzshadow << "\n";
	debug << "min shadow: " << minshadow << "\n";

	std::vector<Polygon_2> input_pgon2d;
	std::map<Point_2,Point_3> vertmap;

	debug << "projecting into to 2d\n";
	for (int i=0;i<pgons.size();++i) {
		Polygon_3 tmppoly3d = pgons[i];
		Polygon_2 tmppoly2d;
		for (int j=0;j<tmppoly3d.size();++j) {
			Point_3 tmp3d = tmppoly3d[j];
			Point_2 tmp2d;
			if ( projection == XY_PROJECTION )
				tmp2d = Point_2( tmp3d.x(), tmp3d.y() );
			else if ( projection == XZ_PROJECTION )
				tmp2d = Point_2( tmp3d.x(), tmp3d.z() );
			else if ( projection == YZ_PROJECTION )
				tmp2d = Point_2( tmp3d.y(), tmp3d.z() );
			vertmap[tmp2d] = tmp3d;
			tmppoly2d.push_back( tmp2d );
		}
		input_pgon2d.push_back( tmppoly2d );
	}

	debug << "checking orientation, fixing if necessary\n";
	CGAL::Orientation original_body_orientation = input_pgon2d[0].orientation();
	if (input_pgon2d[0].orientation() != CGAL::COUNTERCLOCKWISE)
		input_pgon2d[0].reverse_orientation();
	for (int i=1;i<input_pgon2d.size();i++) {
		if (input_pgon2d[i].orientation() != CGAL::CLOCKWISE)
			input_pgon2d[i].reverse_orientation();
	}


	debug << "checking for simplicity, self-intersection, etc\n";
	if (!input_pgon2d[0].is_simple()) return BODY_NOT_SIMPLE;
	for (int i=1;i<input_pgon2d.size();i++) {
		if (input_pgon2d[i].size()<3) return HOLE_LACKS_POINTS;
		if (!input_pgon2d[i].is_simple()) return HOLE_NOT_SIMPLE;
		if (!inside(input_pgon2d[i], input_pgon2d[0])) return HOLE_OUTSIDE_BODY;
		for (int j=1;j<input_pgon2d.size();j++) {
			if (inside(input_pgon2d[i], input_pgon2d[j])) return HOLE_INSIDE_HOLE;
		}
	}


	debug << "tessellate into simple polygons-without-holes\n";
	std::vector<Polygon_2> output_pgons2d;
	if (tesstype==CONSTRAINED_DELAUNAY_TRIANGULATION)
		triangulate( input_pgon2d, output_pgons2d );
	else
		skeletonize( input_pgon2d, output_pgons2d );

	debug << "restore original orientation so the 3d normals will be OK\n";
	for (int i=0;i<output_pgons2d.size();i++) {
		if (output_pgons2d[i].orientation() != original_body_orientation)
			output_pgons2d[i].reverse_orientation();
	}

	debug << "project from 2d back into 3d\n";
	for (int i=0;i<output_pgons2d.size();i++) {
		Polygon_2 pgon2d = output_pgons2d[i];
		Polygon_3 pgon3d;
		for (int j=0;j<pgon2d.size();j++) {
			Point_2 pt2d = pgon2d[j];
			Point_3 pt3d;
			if (vertmap.count( pt2d )) {
				pt3d = vertmap[ pt2d ];
			} else {
				Line_3 line3d;
				if ( projection == XY_PROJECTION ) {
					Point_3 tmp3d = Point_3( pt2d.x(), pt2d.y(), 0 );
					line3d = Line_3( tmp3d, Direction_3(0,0,1) );
				} else if ( projection == XZ_PROJECTION ) {
					Point_3 tmp3d = Point_3( pt2d.x(), 0, pt2d.y() );
					line3d = Line_3( tmp3d, Direction_3(0,1,0) );
				} else if ( projection == YZ_PROJECTION ) {
					Point_3 tmp3d = Point_3( 0, pt2d.x(), pt2d.y() );
					line3d = Line_3( tmp3d, Direction_3(1,0,0) );
				}
				CGAL::Object obj = CGAL::intersection( line3d, polygon_plane3d );
				const Point_3 *point_test = CGAL::object_cast<Point_3>(&obj);
				if (point_test)
					pt3d = *point_test;
				else
					return PROJECTION_2D3D_FAILED;
			}
			pgon3d.push_back( pt3d );
		}
		output_pgons3d.push_back( pgon3d );
	}
	std::cout << debug.str() << std::endl;
	return TESSELLATER_OK;
}

tessellater_status tessellate(
        std::vector<CGAL_Polygon_3> &input_pgon3d,
        std::vector<CGAL_Polygon_3> &output_pgons3d,
        Tessellation tesstype )
{
	std::stringstream debug;
	debug << "tess" << tesstype << std::endl;
	CGAL::Cartesian_converter<CGAL_Kernel3,TessKernel> converter1;
	CGAL::Cartesian_converter<TessKernel,CGAL_Kernel3> converter2;
	std::vector<Polygon_3> tess_pgon3d;
	std::vector<Polygon_3> tess_pgon3d_out;
	for (int i=0;i<input_pgon3d.size();i++) {
		CGAL_Polygon_3 oscad_pgon = input_pgon3d[i];
		Polygon_3 tess_pgon;
		for (int j=0;j<oscad_pgon.size();j++) {
			CGAL_Point_3 oscad_point = oscad_pgon[j];
			debug << oscad_point << ",";
			Point_3 tess_point = converter1( oscad_point );
			tess_pgon.push_back( tess_point );
		}
		debug << "\n";
		tess_pgon3d.push_back( tess_pgon );
	}

	tessellater_status s = tessellate( tess_pgon3d, tess_pgon3d_out, tesstype );
	debug << "tessellated. status: " << s << std::endl;

	output_pgons3d.clear();
	for (int i=0;i<tess_pgon3d_out.size();i++) {
		CGAL_Polygon_3 oscad_pgon;
		Polygon_3 tess_pgon = tess_pgon3d_out[i];
		for (int j=0;j<tess_pgon.size();j++) {
			Point_3 tess_point = tess_pgon[j];
			debug << tess_point << ",";
			CGAL_Point_3 oscad_point = converter2( tess_point );
			oscad_pgon.push_back( oscad_point );
		}
		debug << "\n";
		output_pgons3d.push_back( oscad_pgon );
	}
	std::cout << debug.str() << std::endl;
	return TESSELLATER_OK;
}
