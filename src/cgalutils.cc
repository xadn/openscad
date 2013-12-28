#ifdef ENABLE_CGAL

#include "cgalutils.h"
#include "polyset.h"
#include "printutils.h"
#include "dxftess.h"
#include "cgal.h"

#include <map>


/////// Tessellation begin

/*

This is our custom tessellator of Nef Polyhedron faces. The problem with 
Nef faces is that sometimes the 'default' tessellator of Nef Polyhedron 
doesnt work. This is particularly true with situations where the polygon 
face is not, actually, 'simple', according to CGAL itself. This can 
occur on a bad quality STL import but also for other reasons. The 
resulting Nef face will appear to the average human eye as an ordinary, 
simple polygon... but in reality it has multiple edges that are 
slightly-out-of-alignment and sometimes they backtrack on themselves.

When the triangulator is fed a polygon with self-intersecting edges, 
it's default behavior is to throw an exception. The other terminology 
for this is to say that the 'constraints' in the triangulation are 
'intersecting'. The 'constraints' represent the edges of the polygon. 
The 'triangulation' is the covering of all the polygon points with 
triangles.

How do we allow interseting constraints during triangulation? We use an 
'Itag' for the triangulation, per the CGAL docs. This allows the 
triangulator to run without throwing an exception when it encounters 
self-intersecting polygon edges. The trick here is that when it finds
an intersection, it actually creates a new point. 

The triangulator creates new points in 2d, but they aren't matched to 
any 3d points on our 3d polygon plane. (The plane of the Nef face). How 
to fix this problem? We actually 'project back up' into 3d plane from the
2d point. This is handled in the 'deproject()' function.

But there is another problem. The intersecting-constraints create 
polygons that are not simple. This means that CGAL functions like 
orientation_2() and bounded_side() simply will not work because they all
require input polygons to pass the 'is_simple2()' test. We have to use 
alternatives.

*/

#include <CGAL/Delaunay_mesher_no_edge_refinement_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>

typedef CGAL_Kernel3 Kernel;
typedef typename CGAL::Triangulation_vertex_base_2<Kernel> Vb;
//typedef typename CGAL::Constrained_triangulation_face_base_2<Kernel> Fb;
typedef CGAL::Delaunay_mesh_face_base_2<Kernel> Fb;
typedef typename CGAL::Triangulation_data_structure_2<Vb,Fb> TDS;
typedef CGAL::Exact_intersections_tag ITAG;
typedef typename CGAL::Constrained_Delaunay_triangulation_2<Kernel,TDS,ITAG> CDT;
//typedef typename CGAL::Constrained_Delaunay_triangulation_2<Kernel,TDS> CDT;

typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point CDTPoint;

typedef CGAL::Ray_2<Kernel> CGAL_Ray_2;
typedef CGAL::Line_3<Kernel> CGAL_Line_3;
typedef CGAL::Point_2<Kernel> CGAL_Point_2;
typedef CGAL::Vector_2<Kernel> CGAL_Vector_2;
typedef CGAL::Segment_2<Kernel> CGAL_Segment_2;
typedef std::vector<CGAL_Point_3> CGAL_Polygon_3;
typedef CGAL::Direction_2<Kernel> CGAL_Direction_2;
typedef CGAL::Direction_3<Kernel> CGAL_Direction_3;

typedef enum { XYPLANE, YZPLANE, XZPLANE, FLIP, NONE } projection_t;

// this is how we make 3d points appear as though they were 2d points to
//the tessellation algorithm.
CGAL_Point_2 get_projected_point( CGAL_Point_3 &p3, projection_t projection ) {
	NT3 x,y,tmp;
	if      (projection & XYPLANE) { x = p3.x(); y = p3.y(); }
	else if (projection & XZPLANE) { x = p3.x(); y = p3.z(); }
	else if (projection & YZPLANE) { x = p3.y(); y = p3.z(); }
	if (projection & FLIP) return CGAL_Point_2( y,x );
	return CGAL_Point_2( x,y );
}

// given 2d point, 3d plane, and 3d->2d projection, 'deproject' from
// 2d back onto a point on the 3d plane. true on failure, false on success
bool deproject( CGAL_Point_2 &p2, projection_t &projection, CGAL_Plane_3 &plane, CGAL_Point_3 &p3 )
{
	NT3 x,y;
	CGAL_Line_3 l;
	CGAL_Point_3 p;
	CGAL_Point_2 pf;
	if (projection & FLIP) pf = CGAL_Point_2(p2.y(),p2.x());
	else CGAL_Point_2(p2.x(),p2.y());
        if (projection & XYPLANE) {
		p = CGAL_Point_3( pf.x(), pf.y(), 0 );
		l = CGAL_Line_3( p, CGAL_Direction_3(0,0,1) );
	} else if (projection & XZPLANE) {
		p = CGAL_Point_3( pf.x(), 0, pf.y() );
		l = CGAL_Line_3( p, CGAL_Direction_3(0,1,0) );
	} else if (projection & YZPLANE) {
		p = CGAL_Point_3( 0, pf.x(), pf.y() );
		l = CGAL_Line_3( p, CGAL_Direction_3(1,0,0) );
	}
	CGAL::Object obj = CGAL::intersection( l, plane );
	const CGAL_Point_3 *point_test = CGAL::object_cast<CGAL_Point_3>(&obj);
	if (point_test) {
		p3 = *point_test;
		return false;
	}
	PRINT("deproject failure");
	return true;
}

// simple criteria guarantees the algorithm will terminate (i.e. not lock
// up and freeze the program)
template <class T> class DummyCriteria {
public:
        typedef double Quality;
        class Is_bad {
        public:
                CGAL::Mesh_2::Face_badness operator()(const Quality) const {
                        return CGAL::Mesh_2::NOT_BAD;
                }
                CGAL::Mesh_2::Face_badness operator()(const typename T::Face_handle&, Quality&q) const {
                        q = 1;
                        return CGAL::Mesh_2::NOT_BAD;
                }
        };
        Is_bad is_bad_object() const { return Is_bad(); }
};

NT3 sign( NT3 &n ) { if (n>0) return 1; if (n<0) return -1; return 0; }

/* wedge, also related to 'determinant', 'signed parallelogram area', 
'side', 'turn', 'winding', '2d portion of cross-product', etc etc. 
can tell you whether v1 is 'counterclockwise' or 'clockwise' from v2, 
based on the sign of the result. */
NT3 wedge( CGAL_Vector_2 &v1, CGAL_Vector_2 &v2 ) {
	return v1.x()*v2.y()-v2.x()*v1.y();
}

/* given a point and a possibly non-simple polygon, determine if the
point is inside the polygon or not, using the given winding rule. note
that even_odd is not implemented. */
typedef enum { NONZERO_WINDING, EVEN_ODD } winding_rule_t;
bool inside(CGAL_Point_2 &p1,std::vector<CDTPoint> &pgon, winding_rule_t winding_rule)
{
	NT3 overall_winding = NT3(0);
	CGAL_Point_2 p2;
	CGAL_Ray_2 eastray(p1,CGAL_Direction_2(1,0));
	for (size_t i=0;i<pgon.size();i++) {
		CGAL_Point_2 tail = pgon[i];
		CGAL_Point_2 head = pgon[(i+1)%pgon.size()];
		CGAL_Segment_2 seg( tail, head );
		CGAL::Object obj = intersection( eastray, seg );
		const CGAL_Point_2 *point_test = CGAL::object_cast<CGAL_Point_2>(&obj);
		if (point_test) {
			p2 = *point_test;
			CGAL_Vector_2 v1( p1, p2 );
			CGAL_Vector_2 v2( p2, head );
			NT3 this_winding = wedge( v1, v2 );
			overall_winding += sign(this_winding);
		} else {
			continue;
		}
	}
	if (overall_winding != NT3(0) && winding_rule == NONZERO_WINDING ) return true;
	return false;
}

/* given a single near-planar 3d polygon with holes, tessellate into a 
sequence of polygons without holes. as of writing, this means conversion 
into a sequence of triangles. the plane is the same as the polygon and 
holes. */
bool tessellate_3d_face_with_holes( std::vector<CGAL_Polygon_3> &polygons, std::vector<CGAL_Polygon_3> &triangles, CGAL_Plane_3 &plane )
{
	if (polygons.size()==1 && polygons[0].size()==3) {
		PRINT("input polygon has 3 points. shortcut tessellation.");
		//triangles.push_back( polygons[0] );
		//return false;
	}
	bool err = false;
	CDT cdt;
	std::map<CDTPoint,CGAL_Point_3> vertmap;

	PRINT("finding good projection");
	projection_t goodproj;
        NT3 qxy = plane.a()*plane.a()+plane.b()*plane.b();
        NT3 qyz = plane.b()*plane.b()+plane.c()*plane.c();
        NT3 qxz = plane.a()*plane.a()+plane.c()*plane.c();
        NT3 min = std::min(qxy,std::min(qyz,qxz));
        if (min==qxy) {
		goodproj = XYPLANE;
	} else if (min==qyz) {
		goodproj = YZPLANE;
	} else if (min==qxz) {
		goodproj = XZPLANE;
	}
	if (orientation(plane,a,b,c)==CLOCKWISE) goodproj |= FLIP;

	PRINTB("plane %s",plane );
	PRINTB("proj: %i",goodproj);
	PRINT("Inserting points and edges into Constrained Delaunay Triangulation");
	std::vector< std::vector<CGAL_Point_2> > polygons2d;
	for (size_t i=0;i<polygons.size();i++) {
	        std::vector<Vertex_handle> vhandles;
		std::vector<CGAL_Point_2> polygon2d;
		for (size_t j=0;j<polygons[i].size();j++) {
			CGAL_Point_3 p3 = polygons[i][j];
			CGAL_Point_2 p2 = get_projected_point( p3, goodproj );
			CDTPoint cdtpoint = CDTPoint( p2.x(), p2.y() );
			vertmap[ cdtpoint ] = p3;
			Vertex_handle vh = cdt.push_back( cdtpoint );
			vhandles.push_back( vh );
			polygon2d.push_back( p2 );
		}
		polygons2d.push_back( polygon2d );
		for (size_t k=0;k<vhandles.size();k++) {
			int vindex1 = (k+0);
			int vindex2 = (k+1)%vhandles.size();
			CGAL::Failure_behaviour old_behaviour = CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
			try {
				cdt.insert_constraint( vhandles[vindex1], vhandles[vindex2] );
			} catch (const CGAL::Failure_exception &e) {
				PRINTB("WARNING: Constraint insertion failure %s", e.what());
			}
			CGAL::set_error_behaviour(old_behaviour);
		}
	}

	size_t numholes = polygons2d.size()-1;
	PRINTB("seeding %i holes",numholes);
	std::list<CDTPoint> list_of_seeds;
	for (size_t i=1;i<polygons2d.size();i++) {
		std::vector<CGAL_Point_2> &pgon = polygons2d[i];
		for (size_t j=0;j<pgon.size();j++) {
			CGAL_Point_2 p1 = pgon[(j+0)];
			CGAL_Point_2 p2 = pgon[(j+1)%pgon.size()];
			CGAL_Point_2 p3 = pgon[(j+2)%pgon.size()];
			CGAL_Point_2 mp = CGAL::midpoint(p1,CGAL::midpoint(p2,p3));
			if (inside(mp,pgon,NONZERO_WINDING)) {
				CDTPoint cdtpt( mp.x(), mp.y() );
				list_of_seeds.push_back( cdtpt );
				break;
			}
		}
	}
	std::list<CDTPoint>::iterator li = list_of_seeds.begin();
	for (;li!=list_of_seeds.end();li++) {
		//PRINTB("seed %s",*li);
		double x = CGAL::to_double( li->x() );
		double y = CGAL::to_double( li->y() );
		PRINTB("seed %f,%f",x%y);
	}
	PRINT("seeding done");

	PRINT( "meshing" );
	CGAL::refine_Delaunay_mesh_2_without_edge_refinement( cdt,
		list_of_seeds.begin(), list_of_seeds.end(),
		DummyCriteria<CDT>() );

	PRINT("meshing done");
	// fails w is_simple
	//CGAL::Orientation original_orientation =
	//	CGAL::orientation_2( orienpgon.begin(), orienpgon.end() );

	CDT::Finite_faces_iterator fit;
	for( fit=cdt.finite_faces_begin(); fit!=cdt.finite_faces_end(); fit++ )
	{
		if(fit->is_in_domain()) {
			CDTPoint p1 = cdt.triangle( fit )[0];
			CDTPoint p2 = cdt.triangle( fit )[1];
			CDTPoint p3 = cdt.triangle( fit )[2];
			CGAL_Point_3 cp1,cp2,cp3;
			CGAL_Polygon_3 pgon;
			if (vertmap.count(p1)) cp1 = vertmap[p1];
			else err = deproject( p1, goodproj, plane, cp1 );
			if (vertmap.count(p2)) cp2 = vertmap[p2];
			else err = deproject( p2, goodproj, plane, cp2 );
			if (vertmap.count(p3)) cp3 = vertmap[p3];
			else err = deproject( p3, goodproj, plane, cp3 );
			if (err) PRINT("WARNING: 2d->3d deprojection failure");
			pgon.push_back( cp1 );
			pgon.push_back( cp2 );
			pgon.push_back( cp3 );
                        triangles.push_back( pgon );
                }
        }

	PRINTB("built %i triangles\n",triangles.size());
	return err;
}
/////// Tessellation end

/* Create a PolySet from a Nef Polyhedron 3. return false on success, 
true on failure. The trick to this is that Nef Polyhedron3 faces have 
'holes' in them. . . while PolySet (and many other 3d polyhedron 
formats) do not allow for holes in their faces. The function documents 
the method used to deal with this */
bool createPolySetFromNefPolyhedron3(const CGAL_Nef_polyhedron3 &N, PolySet &ps)
{
	bool err = false;
	CGAL_Nef_polyhedron3::Halffacet_const_iterator hfaceti;
	CGAL_forall_halffacets( hfaceti, N ) {
		CGAL_Plane_3 plane( hfaceti->plane() );
		std::vector<CGAL_Polygon_3> polygons;
		// the 0-mark-volume is the 'empty' volume of space. skip it.
		if (hfaceti->incident_volume()->mark() == 0) continue;
		CGAL_Nef_polyhedron3::Halffacet_cycle_const_iterator cyclei;
		CGAL_forall_facet_cycles_of( cyclei, hfaceti ) {
			CGAL_Nef_polyhedron3::SHalfedge_around_facet_const_circulator c1(cyclei);
			CGAL_Nef_polyhedron3::SHalfedge_around_facet_const_circulator c2(c1);
			CGAL_Polygon_3 polygon;
			CGAL_For_all( c1, c2 ) {
				CGAL_Point_3 p = c1->source()->center_vertex()->point();
				polygon.push_back( p );
			}
			polygons.push_back( polygon );
		}
		/* at this stage, we have a sequence of polygons. the first
		is the "outside edge' or 'body' or 'border', and the rest of the
		polygons are 'holes' within the first. there are several
		options here to get rid of the holes. we choose to go ahead
		and tessellate the polygon in 3d before adding to the PolySet.*/
		std::vector<CGAL_Polygon_3> triangles;
		bool err = tessellate_3d_face_with_holes( polygons, triangles, plane );
		if (!err) for (size_t i=0;i<triangles.size();i++) {
			if (triangles[i].size()!=3) {
				PRINT("WARNING: triangle != 3 points");
				continue;
			}
			ps.append_poly();
			for (size_t j=0;j<3;j++) {
				double x1,y1,z1;
				x1 = CGAL::to_double( triangles[i][j].x() );
				y1 = CGAL::to_double( triangles[i][j].y() );
				z1 = CGAL::to_double( triangles[i][j].z() );
				ps.append_vertex(x1,y1,z1);
			}
		}
	}
	return err;
}

bool createPolySetFromPolyhedron(const CGAL_Polyhedron &p, PolySet &ps)
{
	bool err = false;
	typedef CGAL_Polyhedron::Vertex                                 Vertex;
	typedef CGAL_Polyhedron::Vertex_const_iterator                  VCI;
	typedef CGAL_Polyhedron::Facet_const_iterator                   FCI;
	typedef CGAL_Polyhedron::Halfedge_around_facet_const_circulator HFCC;
		
	for (FCI fi = p.facets_begin(); fi != p.facets_end(); ++fi) {
		HFCC hc = fi->facet_begin();
		HFCC hc_end = hc;
		Vertex v1, v2, v3;
		v1 = *VCI((hc++)->vertex());
		v3 = *VCI((hc++)->vertex());
		do {
			v2 = v3;
			v3 = *VCI((hc++)->vertex());
			double x1 = CGAL::to_double(v1.point().x());
			double y1 = CGAL::to_double(v1.point().y());
			double z1 = CGAL::to_double(v1.point().z());
			double x2 = CGAL::to_double(v2.point().x());
			double y2 = CGAL::to_double(v2.point().y());
			double z2 = CGAL::to_double(v2.point().z());
			double x3 = CGAL::to_double(v3.point().x());
			double y3 = CGAL::to_double(v3.point().y());
			double z3 = CGAL::to_double(v3.point().z());
			ps.append_poly();
			ps.append_vertex(x1, y1, z1);
			ps.append_vertex(x2, y2, z2);
			ps.append_vertex(x3, y3, z3);
		} while (hc != hc_end);
	}
	return err;
}

#undef GEN_SURFACE_DEBUG

class CGAL_Build_PolySet : public CGAL::Modifier_base<CGAL_HDS>
{
public:
	typedef CGAL_HDS::Vertex::Point CGALPoint;

	const PolySet &ps;
	CGAL_Build_PolySet(const PolySet &ps) : ps(ps) { }

	void operator()(CGAL_HDS& hds)
	{
		CGAL_Polybuilder B(hds, true);

		std::vector<CGALPoint> vertices;
		Grid3d<int> vertices_idx(GRID_FINE);

		for (size_t i = 0; i < ps.polygons.size(); i++) {
			const PolySet::Polygon *poly = &ps.polygons[i];
			for (size_t j = 0; j < poly->size(); j++) {
				const Vector3d &p = poly->at(j);
				if (!vertices_idx.has(p[0], p[1], p[2])) {
					vertices_idx.data(p[0], p[1], p[2]) = vertices.size();
					vertices.push_back(CGALPoint(p[0], p[1], p[2]));
				}
			}
		}

		B.begin_surface(vertices.size(), ps.polygons.size());
#ifdef GEN_SURFACE_DEBUG
		printf("=== CGAL Surface ===\n");
#endif

		for (size_t i = 0; i < vertices.size(); i++) {
			const CGALPoint &p = vertices[i];
			B.add_vertex(p);
#ifdef GEN_SURFACE_DEBUG
			printf("%d: %f %f %f\n", i, p.x().to_double(), p.y().to_double(), p.z().to_double());
#endif
		}

		for (size_t i = 0; i < ps.polygons.size(); i++) {
			const PolySet::Polygon *poly = &ps.polygons[i];
			std::map<int,int> fc;
			bool facet_is_degenerated = false;
			for (size_t j = 0; j < poly->size(); j++) {
				const Vector3d &p = poly->at(j);
				int v = vertices_idx.data(p[0], p[1], p[2]);
				if (fc[v]++ > 0)
					facet_is_degenerated = true;
			}
			
			if (!facet_is_degenerated)
				B.begin_facet();
#ifdef GEN_SURFACE_DEBUG
			printf("F:");
#endif
			for (size_t j = 0; j < poly->size(); j++) {
				const Vector3d &p = poly->at(j);
#ifdef GEN_SURFACE_DEBUG
				printf(" %d (%f,%f,%f)", vertices_idx.data(p[0], p[1], p[2]), p[0], p[1], p[2]);
#endif
				if (!facet_is_degenerated)
					B.add_vertex_to_facet(vertices_idx.data(p[0], p[1], p[2]));
			}
#ifdef GEN_SURFACE_DEBUG
			if (facet_is_degenerated)
				printf(" (degenerated)");
			printf("\n");
#endif
			if (!facet_is_degenerated)
				B.end_facet();
		}

#ifdef GEN_SURFACE_DEBUG
		printf("====================\n");
#endif
		B.end_surface();

		#undef PointKey
	}
};

bool createPolyhedronFromPolySet(const PolySet &ps, CGAL_Polyhedron &p)
{
	bool err = false;
	CGAL::Failure_behaviour old_behaviour = CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
	try {
		CGAL_Build_PolySet builder(ps);
		p.delegate(builder);
	}
	catch (const CGAL::Assertion_exception &e) {
		PRINTB("CGAL error in CGAL_Build_PolySet: %s", e.what());
		err = true;
	}
	CGAL::set_error_behaviour(old_behaviour);
	return err;
}

CGAL_Iso_cuboid_3 bounding_box( const CGAL_Nef_polyhedron3 &N )
{
	CGAL_Iso_cuboid_3 result(0,0,0,0,0,0);
	CGAL_Nef_polyhedron3::Vertex_const_iterator vi;
	std::vector<CGAL_Nef_polyhedron3::Point_3> points;
	// can be optimized by rewriting bounding_box to accept vertices
	CGAL_forall_vertices( vi, N )
		points.push_back( vi->point() );
	if (points.size())
		result = CGAL::bounding_box( points.begin(), points.end() );
	return result;
}

CGAL_Iso_rectangle_2e bounding_box( const CGAL_Nef_polyhedron2 &N )
{
	CGAL_Iso_rectangle_2e result(0,0,0,0);
	CGAL_Nef_polyhedron2::Explorer explorer = N.explorer();
	CGAL_Nef_polyhedron2::Explorer::Vertex_const_iterator vi;
	std::vector<CGAL_Point_2e> points;
	// can be optimized by rewriting bounding_box to accept vertices
	for ( vi = explorer.vertices_begin(); vi != explorer.vertices_end(); ++vi )
		if ( explorer.is_standard( vi ) )
			points.push_back( explorer.point( vi ) );
	if (points.size())
		result = CGAL::bounding_box( points.begin(), points.end() );
	return result;
}

void ZRemover::visit( CGAL_Nef_polyhedron3::Halffacet_const_handle hfacet )
{
	log << " <!-- ZRemover Halffacet visit. Mark: " << hfacet->mark() << " -->\n";
	if ( hfacet->plane().orthogonal_direction() != this->up ) {
		log << "  <!-- ZRemover down-facing half-facet. skipping -->\n";
		log << " <!-- ZRemover Halffacet visit end-->\n";
		return;
	}

	// possible optimization - throw out facets that are vertically oriented

	CGAL_Nef_polyhedron3::Halffacet_cycle_const_iterator fci;
	int contour_counter = 0;
	CGAL_forall_facet_cycles_of( fci, hfacet ) {
		if ( fci.is_shalfedge() ) {
			log << " <!-- ZRemover Halffacet cycle begin -->\n";
			CGAL_Nef_polyhedron3::SHalfedge_around_facet_const_circulator c1(fci), cend(c1);
			std::vector<CGAL_Nef_polyhedron2::Explorer::Point> contour;
			CGAL_For_all( c1, cend ) {
				CGAL_Nef_polyhedron3::Point_3 point3d = c1->source()->target()->point();
				CGAL_Nef_polyhedron2::Explorer::Point point2d(CGAL::to_double(point3d.x()),
																											CGAL::to_double(point3d.y()));
				contour.push_back( point2d );
			}
			if (contour.size()==0) continue;

			log << " <!-- is_simple_2:" << CGAL::is_simple_2( contour.begin(), contour.end() ) << " --> \n";

			tmpnef2d.reset( new CGAL_Nef_polyhedron2( contour.begin(), contour.end(), boundary ) );

			if ( contour_counter == 0 ) {
				log << " <!-- contour is a body. make union(). " << contour.size() << " points. -->\n" ;
				*(output_nefpoly2d) += *(tmpnef2d);
			} else {
				log << " <!-- contour is a hole. make intersection(). " << contour.size() << " points. -->\n";
				*(output_nefpoly2d) *= *(tmpnef2d);
			}

			/*log << "\n<!-- ======== output tmp nef: ==== -->\n"
				<< OpenSCAD::dump_svg( *tmpnef2d ) << "\n"
				<< "\n<!-- ======== output accumulator: ==== -->\n"
				<< OpenSCAD::dump_svg( *output_nefpoly2d ) << "\n";*/

			contour_counter++;
		} else {
			log << " <!-- ZRemover trivial facet cycle skipped -->\n";
		}
		log << " <!-- ZRemover Halffacet cycle end -->\n";
	}
	log << " <!-- ZRemover Halffacet visit end -->\n";
}

#endif /* ENABLE_CGAL */

