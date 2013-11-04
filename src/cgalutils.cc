#ifdef ENABLE_CGAL

#include "cgalutils.h"
#include "polyset.h"
#include "printutils.h"
#include "facetess.h"

#include "cgal.h"

#include <map>

PolySet *createPolySetFromPolyhedron(const CGAL_Polyhedron &p)
{
	PolySet *ps = new PolySet();

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
			ps->append_poly();
			ps->append_vertex(x1, y1, z1);
			ps->append_vertex(x2, y2, z2);
			ps->append_vertex(x3, y3, z3);
		} while (hc != hc_end);
	}
	return ps;
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

CGAL_Polyhedron *createPolyhedronFromPolySet(const PolySet &ps)
{
	CGAL_Polyhedron *P = NULL;
	CGAL::Failure_behaviour old_behaviour = CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
	try {
		P = new CGAL_Polyhedron;
		CGAL_Build_PolySet builder(ps);
		P->delegate(builder);
	}
	catch (const CGAL::Assertion_exception &e) {
		PRINTB("CGAL error in CGAL_Build_PolySet: %s", e.what());
		delete P;
		P = NULL;
	}
	CGAL::set_error_behaviour(old_behaviour);
	return P;
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
	for ( vi = explorer.vertices_begin(); vi != explorer.vertices_end(); ++vi )
		if ( explorer.is_standard( vi ) )
			points.push_back( explorer.point( vi ) );
	if (points.size())
		result = CGAL::bounding_box( points.begin(), points.end() );
	return result;
}

void ZRemover::visit( CGAL_Nef_polyhedron3::Halffacet_const_handle hfacet )
{
	PRINTDB( " <!-- ZRemover Halffacet visit. Mark: %i -->",  hfacet->mark());
	if ( hfacet->plane().orthogonal_direction() != this->up ) {
		PRINTD( "  <!-- ZRemover down-facing half-facet. skipping -->");
		PRINTD( " <!-- ZRemover Halffacet visit end-->");
		return;
	}

	// possible optimization - throw out facets that are vertically oriented

	CGAL_Nef_polyhedron3::Halffacet_cycle_const_iterator fci;
	int contour_counter = 0;
	CGAL_forall_facet_cycles_of( fci, hfacet ) {
		if ( fci.is_shalfedge() ) {
			PRINTD( " <!-- ZRemover Halffacet cycle begin -->");
			CGAL_Nef_polyhedron3::SHalfedge_around_facet_const_circulator c1(fci), cend(c1);
			std::vector<CGAL_Nef_polyhedron2::Explorer::Point> contour;
			CGAL_For_all( c1, cend ) {
				CGAL_Nef_polyhedron3::Point_3 point3d = c1->source()->target()->point();
				CGAL_Nef_polyhedron2::Explorer::Point point2d(CGAL::to_double(point3d.x()),
																											CGAL::to_double(point3d.y()));
				contour.push_back( point2d );
			}
			if (contour.size()==0) continue;

			PRINTDB( " <!-- is_simple_2: %i -->", CGAL::is_simple_2( contour.begin(), contour.end() ) );

			tmpnef2d.reset( new CGAL_Nef_polyhedron2( contour.begin(), contour.end(), boundary ) );

			if ( contour_counter == 0 ) {
				PRINTDB( " <!-- contour is a body. make union(). %i points. -->", contour.size());
				*(output_nefpoly2d) += *(tmpnef2d);
			} else {
				PRINTDB( " <!-- contour is a hole. make intersection(). %i points -->", contour.size());
				*(output_nefpoly2d) *= *(tmpnef2d);
			}

			/*PRINTD( "\n<!-- ======== output tmp nef: ==== -->");
				<< OpenSCAD::dump_svg( *tmpnef2d ) << "\n"
				<< "\n<!-- ======== output accumulator: ==== -->");
				<< OpenSCAD::dump_svg( *output_nefpoly2d ) << "\n";*/

			contour_counter++;
		} else {
			PRINTD(" <!-- ZRemover trivial facet cycle skipped -->");
		}
		PRINTD(" <!-- ZRemover Halffacet cycle end -->");
	}
	PRINTD(" <!-- ZRemover Halffacet visit end -->");
}

using namespace OpenSCAD;

class TessBuilder : public CGAL::Modifier_base<CGAL_HDS>
{
public:
	facetess::tesstype faces_tess;
	facetess::tesstype faces_w_holes_tess;
	const CGAL_Nef_polyhedron3 &nef;
	TessBuilder(const CGAL_Nef_polyhedron3 &N);
	void operator()(CGAL_HDS& hds);
	~TessBuilder();
};

TessBuilder::TessBuilder(const CGAL_Nef_polyhedron3 &N):nef(N)
{
	faces_tess=facetess::CGAL_NEF_STANDARD;
	faces_w_holes_tess=facetess::CGAL_NEF_STANDARD;
}

TessBuilder::~TessBuilder()
{
}

// called when Polyhedron.delegate() is called.
void TessBuilder::operator()(CGAL_HDS& hds)
{
	PRINTD("TessBuilder operator()");
	hds.clear();
	CGAL_Polybuilder B(hds, true);
	std::vector<CGAL_HDS::Vertex::Point> vertices;
	CGAL_Nef_polyhedron3::Vertex_const_iterator vi;
	std::map<CGAL_HDS::Vertex::Point,int> vertmap1;
	std::map<CGAL_Nef_polyhedron3::Vertex_const_iterator,int> vertmap1a;
	int vertcount = 0;

	B.begin_surface(nef.number_of_vertices(), nef.number_of_facets());

	CGAL_forall_vertices( vi, nef ) {
		B.add_vertex( vi->point() );
		vertmap1[ vi->point() ] = vertcount;
		vertmap1a[ vi ] = vertcount;
		vertcount++;
	}

	CGAL_Nef_polyhedron3::Halffacet_const_iterator hfaceti;
	CGAL_forall_halffacets( hfaceti, nef ) {
		if (hfaceti->incident_volume()->mark() == 0) continue;
		PRINTD("iterating through next facet...");
		//if (hfaceti->mark() == 0) continue;
		CGAL_Nef_polyhedron3::Halffacet_cycle_const_iterator cyclei;
		//int cycle_count = 0;
		std::vector<CGAL_Polygon_3> pgons;
		CGAL_forall_facet_cycles_of( cyclei, hfaceti ) {
			CGAL_Nef_polyhedron3::SHalfedge_around_facet_const_circulator c1(cyclei);
			CGAL_Nef_polyhedron3::SHalfedge_around_facet_const_circulator c2(c1);
			CGAL_Polygon_3 pgon;
			CGAL_For_all( c1, c2 ) {
				pgon.push_back(c1->source()->source()->point());
				//pgon.push_back(c1->source()->center_vertex());
			}
			pgons.push_back( pgon );
		}
		// first polygon = outer contour. next polygons = hole contours.
		std::vector<CGAL_Polygon_3> pgons_without_holes;
		PRINTDB("pgon contours input: %i", pgons.size());
		if (pgons.size()>1)
			facetess::tessellate( pgons, pgons_without_holes, faces_tess );
		else
			facetess::tessellate( pgons, pgons_without_holes, faces_w_holes_tess );

		PRINTDB("pgon contours output: %i", pgons_without_holes.size());
		for (size_t i=0;i<pgons_without_holes.size();i++) {
			CGAL_Polygon_3 pgon = pgons_without_holes[i];
			for (size_t j=0;j<pgon.size();j++) {
				CGAL_Point_3 point = pgon[j];
				// add new points created by tessellation
				if (vertmap1.count(point)==0) {
					B.add_vertex( point );
					vertmap1[ point ] = vertcount;
					vertcount++;
				}
			}
			B.begin_facet();
			for (size_t j=0;j<pgon.size();j++) {
				CGAL_Point_3 point = pgon[j];
				int vert_index = vertmap1[ point ] ;
				B.add_vertex_to_facet( vert_index );
			}
			B.end_facet();
		}
	}
	B.end_surface();
	hds.normalize_border();
	PRINTD("TessBuilder operator() finished");
}

/* Convert Nef Polyhedron3 to Polyhedron3 with custom face tessellation.
   Faces-without-holes can be tessellated separately from faces-with-holes.
*/
void nef3_to_polyhedron( const CGAL_Nef_polyhedron3 &N, CGAL_Polyhedron &P,
	facetess::tesstype faces_tess, facetess::tesstype faces_with_holes_tess )
{
	if (faces_tess == facetess::CGAL_NEF_STANDARD && faces_with_holes_tess == facetess::CGAL_NEF_STANDARD ) {
  		N.convert_to_Polyhedron( P );
	}
	else {
		P.clear();
		TessBuilder builder( N );
		builder.faces_tess = faces_tess;
		builder.faces_w_holes_tess = faces_with_holes_tess;
		P.delegate( builder );
	}
}

#endif /* ENABLE_CGAL */

