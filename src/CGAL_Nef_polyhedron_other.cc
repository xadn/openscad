
#ifdef ENABLE_CGAL

#include "printutils.h"
#include "CGAL_Nef_polyhedron.h"
#include "tess3d.h"

// this class helps build a Polyhedron from a Nef Polyhedron
// See the Polyhedron Builder cgal manual, as well as cgalutils.cc
// CGAL_Build_Polyset.
class Builder : public CGAL::Modifier_base<CGAL_HDS>
{
public:
	std::stringstream debug;
	Tessellation faces_tess;
	Tessellation faces_w_holes_tess;
 	const CGAL_Nef_polyhedron3 &nef;
	Builder(const CGAL_Nef_polyhedron3 &nef) : nef(nef) { }
	void operator()(CGAL_HDS& hds);
};

// called when delegate() is called.
void Builder::operator()(CGAL_HDS& hds)
{
	debug << "builder operator()\n";
	CGAL_Polybuilder B(hds, true);
	std::vector<CGAL_HDS::Vertex::Point> vertices;
	CGAL_Nef_polyhedron3::Vertex_const_iterator vi;
	std::map<CGAL_HDS::Vertex::Point,int> vertmap1;
	int vertcount = 0;

	B.begin_surface(nef.number_of_vertices(), nef.number_of_facets());

	CGAL_forall_vertices( vi, nef ) {
		B.add_vertex( vi->point() );
		vertmap1[ vi->point() ] = vertcount;
		vertcount++;
	}

	CGAL_Nef_polyhedron3::Halffacet_const_iterator hfaceti;
	CGAL_forall_halffacets( hfaceti, nef ) {
		debug << "iterating through next facet...\n";
		if (hfaceti->incident_volume()->mark() == 0) continue;
		//if (hfaceti->mark() == 0) continue;
		CGAL_Nef_polyhedron3::Halffacet_cycle_const_iterator cyclei;
		int cycle_count = 0;
		std::vector<CGAL_Polygon_3> pgons;
		CGAL_forall_facet_cycles_of( cyclei, hfaceti ) {
			CGAL_Nef_polyhedron3::SHalfedge_around_facet_const_circulator c1(cyclei);
			CGAL_Nef_polyhedron3::SHalfedge_around_facet_const_circulator c2(c1);
			CGAL_Polygon_3 pgon;
			CGAL_For_all( c1, c2 ) {
				pgon.push_back(c1->source()->source()->point());
			}
			pgons.push_back( pgon );
		}
		// first polygon = outer contour. next polygons = hole contours.
		std::vector<CGAL_Polygon_3> pgons_without_holes;
		std::cout << "pgon contours input: " << pgons.size() << std::endl;
		if (pgons.size()>1)
			tessellate( pgons, pgons_without_holes, faces_tess );
		else
			tessellate( pgons, pgons_without_holes, faces_w_holes_tess );

		std::cout << "pgon contours output: " << pgons_without_holes.size() << "\n";
		for (int i=0;i<pgons_without_holes.size();i++) {
			CGAL_Polygon_3 pgon = pgons_without_holes[i];
			for (int j=0;j<pgon.size();j++) {
				CGAL_Point_3 point = pgon[j];
				// add new points created by tessellation
				if (vertmap1.count(point)==0) {
					B.add_vertex( point );
					vertmap1[ point ] = vertcount;
					vertcount++;
				}
			}
			B.begin_facet();
			for (int j=0;j<pgon.size();j++) {
				CGAL_Point_3 point = pgon[j];
				int vert_index = vertmap1[ point ] ;
				B.add_vertex_to_facet( vert_index );
			}
			B.end_facet();
		}
	}
	B.end_surface();
}

/* Convert Nef Polyhedron3 to Polyhedron3 with custom face tessellation.
   Faces-without-holes can be tessellated separately from faces-with-holes.
*/
void CGAL_Nef_polyhedron::convertToPolyhedron( CGAL_Polyhedron &P,
	Tessellation faces_tess, Tessellation faces_with_holes_tess ) const
{
	std::cout << "convert " << faces_tess << " " << faces_with_holes_tess << "\n";
	assert(dim==3);
	if (faces_tess == CGAL_NEF_STANDARD && faces_with_holes_tess == CGAL_NEF_STANDARD ) {
  	 	this->p3->convert_to_Polyhedron( P );
	}
	else {
		P.clear();
		Builder builder( *this->p3 );
		builder.faces_tess = faces_tess;
		builder.faces_w_holes_tess = faces_with_holes_tess;
		P.delegate( builder );
	}
}

#endif /* ENABLE_CGAL */



