/*
 *  OpenSCAD (www.openscad.org)
 *  Copyright (C) 2009-2011 Clifford Wolf <clifford@clifford.at> and
 *                          Marius Kintel <marius@kintel.net>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  As a special exception, you have permission to link this program
 *  with the CGAL library and distribute executables, as long as you
 *  follow the requirements of the GNU GPL in regard to all of the
 *  software in the executable aside from CGAL.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "export.h"
#include "printutils.h"
#include "polyset.h"
#include "dxfdata.h"
#include <vector>
#include <string>
#include <boost/foreach.hpp>

#ifdef ENABLE_CGAL
#include "CGAL_Nef_polyhedron.h"
#include "cgal.h"
#include "cgalutils.h"

#include <Eigen/Geometry>

struct stl_point {
	std::string x, y, z;
	std::string const tostr() const;
};

struct stl_triangle {
	stl_point p1, p2, p3, normal;
};

std::ostream&  operator<<( std::ostream &stream, const stl_point &p ) {
	stream << p.x << " " << p.y << " " << p.z;
	return stream;
}

std::string const stl_point::tostr() const
{
	std::stringstream s;
	s << *this;
	return s.str();
}

bool operator!=( const stl_point &p, const stl_point &p2 ) {
	return (p.x!=p2.x || p.y!=p2.y || p.z!=p2.z);
}

// for map() only. not geometric nor algebraic sense of less-than
bool operator<( const stl_point &p, const stl_point &p2 ) {
	return (p.tostr() < p2.tostr());
}

stl_point cgal_point_to_stl_point( const CGAL_Point_3 &p )
{
	stl_point pt;
	pt.x = boost::lexical_cast<std::string>( CGAL::to_double( p.x() ) );
	pt.y = boost::lexical_cast<std::string>( CGAL::to_double( p.y() ) );
	pt.z = boost::lexical_cast<std::string>( CGAL::to_double( p.z() ) );
	return pt;
}

std::vector<stl_triangle> get_cgal_polyhedron_triangles( const CGAL_Polyhedron &P )
{
	std::vector<stl_triangle> triangles;
/*
	CGAL::Failure_behaviour old_behaviour = CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
	try {
	CGAL_Polyhedron P;
	root_N->p3->convert_to_Polyhedron(P);
*/

	typedef CGAL_Polyhedron::Vertex                                 Vertex;
	typedef CGAL_Polyhedron::Vertex_const_iterator                  VCI;
	typedef CGAL_Polyhedron::Facet_const_iterator                   FCI;
	typedef CGAL_Polyhedron::Halfedge_around_facet_const_circulator HFCC;

	setlocale(LC_NUMERIC, "C"); // Ensure radix is . (not ,) in output

	for (FCI fi = P.facets_begin(); fi != P.facets_end(); ++fi) {
		HFCC hc = fi->facet_begin();
		HFCC hc_end = hc;
		Vertex v1, v2, v3;
		v1 = *VCI((hc++)->vertex());
		v3 = *VCI((hc++)->vertex());
		do {
			v2 = v3;
			v3 = *VCI((hc++)->vertex());
			CGAL_Point_3 p1 = v1.point();
			CGAL_Point_3 p2 = v2.point();
			CGAL_Point_3 p3 = v3.point();
			stl_triangle tri;
			tri.p1 = cgal_point_to_stl_point( p1 );
			tri.p2 = cgal_point_to_stl_point( p2 );
			tri.p3 = cgal_point_to_stl_point( p3 );
			tri.normal.x = "1";
			tri.normal.y = "0";
			tri.normal.z = "0";
			if (tri.p1 != tri.p2 && tri.p1 != tri.p3 && tri.p2 != tri.p3) {
				// The above condition ensures that there are 3 distinct vertices, but
				// they may be collinear. If they are, the unit normal is meaningless
				// so the default value of "1 0 0" can be used. If the vertices are not
				// collinear then the unit normal must be calculated from the
				// components.
				if (!CGAL::collinear(p1, p2, p3)) {
					CGAL_Polyhedron::Traits::Vector_3 normal = CGAL::normal(p1, p2, p3);
					double nx = CGAL::sign(normal.x()) * sqrt(CGAL::to_double(normal.x()*normal.x()/normal.squared_length()));
					double ny = CGAL::sign(normal.y()) * sqrt(CGAL::to_double(normal.y()*normal.y()/normal.squared_length()));
					double nz = CGAL::sign(normal.z()) * sqrt(CGAL::to_double(normal.z()*normal.z()/normal.squared_length()));
					tri.normal.x = boost::lexical_cast<std::string>( nx );
					tri.normal.y = boost::lexical_cast<std::string>( ny );
					tri.normal.z = boost::lexical_cast<std::string>( nz );
				}
				triangles.push_back(tri);
			}
		} while (hc != hc_end);
	}

	setlocale(LC_NUMERIC, "");      // Set default locale

	return triangles;
}

std::vector<stl_triangle> get_cgal_nef_poly_triangles( const CGAL_Nef_polyhedron3 &N, OpenSCAD::facetess::tesstype tess )
{
	if (tess!=OpenSCAD::facetess::CGAL_NEF_STANDARD && tess!=OpenSCAD::facetess::CONSTRAINED_DELAUNAY_TRIANGULATION) {
		PRINT("WARNING: This format needs triangle tessellation. Using CGAL default.");
		tess = OpenSCAD::facetess::CGAL_NEF_STANDARD; // force use of tessellator that generates triangles
	}
	std::vector<stl_triangle> triangles;
	CGAL::Failure_behaviour old_behaviour = CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
	try {
		CGAL_Polyhedron P;
		nef3_to_polyhedron( N, P, tess, tess );
		triangles = get_cgal_polyhedron_triangles( P );
	}
	catch (CGAL::Assertion_exception e) {
		PRINTB("CGAL error in CGAL_Nef_polyhedron3::convert_to_Polyhedron(): %s", e.what());
	}
	CGAL::set_error_behaviour(old_behaviour);
	return triangles;
}


/*!
  Saves the current 3D CGAL Nef polyhedron as STL to the given file.
  The file must be open.
 */
void export_stl(CGAL_Nef_polyhedron *root_N, std::ostream &output, OpenSCAD::facetess::tesstype tess)
{
	output << "solid OpenSCAD_Model\n";
	assert( root_N );
	assert( root_N->p3 );
	std::vector<stl_triangle> triangles = get_cgal_nef_poly_triangles( *root_N->p3, tess );
	BOOST_FOREACH( stl_triangle &t, triangles) {
		output  << "\n facet normal " << t.normal.x << " " << t.normal.y << " " << t.normal.z
			<< "\n    outer loop"
			<< "\n      vertex " << t.p1.x << " " << t.p1.y << " " << t.p1.z
			<< "\n      vertex " << t.p2.x << " " << t.p2.y << " " << t.p2.z
			<< "\n      vertex " << t.p3.x << " " << t.p3.y << " " << t.p3.z
			<< "\n    endloop"
			<< "\n  endfacet"
			<< "\n";
	}
	output << "endsolid OpenSCAD_Model\n";
}

/*!
Saves the current 3D CGAL Nef polyhedron volumes to the given stream,
which is assumed to already have had vertexes written to it.
CGAL Nef polyhedra Volumes are converted to AMF volumes.
*/
void export_amf_volumes(CGAL_Nef_polyhedron *root_N, std::ostream &output,
     std::map<stl_point,int> &vertexmap1, OpenSCAD::facetess::tesstype tess)
{
	// Volumes only work if we regularize the polyhedron first. See
	// https://github.com/noelwarr/rgal/blob/master/cpp/rb_Nef_polyhedron_3.cpp
	CGAL_Nef_polyhedron3::Volume_const_iterator vol_i;
	CGAL_Nef_polyhedron3 N = *(root_N->p3);
	CGAL_Nef_polyhedron3 reg_nef_poly = N.regularization();
 	CGAL_forall_volumes(vol_i,reg_nef_poly) {
		if ((*vol_i).mark()) { // use inner volumes, not outer volumes
			output << " <volume>\n";
			CGAL_Polyhedron P;
			// Use 'shell' visitor pattern.
			// convert_inner_shell_to_polyhedron is buggy.
			CGAL_Nef_polyhedron3::Shell_entry_const_iterator shell_i;
			CGAL_forall_shells_of(shell_i, vol_i) {
				CGAL_Nef_polyhedron3::SFace_const_handle sfch(shell_i);
				CGAL_Nef_polyhedron3 sub_nef_poly(reg_nef_poly,sfch);
				std::vector<stl_triangle> triangles = get_cgal_nef_poly_triangles( sub_nef_poly, tess );
				BOOST_FOREACH(stl_triangle &t, triangles) {
					output  << "  <triangle>\n"
						<< "   <v1>" << vertexmap1[ t.p1 ] << "</v1>\n"
						<< "   <v2>" << vertexmap1[ t.p2 ] << "</v2>\n"
						<< "   <v3>" << vertexmap1[ t.p3 ] << "</v3>\n"
						<< "  </triangle>\n";
				}
				output << " </volume>\n\n";
			} // forall shells
		} // inner volume
	} // forall volumes
}

/*!
Saves the current 3D CGAL Nef polyhedron as AMF to the given file.
The file must be open.
*/
void export_amf(CGAL_Nef_polyhedron *root_N, std::ostream &output, OpenSCAD::facetess::tesstype tess)
{
	// based on code originally by LogXen
	std::vector<stl_triangle> triangles = get_cgal_nef_poly_triangles( *root_N->p3, tess );

	std::map<stl_point,int> vertexmap1;
	std::map<int,stl_point> vertexmap2;
	int i = 0;
	BOOST_FOREACH( stl_triangle &t, triangles ) {
		if (! vertexmap1.count( t.p1 )) {
			vertexmap1[t.p1] = i;
			vertexmap2[i] = t.p1;
			i++;
		}
		if (! vertexmap1.count( t.p2 )) {
			vertexmap1[t.p2] = i;
			vertexmap2[i] = t.p2;
			i++;
		}
		if (! vertexmap1.count( t.p3 )) {
			vertexmap1[t.p3] = i;
			vertexmap2[i] = t.p3;
			i++;
		}
	}

	output  << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        	<< "<amf unit=\"millimeter\">\n"
        	<< " <object id=\"0\">\n"
        	<< " <mesh>\n"
		<< " <vertices>\n";
	std::map<int,stl_point>::const_iterator vi;
	for ( vi = vertexmap2.begin(); vi != vertexmap2.end(); ++vi ){
 		output  << "  <vertex><coordinates>\n"
			<< "   <x>" << vi->second.x << "</x>\n"
			<< "   <y>" << vi->second.y << "</y>\n"
			<< "   <z>" << vi->second.z << "</z>\n"
			<< "  </coordinates></vertex>\n";
	}
	output << " </vertices>\n\n";
	export_amf_volumes( root_N, output, vertexmap1, tess );
	output  << " </mesh>\n"
		<< " </object>\n"
		<< "</amf>\n";
}

void export_off(CGAL_Nef_polyhedron *root_N, std::ostream &output, OpenSCAD::facetess::tesstype tess)
{
	assert( root_N );
	assert( root_N->p3 );
	CGAL::Failure_behaviour old_behaviour = CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
	try {
		CGAL_Polyhedron P;
		nef3_to_polyhedron( *(root_N->p3), P, tess, tess );
		output << P;
	}
	catch (const CGAL::Assertion_exception &e) {
		PRINTB("CGAL error in CGAL_Nef_polyhedron3::convert_to_Polyhedron(): %s", e.what());
	}
	CGAL::set_error_behaviour(old_behaviour);
}

/*!
 Saves the current 2D CGAL Nef polyhedron as DXF to the given absolute filename.
 */
void export_dxf(CGAL_Nef_polyhedron *root_N, std::ostream &output)
{
	setlocale(LC_NUMERIC, "C"); // Ensure radix is . (not ,) in output
	// Some importers (e.g. Inkscape) needs a BLOCKS section to be present
	output << "  0\n"
				 <<	"SECTION\n"
				 <<	"  2\n"
				 <<	"BLOCKS\n"
				 <<	"  0\n"
				 << "ENDSEC\n"
				 << "  0\n"
				 << "SECTION\n"
				 << "  2\n"
				 << "ENTITIES\n";

	DxfData *dd =root_N->convertToDxfData();
	for (size_t i=0; i<dd->paths.size(); i++)
	{
		for (size_t j=1; j<dd->paths[i].indices.size(); j++) {
			const Vector2d &p1 = dd->points[dd->paths[i].indices[j-1]];
			const Vector2d &p2 = dd->points[dd->paths[i].indices[j]];
			double x1 = p1[0];
			double y1 = p1[1];
			double x2 = p2[0];
			double y2 = p2[1];
			output << "  0\n"
						 << "LINE\n";
			// Some importers (e.g. Inkscape) needs a layer to be specified
			output << "  8\n"
						 << "0\n"
						 << " 10\n"
						 << x1 << "\n"
						 << " 11\n"
						 << x2 << "\n"
						 << " 20\n"
						 << y1 << "\n"
						 << " 21\n"
						 << y2 << "\n";
		}
	}

	output << "  0\n"
				 << "ENDSEC\n";

	// Some importers (e.g. Inkscape) needs an OBJECTS section with a DICTIONARY entry
	output << "  0\n"
				 << "SECTION\n"
				 << "  2\n"
				 << "OBJECTS\n"
				 << "  0\n"
				 << "DICTIONARY\n"
				 << "  0\n"
				 << "ENDSEC\n";

	output << "  0\n"
				 <<"EOF\n";

	delete dd;
	setlocale(LC_NUMERIC, "");      // Set default locale
}

#endif // ENABLE_CGAL

#ifdef DEBUG
void export_stl(const PolySet &ps, std::ostream &output)
{
	output << "solid OpenSCAD_PolySet\n";
	BOOST_FOREACH(const PolySet::Polygon &p, ps.polygons) {
		output << "facet\n";
		output << "outer loop\n";
		BOOST_FOREACH(const Vector3d &v, p) {
			output << "vertex " << v[0] << " " << v[1] << " " << v[2] << "\n";
		}
		output << "endloop\n";
		output << "endfacet\n";
	}
	output << "endsolid OpenSCAD_PolySet\n";
}
#endif
