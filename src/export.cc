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

#define QUOTE(x__) # x__
#define QUOTED(x__) QUOTE(x__)

#include <boost/algorithm/string.hpp>
#include "export.h"
#include "printutils.h"
#include "polyset.h"
#include "dxfdata.h"
#include "importamfhandlers.h"
#ifdef ENABLE_CGAL
#include "CGAL_Nef_polyhedron.h"
#include "cgal.h"

std::string amf_dump_final(std::map<std::string,size_t> vertex_map, std::vector<amf_triangle> triangles)
{
	std::stringstream output;
	output << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\r\n"
		<< "<amf unit=\"millimeter\">\r\n"
		<< " <metadata type=\"cad\">OpenSCAD " << QUOTED(OPENSCAD_VERSION)
#ifdef OPENSCAD_COMMIT
		<< " (git " << QUOTED(OPENSCAD_COMMIT) << ")"
#endif
		<< "</metadata>\r\n"
		<< " <object id=\"0\">\r\n"
		<< "  <mesh>\r\n"
		<< "   <vertices>\r\n";
	size_t index = 0;
	std::map<std::string, size_t>::iterator i;
	for( i=vertex_map.begin(); i!=vertex_map.end(); i++ ) {
		std::string s = i->first;
		vertex_map[s] = index++;
		std::vector<std::string> coords;
		boost::split( coords, s, boost::is_any_of(" ") );
		assert( coords.size() == 3 ); // if not , somethings very wrong
		output << ""
			<< "    <vertex><coordinates>\r\n"
			<< "     <x>" << coords[0] << "</x>\r\n"
			<< "     <y>" << coords[1] << "</y>\r\n"
			<< "     <z>" << coords[2] << "</z>\r\n"
			<< "    </coordinates></vertex>\r\n";
	}
	output << "   </vertices>\r\n";
	output << "   <volume>\r\n";
	for(size_t i=0;i<triangles.size();i++) {
		amf_triangle t = triangles[i];
		output << "    <triangle>\r\n";
		size_t index, indexx;
		indexx = vertex_map[t.vs1];
		output << "     <v1>" << indexx << "</v1>\r\n";
		indexx = vertex_map[t.vs2];
		output << "     <v2>" << indexx << "</v2>\r\n";
		indexx = vertex_map[t.vs3];
		output << "     <v3>" << indexx << "</v3>\r\n";
		output << "    </triangle>\r\n";
	}
	output << "   </volume>\r\n"
		<< "  </mesh>\r\n"
		<< " </object>\r\n"
		<< "</amf>\r\n";
	return output.str();
}


std::string stl_dump_vertex( CGAL_Polyhedron::Vertex v[3], std::string vs[3])
{
	std::stringstream output;
	// If vertices are collinear, the unit normal is meaningless
	// so the default value of "1 0 0" can be used. If the vertices are not
	// collinear then the unit normal must be calculated from the
	// components.
	if (!CGAL::collinear(v[0].point(),v[1].point(),v[2].point())) {
		CGAL_Polyhedron::Traits::Vector_3 normal = CGAL::normal(v[0].point(),v[1].point(),v[2].point());
		output << "  facet normal "
			<< CGAL::sign(normal.x()) * sqrt(CGAL::to_double(normal.x()*normal.x()/normal.squared_length()))
			<< " "
			<< CGAL::sign(normal.y()) * sqrt(CGAL::to_double(normal.y()*normal.y()/normal.squared_length()))
			<< " "
			<< CGAL::sign(normal.z()) * sqrt(CGAL::to_double(normal.z()*normal.z()/normal.squared_length()))
			<< "\n";
	} else {
		output << "  facet normal 1 0 0\n";
	}
	output << "    outer loop\n"
		<< "      vertex " << vs[0] << "\n"
		<< "      vertex " << vs[1] << "\n"
		<< "      vertex " << vs[2] << "\n"
		<< "    endloop\n"
		<< "  endfacet\n";
	return output.str();
}

/* core export functionality shared between AMF and STL. */
std::string export_core(CGAL_Nef_polyhedron *root_N, std::string exptype)
{
	std::stringstream output;
	CGAL::Failure_behaviour old_behaviour = CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
	try {
	CGAL_Polyhedron P;
	root_N->p3->convert_to_Polyhedron(P);

	typedef CGAL_Polyhedron::Vertex Vertex;
	typedef CGAL_Polyhedron::Vertex_const_iterator  VCI;
	typedef CGAL_Polyhedron::Facet_const_iterator   FCI;
	typedef CGAL_Polyhedron::Halfedge_around_facet_const_circulator HFCC;

	setlocale(LC_NUMERIC, "C"); // Ensure radix is . (not ,) in output

	if (exptype == "stl") output << "solid OpenSCAD_Model\n";

	std::map<std::string,size_t> amf_vertex_map;
	std::vector<amf_triangle> amf_triangles;

	for (FCI fi = P.facets_begin(); fi != P.facets_end(); ++fi) {
		HFCC hc = fi->facet_begin();
		HFCC hc_end = hc;
		Vertex v[3];
		v[0] = *VCI((hc++)->vertex());
		v[2] = *VCI((hc++)->vertex());
		do {
			v[1] = v[2];
			v[2] = *VCI((hc++)->vertex());
			double x[3],y[3],z[3];
			x[0] = CGAL::to_double(v[0].point().x());
			y[0] = CGAL::to_double(v[0].point().y());
			z[0] = CGAL::to_double(v[0].point().z());
			x[1] = CGAL::to_double(v[1].point().x());
			y[1] = CGAL::to_double(v[1].point().y());
			z[1] = CGAL::to_double(v[1].point().z());
			x[2] = CGAL::to_double(v[2].point().x());
			y[2] = CGAL::to_double(v[2].point().y());
			z[2] = CGAL::to_double(v[2].point().z());
			std::string vs[3];
			for (int i=0;i<3;i++) {
				std::stringstream stream;
				stream << x[i] << " " << y[i] << " " << z[i];
				vs[i] = stream.str();
				if (exptype=="amf") amf_vertex_map[vs[i]]=0;
			}

			if (vs[0] != vs[1] && vs[0] != vs[2] && vs[1] != vs[2]) {
				// The above condition ensures that there are 3 distinct vertices
				if (exptype=="stl") {
					output << stl_dump_vertex( v, vs );
				} else if (exptype=="amf") {
					amf_triangle tri = {vs[0], vs[1], vs[2]};
					amf_triangles.push_back(tri);
				}
			}
		} while (hc != hc_end);
	}

	if (exptype=="amf")
		output << amf_dump_final( amf_vertex_map, amf_triangles );
	else if (exptype=="stl")
		output << "endsolid OpenSCAD_Model\n";

	setlocale(LC_NUMERIC, "");  // Set default locale

	} catch (CGAL::Assertion_exception e) {
		PRINTB("CGAL error in CGAL_Nef_polyhedron3::convert_to_Polyhedron(): %s", e.what());
	}
	CGAL::set_error_behaviour(old_behaviour);
	return output.str();
}

/*!
	Saves the current 3D CGAL Nef polyhedron as STL to the given file.
	The file must be open.
 */
void export_stl(CGAL_Nef_polyhedron *root_N, std::ostream &output)
{
	output << export_core( root_N, "stl");
}

/*!
Saves the current 3D CGAL Nef polyhedron as AMF to the given file.
The file must be open.
 */
void export_amf(CGAL_Nef_polyhedron *root_N, std::ostream &output)
{
	output << export_core( root_N, "amf" );

}

void export_off(CGAL_Nef_polyhedron *root_N, std::ostream &output)
{
	CGAL::Failure_behaviour old_behaviour = CGAL::set_error_behaviour(CGAL::THROW_EXCEPTION);
	try {
		CGAL_Polyhedron P;
		root_N->p3->convert_to_Polyhedron(P);
		output << P;
	}
	catch (CGAL::Assertion_exception e) {
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
	setlocale(LC_NUMERIC, "");  // Set default locale
}

#endif

#ifdef DEBUG
#include <boost/foreach.hpp>
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
