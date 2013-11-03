#ifndef EXPORT_H_
#define EXPORT_H_

#include <iostream>
#include "Tree.h"
#include "Camera.h"
#ifdef ENABLE_CGAL
#include "facetess.h"
#include "CGAL_Nef_polyhedron.h"

using namespace OpenSCAD::facetess;
void export_stl(CGAL_Nef_polyhedron *root_N, std::ostream &output, tesstype tess);
void export_amf(CGAL_Nef_polyhedron *root_N, std::ostream &output, tesstype tess);
void export_off(CGAL_Nef_polyhedron *root_N, std::ostream &output, tesstype tess);
void export_dxf(CGAL_Nef_polyhedron *root_N, std::ostream &output);
void export_png_with_cgal(CGAL_Nef_polyhedron *root_N, Camera &c, std::ostream &output);
void export_png_with_opencsg(Tree &tree, Camera &c, std::ostream &output);
void export_png_with_throwntogether(Tree &tree, Camera &c, std::ostream &output);

#endif

#ifdef DEBUG
void export_stl(const class PolySet &ps, std::ostream &output);
#endif

#endif
