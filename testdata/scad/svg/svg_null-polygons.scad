linear_extrude(height=1) import_svg("../../svg/null-polygons.svg");
translate([0,20,0]) linear_extrude("../../svg/null-polygons.svg", height=1);
