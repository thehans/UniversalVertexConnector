// Universal Vertex Connector
// by Hans Loeblich
// License: CC-BY-SA (Creative Commons - Attribution - Share Alike license)

use <conway.scad> // https://github.com/KitWallace/Polyhedra/blob/master/openSCAD/conway.scad
use <generated_archimedean_solids.scad>
use <generated_johnson_solids.scad>

/* [View/Export Options] */
view_type = 1; // [0:Single Vertex - render/export the given vertex_index, 1:All Vertices - preview mockup, 2:Plain Polyhedron - no vertex connectors]
// Place polyhedron on largest face
place_face = false;
// Display edges in "All Vertices" view
show_edges = true;
// Vertex index to display / export (for Single vertex view only).
vertex_index = 0;

/* [Polyhedron Selection] */
// Only one of the select boxes/sliders below is enabled at a time based on this category:
Polyhedron_Category = 1; // [1:PLATONIC SOLID,2:ARCHIMEDEAN SOLID,3:JOHNSON SOLID,4:PRISM,5:ANTIPRISM]
Select_Platonic_Solid = 1; // [1:Tetrahedron,2:Cube,3:Octahedron,4:Dodecahedron,5:Icosahedron]
Select_Archimedean_Solid = 1; // [1:Cuboctahedron,2:Icosidodecahedron,3:L Snub Cube,4:Rhombicosidodecahedron,5:Rhombicuboctahedron,6:R Snub Dodecahedron,7:Truncated Cube,8:Truncated Cuboctahedron,9:Truncated Dodecahedron,10:Truncated Icosahedron,11:Truncated Icosidodecahedron,12:Truncated Octahedron,13:Truncated Tetrahedron]
Select_Johnson_Solid = 1; // [1:J1 Square pyramid,2:J2 Pentagonal pyramid,3:J3 Triangular cupola,4:J4 Square cupola,5:J5 Pentagonal cupola,6:J6 Pentagonal rotunda,7:J7 Elongated triangular pyramid,8:J8 Elongated square pyramid,9:J9 Elongated pentagonal pyramid,10:J10 Gyroelongated square pyramid,11:J11 Gyroelongated pentagonal pyramid,12:J12 Triangular bipyramid,13:J13 Pentagonal bipyramid,14:J14 Elongated triangular bipyramid,15:J15 Elongated square bipyramid,16:J16 Elongated pentagonal bipyramid,17:J17 Gyroelongated square bipyramid,18:J18 Elongated triangular cupola,19:J19 Elongated square cupola,20:J20 Elongated pentagonal cupola,21:J21 Elongated pentagonal rotunda,22:J22 Gyroelongated triangular cupola,23:J23 Gyroelongated square cupola,24:J24 Gyroelongated pentagonal cupola,25:J25 Gyroelongated pentagonal rotunda,26:J26 Gyrobifastigium,27:J27 Triangular orthobicupola,28:J28 Square orthobicupola,29:J29 Square gyrobicupola,30:J30 Pentagonal orthobicupola,31:J31 Pentagonal gyrobicupola,32:J32 Pentagonal orthocupolarotunda,33:J33 Pentagonal gyrocupolarotunda,34:J34 Pentagonal orthobirotunda,35:J35 Elongated triangular orthobicupola,36:J36 Elongated triangular gyrobicupola,37:J37 Elongated square gyrobicupola,38:J38 Elongated pentagonal orthobicupola,39:J39 Elongated pentagonal gyrobicupola,40:J40 Elongated pentagonal orthocupolarotunda,41:J41 Elongated pentagonal gyrocupolarotunda,42:J42 Elongated pentagonal orthobirotunda,43:J43 Elongated pentagonal gyrobirotunda,44:J44 Gyroelongated triangular bicupola,45:J45 Gyroelongated square bicupola,46:J46 Gyroelongated pentagonal bicupola,47:J47 Gyroelongated pentagonal cupolarotunda,48:J48 Gyroelongated pentagonal birotunda,49:J49 Augmented triangular prism,50:J50 Biaugmented triangular prism,51:J51 Triaugmented triangular prism,52:J52 Augmented pentagonal prism,53:J53 Biaugmented pentagonal prism,54:J54 Augmented hexagonal prism,55:J55 Parabiaugmented hexagonal prism,56:J56 Metabiaugmented hexagonal prism,57:J57 Triaugmented hexagonal prism,58:J58 Augmented dodecahedron,59:J59 Parabiaugmented dodecahedron,60:J60 Metabiaugmented dodecahedron,61:J61 Triaugmented dodecahedron,62:J62 Metabidiminished icosahedron,63:J63 Tridiminished icosahedron,64:J64 Augmented tridiminished icosahedron,65:J65 Augmented truncated tetrahedron,66:J66 Augmented truncated cube,67:J67 Biaugmented truncated cube,68:J68 Augmented truncated dodecahedron,69:J69 Parabiaugmented truncated dodecahedron,70:J70 Metabiaugmented truncated dodecahedron,71:J71 Triaugmented truncated dodecahedron,72:J72 Gyrate rhombicosidodecahedron,73:J73 Parabigyrate rhombicosidodecahedron,74:J74 Metabigyrate rhombicosidodecahedron,75:J75 Trigyrate rhombicosidodecahedron,76:J76 Diminished rhombicosidodecahedron,77:J77 Paragyrate diminished rhombicosidodecahedron,78:J78 Metagyrate diminished rhombicosidodecahedron,79:J79 Bigyrate diminished rhombicosidodecahedron,80:J80 Parabidiminished rhombicosidodecahedron,81:J81 Metabidiminished rhombicosidodecahedron,82:J82 Gyrate bidiminished rhombicosidodecahedron,83:J83 Tridiminished rhombicosidodecahedron,84:J84 Snub disphenoid,85:J85 Snub square antiprism,86:J86 Sphenocorona,87:J87 Augmented sphenocorona,88:J88 Sphenomegacorona,89:J89 Hebesphenomegacorona,90:J90 Disphenocingulum,91:J91 Bilunabirotunda,92:J92 Triangular hebesphenorotunda]
// Select "Prism" or "Antiprism" Polyhedron Category First!
Sides = 3; // [3:100]

/* [Connector settings] */
connector_type = 1; // [0:female, 1:male]
// Outer Diameter of edge material
OD = 11.8;
// Thickness.  ID is calculated as OD-2*Th (female: thickness of printed walls, male: thickness of edge material/pipe wall) 
Th = 0.3;
// Length of pipes used (assumes all equal length when rendering edges)
edge_length = 224;
// The "leg base" spaces out edge around the vertex to keep the ends from intersecting.  leg_base_length is calculated as OD * leg_base_scale.  Minimum value with no intersection depends on the polyhedron's smallest face(triangle: 0.87, square: 0.50, pentagon: 0.38)
leg_base_scale = 0.87; // [0.2:0.01:1.0]
// Max leg length
leg_insert_length = 10.0;
// Minimum insert length when flattening
min_insert_length = 5;
// Options for flattening a side of the connector to lay on a 3D printer bed.
flatten_type = 2; // [0:none, 1:peak, 2:legs]
// Options for calculating a "vertex normal" that is used as plane normal for the flatten operation, and for orientation in single vertex view
normal_type = 1; // [0: simple average, 1: 3pt plane, 2: best fit plane *WIP* *BROKEN*]

/* [Hidden] */
$fn = 60; // helps to use something with 3,4,5 as factors in my half-baked theory


// CONSTANT ENUMS: DO NOT CHANGE
// View type enumeration
VIEW_SINGLE = 0;
VIEW_ALL = 1;
VIEW_PLAIN = 2;

// connector_type enumeration
TYPE_FEMALE = 0;
TYPE_MALE = 1;

// index keys for edge object properties
EDGE_VECTOR = 0; // unit direction vector
EDGE_LENGTH = 1; // edge length
EDGE_DELTA_P = 2; // difference between points
EDGE_P_I1 = 3; // point1 index
EDGE_P_I2 = 4; // point2 index

// flatten options for module flatten_vertex(r, vn, a, flatten 
FLATTEN_NONE = 0;
FLATTEN_PEAK = 1;
FLATTEN_LEGS = 2;

// Options for calculating a "vertex normal" for orienting the object, also affects the flatten_vertex operation
NORMAL_SIMPLE_AVG = 0; // works well only on very symmetrical vertices
NORMAL_3PT_PLANE = 1; // perfect for vertices with 3 edges
NORMAL_PLANE_FIT = 2; // point-based plane approximation for vertices with edges > 3

POLYHEDRON_PLATONIC = 1;
POLYHEDRON_ARCHIMEDEAN = 2;
POLYHEDRON_JOHNSON = 3;
POLYHEDRON_PRISM = 4;
POLYHEDRON_ANTIPRISM = 5;


// Calculated Values
ID = OD - 2 * Th;
leg_base_length = OD*leg_base_scale; // if female type then use ID?


// Start Script
solid = SelectPolyhedron(Polyhedron_Category, Select_Platonic_Solid, Select_Archimedean_Solid, Select_Johnson_Solid, Sides);
p = (place_face ? place(solid) : solid);
poly = [p_vertices(p),p_faces(p)];

if (view_type == VIEW_SINGLE) {
  points = poly[0];
  npoints = len(points);
  if (vertex_index >= 0 && vertex_index < npoints) {
    vertex_connector_from_poly(poly, vertex_index, leg_base_length, leg_insert_length, OD, Th, connector_type, flatten_type);
  } else {
    color("red") {
      text(str("vertex_index " , vertex_index, " out of range"), halign="center");
      translate([0,-20]) text(str("valid range is 0-", npoints-1," (", p_name(p), ")"), halign="center");
    }
  }
} else if (view_type == VIEW_ALL) {
  echo(poly);
  // mockup entire structure of vertex connectors, optionally rendering edges
  show_all_connectors(poly, leg_base_length, leg_insert_length, OD, Th, connector_type, flatten_type, show_edges);
} else if (view_type == VIEW_PLAIN) {
  points = poly[0];
  faces = poly[1];
  face0 = faces[0];
  poly_edge_scale = 1/norm(points[face0[0]]-points[face0[1]]); // assumes all edges are same length, based on first edge in list
  scale((edge_length+2*leg_base_length)* poly_edge_scale) polyhedron(points, faces); // preview small r=1? polyhedron
}

function SelectPolyhedron(Polyhedron_Category, Platonic_Solid, Archimedean_Solid, Johnson_Solid, Sides) =
  (Polyhedron_Category == POLYHEDRON_PLATONIC ? 
    PlatonicSolid(Platonic_Solid) :
    (Polyhedron_Category == POLYHEDRON_ARCHIMEDEAN ?
      ArchimedeanSolid(Archimedean_Solid) :
      (Polyhedron_Category == POLYHEDRON_JOHNSON ? 
        JohnsonSolid(Johnson_Solid) :
        (Polyhedron_Category == POLYHEDRON_PRISM ? 
          P(Sides) : 
          (Polyhedron_Category == POLYHEDRON_ANTIPRISM ? 
            A(Sides) :
            undef
          )
        )
      )
    )
  );

// Choose from platonic solids, enumerated for n in range [1:5]
function PlatonicSolid(n) = [T(),C(),O(),D(),I()][n-1];

// Universal vertex connector functions and modules 
// v: direction vector
// l: length
// dp: difference betwen points (points[p_i1]-points[p_i2])
// p_i1, p_i2: point index
function edge(v, l, dp, p_i1, p_i2) = [v, l, dp, p_i1, p_i2];

function get_edges(poly) = 
  let(
    points = poly[0], 
    faces = poly[1]
  )
  [ for(face = faces) let(n = len(face)) for(i = [0:1:n-1]) 
      let(
        p_i1 = face[i % n], 
        p_i2 = face[(i + 1) % n],
        p1 = points[p_i1],
        p2 = points[p_i2],
        dp = p2 - p1,
        l = norm(dp),
        v = dp/l
      ) 
      edge(v, l, dp, p_i1, p_i2)
  ];

// given the edges of a vertex, return a weighted sum unit normal
function vertex_normal(edges) = 
  (normal_type==NORMAL_SIMPLE_AVG ? simple_avg(edges) : 
    (normal_type==NORMAL_3PT_PLANE ? plane_normal(edges) :
      (normal_type==NORMAL_PLANE_FIT ? plane_normal_from_points([for (edge=edges) edge[EDGE_VECTOR]]) :
        undef
      ) 
    )
  );

function simple_avg(edges) = -normalize( sumvk(edges, EDGE_VECTOR, len(edges) - 1) );

function poly_face_normal(poly, face_index) = 
  let(
    points = poly[0],
    faces = poly[1],
    face = faces[face_index],
    edge1 = points[face[1]]-points[face[0]],
    edge2 = points[face[2]]-points[face[1]]
  ) 
  normalize(cross(edge1, edge2));



// Only uses first 3 edges, treating leg endpoints as points defining a plane
// not sure about optimal solutions for >3 legs and (lack of)symmetry such that averaging doesn't work well
function plane_normal(edges) = 
  let(
    P0 = edges[0][EDGE_DELTA_P],
    P1 = edges[1][EDGE_DELTA_P],
    P2 = edges[2][EDGE_DELTA_P],
    DP01 = P1-P0,
    DP02 = P2-P0,
    n = cross(DP02, DP01)
  ) //echo(edges)
  normalize(n) * (n*P0 < 0 ? 1 : -1);

// Fit plane to points:
function plane_normal_from_points(points) = 
  let(
    l = len(points),
    psum = sumv(points, l-1),
    centroid = psum / l,
    l_inv = 1/l,
    xs = [for (p=points) p.x-centroid.x],
    ys = [for (p=points) p.y-centroid.y],
    zs = [for (p=points) p.z-centroid.z],
    xx = xs*xs*l_inv,
    xy = xs*ys*l_inv,
    xz = xs*zs*l_inv,
    yy = ys*ys*l_inv,
    yz = ys*zs*l_inv,
    zz = zs*zs*l_inv,
    det_x = yy*zz - yz*yz,
    det_y = xx*zz - xz*xz,
    det_z = xx*yy - xy*xy,
    xaxis_dir = [        det_x, xz*yz - xy*zz, xy*yz - xz*yy],
    yaxis_dir = [xz*yz - xy*zz,         det_y, xy*xz - yz*xx],
    zaxis_dir = [xy*yz - xy*yy, xy*xz - yz*xx,         det_z],
    xweight = det_x*det_x,
    weighted_dir1 = xweight*xaxis_dir,
    yweight = det_y*det_y * ((weighted_dir1*yaxis_dir < 0) ? -1 : 1),
    weighted_dir2 = weighted_dir1 + yweight*yaxis_dir,
    zweight = det_z*det_z * ((weighted_dir2*zaxis_dir < 0) ? -1 : 1),
    normal = weighted_dir2 + zweight*zaxis_dir
  )
  -normalize(normal);


// sum a vector
function sumv(v,i,s=0) = (i==s ? v[i] : v[i] + sumv(v,i-1,s));

// sum the key element k from all vectors in v (a vector of vectors), from index s to i
function sumvk(v,k,i,s=0) = (i==s ? v[i][k] : v[i][k] + sumvk(v,k,i-1,s));

function normalize(v) = v/norm(v);

// return a(any) unit vector which is orthogonal to the input
function get_orthogonal_vector(v) = v.x == 0 && v.y == 0 ? [0,1,0] : normalize([-v.y, v.x, 0]);

function angle_between_vectors(v1,v2) = acos((v1*v2)/(norm(v1)*norm(v2)));

// rotate from vector v1 to v2
module rotatev(v1, v2) {
    //echo("v1", v1, "v2", v2);
    a = angle_between_vectors(v1, v2);
    v = (a == 180) ? get_orthogonal_vector(v1) : cross(v1, v2);
    //echo("a", a, "v", v);
    rotate(a=a, v=v) children();
}

module show_all_connectors(poly, leg_base_length, leg_insert_length, OD, Th, connector_type, flatten_type, show_edges) {
  points = poly[0];
  distance = edge_length + leg_base_length*2;
  edges = get_edges(poly);
  poly_edge_size = edges[0][EDGE_LENGTH]; // assumes all edges are same length, based on first edge in list
  for (vertex_index = [0:1:len(points)-1]) {
    filtered_edges = [for(edge = edges) if(edge[EDGE_P_I1] == vertex_index) edge];
    translate(points[vertex_index]*distance/poly_edge_size) {
      vertex_connector(filtered_edges, leg_base_length, leg_insert_length, OD, Th, connector_type, flatten_type, show_edges);
      if (show_edges) {
        for (edge = filtered_edges) {
          v2 = edge[EDGE_VECTOR];
          // edges based on faces are duplicated, so only draw one of the pair
          if (edge[EDGE_P_I1] < edge[EDGE_P_I2]) 
            rotatev([0,0,1],v2) 
              translate([0,0,leg_base_length]) 
                %cylinder(d=OD, h=edge_length);
        }
      }
    }
  }

}

// for manually exporting single specific vertex connector (for uniform polyhedra vertex_index=0 is all that is needed)
// orient with vertex normal pointing up
module vertex_connector_from_poly(poly, vertex_index, leg_base_length, leg_insert_length, OD, Th, connector_type, flatten_type=0) {
  edges = get_edges(poly);
  filtered_edges = [for(edge = edges) if(edge[EDGE_P_I1] == vertex_index) edge];
  vn = vertex_normal(filtered_edges);
  // vertex connector
  rotatev(vn, [0,0,1]) 
    vertex_connector(filtered_edges, leg_base_length, leg_insert_length, OD, Th, connector_type, flatten_type, show_edges);
}

// create flat spot on vertex connector to ease printing
// putting conditional here saves a lot of copy/paste inside vertex_connector
module flatten_vertex(r, vn, a, flatten=FLATTEN_LEGS) {
  leg_total_min = leg_base_length + min_insert_length;
  clip_offset = leg_total_min * cos(a) + r*sin(a);
  //echo("r", r, "vn", vn, "a", a, "leg_base_length", leg_base_length, "min_insert_length", min_insert_length, 
  //      "leg_total_min", leg_total_min,  "clip_offset", clip_offset, "norm(vn)", norm(vn));
  s = leg_base_length + leg_insert_length + 2*r;
  if (flatten == FLATTEN_PEAK) {
    difference() {
      children();
      rotatev([0,0,1], vn)
        translate([-s,-s,r/sqrt(3)]) cube(2*s);
    }
  } else if (flatten == FLATTEN_LEGS) {
    difference() {
      children();
      rotatev([0,0,-1], vn)
        translate([-s,-s,clip_offset]) cube(2*s);
    }
  } else {
    children();
  }
}




module vertex_connector(edges, leg_base_length, leg_insert_length, OD, Th, connector_type, flatten_type=FLATTEN_LEGS, show_edges=false) {
  //echo(edges);
  v1 = [0,0,1];
  vn = vertex_normal(edges);
  angles = [for(edge = edges) angle_between_vectors(-vn, edge[EDGE_VECTOR]) ];
  minangle = min(angles);
  maxangle = max(angles);
  //echo("vn", vn, "minangle", minangle, "maxangle", maxangle);
  
  if (connector_type == TYPE_FEMALE) {
    OD1 = OD+2*Th; // OD of printed part
    ID1 = OD; // ID of printed part
    flatten_vertex(ID1/2, vn, minangle, flatten_type) difference() {
      union() {
        sphere(d=OD1);
        for (edge = edges) {
          v2 = edge[EDGE_VECTOR];
          rotatev(v1,v2) 
            cylinder(d=OD1, h=leg_base_length+leg_insert_length);
        }
      }

      union() for (edge = edges) {
        v2 = edge[EDGE_VECTOR];
        rotatev(v1,v2) translate([0,0,leg_base_length]) 
          cylinder(d=ID1, h=leg_insert_length+1);
      }
    }
  } else if (connector_type == TYPE_MALE) {
    flatten_vertex(ID/2, vn, minangle, flatten_type) union() {
      sphere(d=OD);
      for (edge = edges) {
        v2 = edge[EDGE_VECTOR];
        rotatev(v1,v2) {
          y0 = 0;
          y1 = leg_base_length;
          y2 = y1 + leg_insert_length;
          points = [
            [0,y0], [OD/2, y0], 
            [OD/2, y1], [ID/2, y1],
            [ID/2, y2], [0, y2]
          ];
          rotate_extrude() polygon(points);
        }
      }
    } // else if (type == TYPE_) {}

    // TODO both? (outer wall + inner insert) 
    // TODO ball end type (for molecular models?) (female socket, also possibly "both" style)
  }
}
