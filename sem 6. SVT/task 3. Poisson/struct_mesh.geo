// Gmsh project created on Wed May 24 00:19:21 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 1, 0, 1.0};
//+
Point(2) = {1, 1, 0, 1.0};
//+
Point(3) = {0, 0, 0, 1.0};
//+
Point(4) = {1, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 4};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Transfinite Curve {4} = 10 Using Progression 1;
//+
Transfinite Curve {1} = 10 Using Progression 1;
//+
Transfinite Curve {2} = 10 Using Progression 1;
//+
Transfinite Curve {3} = 10 Using Progression 1;
//+
Transfinite Curve {4} = 8 Using Progression 1;
//+
Transfinite Curve {1} = 8 Using Progression 1;
//+
Transfinite Curve {2} = 8 Using Progression 1;
//+
Transfinite Curve {3} = 8 Using Progression 1;
//+
Transfinite Curve {4} = 5 Using Progression 1;
//+
Transfinite Curve {3} = 5 Using Progression 1;
//+
Transfinite Curve {2} = 5 Using Progression 1;
//+
Transfinite Curve {1} = 5 Using Progression 1;
//+
Transfinite Curve {2} = 5 Using Progression 1;
//+
Transfinite Curve {4} = 5 Using Progression 1;
//+
Transfinite Curve {4} = 6 Using Progression 1;
//+
Transfinite Curve {2} = 6 Using Progression 1;
//+
Transfinite Curve {3} = 6 Using Progression 1;
//+
Transfinite Curve {1} = 6 Using Progression 1;
//+
Transfinite Curve {4} = 10 Using Progression 1;
//+
Transfinite Curve {2} = 10 Using Progression 1;
//+
Transfinite Curve {1} = 10 Using Progression 1;
//+
Transfinite Curve {3} = 10 Using Progression 1;
//+
Transfinite Curve {4} = 40 Using Progression 1;
//+
Transfinite Curve {1} = 40 Using Progression 1;
//+
Transfinite Curve {2} = 40 Using Progression 1;
//+
Transfinite Curve {3} = 40 Using Progression 1;
//+
Transfinite Curve {4} = 2 Using Progression 1;
//+
Transfinite Curve {1} = 2 Using Progression 1;
//+
Transfinite Curve {2} = 2 Using Progression 1;
//+
Transfinite Curve {3} = 2 Using Progression 1;
