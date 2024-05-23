// Gmsh project created on Wed May 24 14:41:07 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 1, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {1, 0, 0, 1.0};
//+
Line(1) = {2, 3};
//+
Line(2) = {3, 4};
//+
Line(3) = {1, 1};
//+
Line(3) = {1, 4};
//+
Line(4) = {1, 3};
//+
Line(5) = {2, 1};
//+
Curve Loop(1) = {4, -1, 5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {2, -3, 4};
//+
Plane Surface(2) = {2};
//+
Transfinite Curve {5} = 20 Using Progression 1;
//+
Transfinite Curve {2} = 20 Using Progression 1;
//+
Transfinite Curve {1} = 20 Using Progression 1;
//+
Transfinite Curve {3} = 20 Using Progression 1;
//+
Transfinite Curve {4} = 40 Using Progression 1;
//+
Transfinite Curve {4} = 14 Using Progression 1;
//+
Transfinite Curve {4} = 14 Using Progression 1;
//+
Transfinite Curve {4} = 14 Using Progression 1;
//+
Transfinite Curve {4, 4} = 14 Using Progression 1;
//+
Transfinite Curve {4} = 10 Using Progression 1;
//+
Transfinite Curve {4} = 10 Using Progression 1;
//+
Transfinite Curve {4} = 10 Using Progression 1;
//+
Transfinite Curve {4} = 10 Using Progression 1;
//+
Transfinite Curve {4} = 10 Using Progression 1;
//+
Transfinite Curve {4} = 10 Using Progression 1;
//+
Transfinite Curve {4} = 15 Using Progression 1;
//+
Transfinite Curve {4} = 15 Using Progression 1;
//+
Transfinite Curve {4} = 15 Using Progression 1;
//+
Transfinite Curve {4} = 20 Using Progression 1;
//+
Transfinite Curve {4} = 30 Using Progression 1;
//+
Transfinite Curve {5} = 10 Using Progression 1;
//+
Transfinite Curve {1} = 10 Using Progression 1;
//+
Transfinite Curve {2} = 10 Using Progression 1;
//+
Transfinite Curve {3} = 10 Using Progression 1;
//+
Transfinite Curve {5} = 10 Using Progression 0.8;
//+
Transfinite Curve {2} = 10 Using Progression 0.8;
//+
Transfinite Curve {3} = 10 Using Progression 1.25;
//+
Transfinite Curve {1} = 10 Using Progression 1.25;
//+
Transfinite Curve {1, 1} = 10 Using Progression 0.8;
//+
Transfinite Curve {2} = 10 Using Progression 1.25;
