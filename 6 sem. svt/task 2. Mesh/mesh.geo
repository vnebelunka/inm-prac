//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 1, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {1, 0, 0, 1.0};
//+
Line(1) = {2, 1};
//+
Line(2) = {2, 3};
//+
Line(3) = {4, 3};
//+
Line(4) = {1, 4};
//+
Curve Loop(1) = {4, 3, -2, 1};
//+
Plane Surface(1) = {1};
//+
Transfinite Curve {1} = 10 Using Progression 1;
//+
Line(5) = {2, 4};
//+
Physical Surface(1) -= {1};
//+
Physical Surface(1) -= {1};
//+
Physical Surface(1) -= {1};
//+
Physical Surface(1) -= {1};
//+
Physical Surface(1) -= {1};
//+
Physical Surface(1) -= {1};
//+
Physical Surface(1) -= {1};
//+
Physical Surface(1) -= {1};
//+
Physical Surface(1) -= {1};
//+
Physical Surface(1) -= {1};
//+
Physical Surface(1) -= {1};
//+
Physical Surface(1) -= {1};
//+
Physical Surface(1) -= {1};
//+
Physical Surface(1) -= {1};
//+
Physical Surface(1) -= {1};
//+
Physical Surface(1) -= {1};
//+
Physical Surface(1) -= {1};
//+
Physical Surface(1) -= {1};
//+
Physical Surface(1) -= {1};
//+
Physical Surface(1) -= {1};
//+
Physical Surface(1) -= {1};
//+
Curve Loop(2) = {1, 4, -5};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {2, -3, -5};
//+
Plane Surface(3) = {3};
//+
Transfinite Curve {5} = 10 Using Progression 1;
//+
Transfinite Curve {1, 1} = 2 Using Progression 1;
//+
Transfinite Curve {1} = 10 Using Progression 2;
//+
Transfinite Curve {1} = 10 Using Progression 1.5;
//+
Transfinite Curve {4} = 10 Using Progression 1.5;
//+
Transfinite Curve {3} = 10 Using Progression 1.5;
//+
Transfinite Curve {2} = 10 Using Progression 1.5;
//+
Transfinite Curve {4, 4} = 10 Using Progression 0.5;
//+
Transfinite Curve {5} = 10 Using Progression 1.5;
//+
Transfinite Curve {5} = 10 Using Progression 1.5;
//+
Point(5) = {0.5, 0.5, 0, 1.0};
//+
Line(6) = {2, 5};
//+
Line(7) = {4, 5};
//+
Curve Loop(4) = {1, 4, 7, -6};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {2, -3, 7, -6};
//+
Plane Surface(5) = {5};
//+
Transfinite Curve {1} = 10 Using Progression 1.5;
//+
Transfinite Curve {6} = 10 Using Progression 1.5;
//+
Transfinite Curve {7} = 10 Using Progression 1.5;
//+
Transfinite Curve {4} = 10 Using Progression 1.5;
//+
Transfinite Curve {3} = 10 Using Progression 1.5;
//+
Transfinite Curve {2} = 10 Using Progression 1.5;
//+
Transfinite Curve {4} = 10 Using Progression 0.5;
//+
Plane Surface(6) = {2};
//+
Plane Surface(7) = {3};
//+
Plane Surface(8) = {3};
//+
Recursive Delete {
  Curve{5}; 
}
//+
Recursive Delete {
  Curve{5}; 
}
//+
Line(8) = {2, 5};
//+
Line(9) = {5, 4};
//+
Transfinite Curve {1} = 10 Using Progression 1.5;
//+
Transfinite Curve {2} = 10 Using Progression 1.5;
//+
Transfinite Curve {3} = 10 Using Progression 1.5;
//+
Transfinite Curve {4} = 10 Using Progression 1.5;
//+
Transfinite Curve {9} = 10 Using Progression 1.5;
//+
Transfinite Curve {8} = 10 Using Progression 1.5;
//+
Curve Loop(6) = {1, 4, -9, -8};
//+
Plane Surface(9) = {6};
//+
Curve Loop(7) = {3, -2, 8, 9};
//+
Plane Surface(10) = {7};
//+
Transfinite Curve {4} = 10 Using Progression 0.5;
//+
Transfinite Curve {9} = 10 Using Progression 0.5;
//+
Transfinite Curve {9} = 10 Using Progression 1.5;
//+
Transfinite Curve {9} = 10 Using Progression 1.5;
//+
Transfinite Curve {9} = 10 Using Progression 1.5;
//+
Transfinite Curve {9} = 10 Using Progression 0.5;
