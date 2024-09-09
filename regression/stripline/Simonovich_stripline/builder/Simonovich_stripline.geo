// Control
//   build=Simonovich_stripline
//   check_limits=true
// EndControl
// Strip
//    name=Simonovich_stripline
//    use_symmetry=false
//    upper_material=prepreg
//    upper_thickness=0.000301
//    lower_material=core
//    lower_thickness=0.0003048
//    left_side_gap=0.00486
//    right_side_gap=0.00486
//    trace_thickness=3.175e-05
//    trace_width=0.0002794
//    trace_etch_angle=60
//    default_conductor_material=PEC
//    trace_material_bottom=copper_core
//    trace_material_top=copper_prepreg
//    trace_material_sides=copper_prepreg
//    upper_groundplane_material=copper_prepreg
//    lower_groundplane_material=copper_core
//    length=0.001
// EndStrip
Point(1) = {-0.0049997,0,0,0.00015145};
Point(2) = {-0.0004445,0,0,0.00015145};
Point(3) = {-0.0001397,0,0,0.0003048};
Point(4) = {0,0,0,0.0003048};
Point(5) = {0.0001397,0,0,0.0003048};
Point(6) = {0.0004445,0,0,0.00015145};
Point(7) = {0.0049997,0,0,0.00015145};
Point(8) = {-0.0049997,0.0003048,0,0.00015145};
Point(9) = {-0.0002032,0.0003048,0,0.00015145};
Point(10) = {-0.00017145,0.0003048,0,3.175e-05};
Point(11) = {-0.0001397,0.0003048,0,3.175e-05};
Point(12) = {-0.00010795,0.0003048,0,3.175e-05};
Point(13) = {-7.62e-05,0.0003048,0,0.0003048};
Point(14) = {0,0.0003048,0,0.0003048};
Point(15) = {7.62e-05,0.0003048,0,0.0003048};
Point(16) = {0.00010795,0.0003048,0,3.175e-05};
Point(17) = {0.0001397,0.0003048,0,3.175e-05};
Point(18) = {0.00017145,0.0003048,0,3.175e-05};
Point(19) = {0.0002032,0.0003048,0,0.00015145};
Point(20) = {0.0049997,0.0003048,0,0.00015145};
Point(21) = {-0.000121369128953229,0.00033655,0,3.175e-05};
Point(22) = {-8.96191289532294e-05,0.00033655,0,3.175e-05};
Point(23) = {-5.78691289532294e-05,0.00033655,0,0.0003048};
Point(24) = {0,0.00033655,0,0.0003048};
Point(25) = {5.78691289532294e-05,0.00033655,0,0.0003048};
Point(26) = {8.96191289532294e-05,0.00033655,0,3.175e-05};
Point(27) = {0.000121369128953229,0.00033655,0,3.175e-05};
Point(28) = {-0.0049997,0.0006058,0,0.00015145};
Point(29) = {0,0.0006058,0,0.00015145};
Point(30) = {0.0049997,0.0006058,0,0.00015145};
Point(31) = {-0.0049997,0,0.001,0.00015145};
Point(32) = {-0.0004445,0,0.001,0.00015145};
Point(33) = {-0.0001397,0,0.001,0.0003048};
Point(34) = {0,0,0.001,0.0003048};
Point(35) = {0.0001397,0,0.001,0.0003048};
Point(36) = {0.0004445,0,0.001,0.00015145};
Point(37) = {0.0049997,0,0.001,0.00015145};
Point(38) = {0.0049997,0.0003048,0.001,0.00015145};
Point(39) = {0.0002032,0.0003048,0.001,0.00015145};
Point(40) = {0.00017145,0.0003048,0.001,3.175e-05};
Point(41) = {0.0001397,0.0003048,0.001,3.175e-05};
Point(42) = {0.00010795,0.0003048,0.001,3.175e-05};
Point(43) = {7.62e-05,0.0003048,0.001,0.0003048};
Point(44) = {0,0.0003048,0.001,0.0003048};
Point(45) = {-7.62e-05,0.0003048,0.001,0.0003048};
Point(46) = {-0.00010795,0.0003048,0.001,3.175e-05};
Point(47) = {-0.0001397,0.0003048,0.001,3.175e-05};
Point(48) = {-0.00017145,0.0003048,0.001,3.175e-05};
Point(49) = {-0.0002032,0.0003048,0.001,0.00015145};
Point(50) = {-0.0049997,0.0003048,0.001,0.00015145};
Point(51) = {-0.000121369128953229,0.00033655,0.001,3.175e-05};
Point(52) = {-8.96191289532294e-05,0.00033655,0.001,3.175e-05};
Point(53) = {-5.78691289532294e-05,0.00033655,0.001,0.0003048};
Point(54) = {0,0.00033655,0.001,0.0003048};
Point(55) = {5.78691289532294e-05,0.00033655,0.001,0.0003048};
Point(56) = {8.96191289532294e-05,0.00033655,0.001,3.175e-05};
Point(57) = {0.000121369128953229,0.00033655,0.001,3.175e-05};
Point(58) = {0.0049997,0.0006058,0.001,0.00015145};
Point(59) = {0,0.0006058,0.001,0.00015145};
Point(60) = {-0.0049997,0.0006058,0.001,0.00015145};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {1,8};
Line(8) = {7,20};
Line(9) = {8,9};
Line(10) = {9,10};
Line(11) = {10,11};
Line(12) = {11,12};
Line(13) = {12,13};
Line(14) = {13,14};
Line(15) = {14,15};
Line(16) = {15,16};
Line(17) = {16,17};
Line(18) = {17,18};
Line(19) = {18,19};
Line(20) = {19,20};
Line(21) = {11,21};
Line(22) = {17,27};
Line(23) = {21,22};
Line(24) = {22,23};
Line(25) = {23,24};
Line(26) = {24,25};
Line(27) = {25,26};
Line(28) = {26,27};
Line(29) = {8,28};
Line(30) = {20,30};
Line(31) = {28,29};
Line(32) = {29,30};
Line(33) = {1,31};
Line(34) = {31,32};
Line(35) = {32,2};
Line(36) = {32,33};
Line(37) = {33,3};
Line(38) = {33,34};
Line(39) = {34,4};
Line(40) = {34,35};
Line(41) = {35,5};
Line(42) = {35,36};
Line(43) = {36,6};
Line(44) = {36,37};
Line(45) = {37,7};
Line(46) = {37,38};
Line(47) = {38,20};
Line(48) = {19,39};
Line(49) = {39,38};
Line(50) = {18,40};
Line(51) = {40,39};
Line(52) = {17,41};
Line(53) = {41,40};
Line(54) = {16,42};
Line(55) = {42,41};
Line(56) = {15,43};
Line(57) = {43,42};
Line(58) = {14,44};
Line(59) = {44,43};
Line(60) = {13,45};
Line(61) = {45,44};
Line(62) = {12,46};
Line(63) = {46,45};
Line(64) = {11,47};
Line(65) = {47,46};
Line(66) = {10,48};
Line(67) = {48,47};
Line(68) = {9,49};
Line(69) = {49,48};
Line(70) = {8,50};
Line(71) = {50,49};
Line(72) = {31,50};
Line(73) = {47,51};
Line(74) = {51,21};
Line(75) = {51,52};
Line(76) = {52,22};
Line(77) = {52,53};
Line(78) = {53,23};
Line(79) = {53,54};
Line(80) = {54,24};
Line(81) = {54,55};
Line(82) = {55,25};
Line(83) = {55,56};
Line(84) = {56,26};
Line(85) = {56,57};
Line(86) = {57,27};
Line(87) = {41,57};
Line(88) = {38,58};
Line(89) = {58,30};
Line(90) = {29,59};
Line(91) = {59,58};
Line(92) = {28,60};
Line(93) = {60,59};
Line(94) = {50,60};
Curve Loop(1) = {1,2,3,4,5,6,8,-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-7};
Plane Surface(1) = {1};
Curve Loop(2) = {9,10,11,21,23,24,25,26,27,28,-22,18,19,20,30,-32,-31,-29};
Plane Surface(2) = {2};
Curve Loop(3) = {34,36,38,40,42,44,46,-49,-51,-53,-55,-57,-59,-61,-63,-65,-67,-69,-71,-72};
Plane Surface(3) = {3};
Curve Loop(4) = {1,-33,-34,-35};
Plane Surface(4) = {4};
Curve Loop(5) = {2,35,-36,-37};
Plane Surface(5) = {5};
Curve Loop(6) = {3,37,-38,-39};
Plane Surface(6) = {6};
Curve Loop(7) = {4,39,-40,-41};
Plane Surface(7) = {7};
Curve Loop(8) = {5,41,-42,-43};
Plane Surface(8) = {8};
Curve Loop(9) = {6,43,-44,-45};
Plane Surface(9) = {9};
Curve Loop(10) = {8,45,-46,-47};
Plane Surface(10) = {10};
Curve Loop(11) = {20,-48,-49,-47};
Plane Surface(11) = {11};
Curve Loop(12) = {19,-50,-51,48};
Plane Surface(12) = {12};
Curve Loop(13) = {18,-52,-53,50};
Plane Surface(13) = {13};
Curve Loop(14) = {17,-54,-55,52};
Plane Surface(14) = {14};
Curve Loop(15) = {16,-56,-57,54};
Plane Surface(15) = {15};
Curve Loop(16) = {15,-58,-59,56};
Plane Surface(16) = {16};
Curve Loop(17) = {14,-60,-61,58};
Plane Surface(17) = {17};
Curve Loop(18) = {13,-62,-63,60};
Plane Surface(18) = {18};
Curve Loop(19) = {12,-64,-65,62};
Plane Surface(19) = {19};
Curve Loop(20) = {11,-66,-67,64};
Plane Surface(20) = {20};
Curve Loop(21) = {10,-68,-69,66};
Plane Surface(21) = {21};
Curve Loop(22) = {9,-70,-71,68};
Plane Surface(22) = {22};
Curve Loop(23) = {7,-33,-72,70};
Plane Surface(23) = {23};
Curve Loop(24) = {71,69,67,73,75,77,79,81,83,85,-87,53,51,49,88,-91,-93,-94};
Plane Surface(24) = {24};
Curve Loop(25) = {21,-64,-73,-74};
Plane Surface(25) = {25};
Curve Loop(26) = {23,74,-75,-76};
Plane Surface(26) = {26};
Curve Loop(27) = {24,76,-77,-78};
Plane Surface(27) = {27};
Curve Loop(28) = {25,78,-79,-80};
Plane Surface(28) = {28};
Curve Loop(29) = {26,80,-81,-82};
Plane Surface(29) = {29};
Curve Loop(30) = {27,82,-83,-84};
Plane Surface(30) = {30};
Curve Loop(31) = {28,84,-85,-86};
Plane Surface(31) = {31};
Curve Loop(32) = {22,-52,-87,-86};
Plane Surface(32) = {32};
Curve Loop(33) = {30,47,-88,-89};
Plane Surface(33) = {33};
Curve Loop(34) = {32,-90,-91,-89};
Plane Surface(34) = {34};
Curve Loop(35) = {31,-92,-93,90};
Plane Surface(35) = {35};
Curve Loop(36) = {29,-70,-94,92};
Plane Surface(36) = {36};
Surface Loop(1) = {1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
Volume(1) = {1};
Surface Loop(2) = {2,24,22,21,20,25,26,27,28,29,30,31,32,13,12,11,33,34,35,36};
Volume(2) = {2};
Physical Volume("core") = {1};
Physical Volume("prepreg") = {2};