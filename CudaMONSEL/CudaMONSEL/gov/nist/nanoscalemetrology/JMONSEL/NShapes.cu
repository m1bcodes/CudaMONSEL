//#include "gov\nist\nanoscalemetrology\JMONSEL\NShapes.cuh"
//#include "gov\nist\nanoscalemetrology\JMONSEL\NormalMultiPlaneShape.cuh"
//#include "gov\nist\microanalysis\NISTMonte\MultiPlaneShape.cuh"
//#include "gov\nist\nanoscalemetrology\JMONSEL\NormalIntersectionShape.cuh"
//
//static VectorXd normalize(double vec[])
//{
//   double res[3];
//   double norm = ::sqrt((vec[0] * vec[0]) + (vec[1] * vec[1]) + (vec[2] * vec[2]));
//   res[0] = vec[0] / norm;
//   res[1] = vec[1] / norm;
//   res[2] = vec[2] / norm;
//   return VectorXd(res, res+3);
//}
//
//// returns a vector pointing in the opposite direction of v
//static VectorXd invert(double v[])
//{
//   double res[3];
//   res[0] = -v[0];
//   res[1] = -v[1];
//   res[2] = -v[2];
//   return VectorXd(res, res + 3);
//}
//
//// Add a plane offset by dist*normal from pt
////static void addOffsetPlane(NormalMultiPlaneShapeT& shape, double normal[], double pt[], double dist)
////{
////   normal = normalize(normal).data();
////   shape.addPlane(normal, new double[] {
////      pt[0] + (normal[0] * dist),
////         pt[1] + (normal[1] * dist),
////         pt[2] + (normal[2] * dist)
////   });
////}
//
//static void addPlane(NormalMultiPlaneShapeT& shape, PlaneT& plane)
//{
//   shape.addPlane(plane);
//}
//
////NormalMultiPlaneShapeT createNormalFilm(double normal[], double pt1[], double thickness)
////{
////   NormalMultiPlaneShapeT mp;
////   mp.addPlane(normal, pt1);
////   addOffsetPlane(mp, invert(normal), pt1, thickness);
////   return mp;
////}
//
//NormalShapeT* createLine(double topz, // z of the top face
//   double width, // line width
//   double length, // length of line
//   double thetal, // angle of right sidewall
//   double thetar, // angle of left sidewall
//   double radl, // radius of top right corner
//   double radr // radius of top left corner
//   ) {
//   // Parameter checks
//   if (radr < 0.)
//      radr = 0.;
//   if (radl < 0.)
//      radl = 0.;
//
//   /*
//   * The line will be the intersection of 3 pieces, a right side, a left
//   * side, and an "enclosure". The right side is a multiplane shape with 2
//   * planes, one representing the sidewall and the other a plane that joins
//   * the top of the line to the place where the cylinder touches the
//   * sidewall. The left side does the same for the other side. The enclosure
//   * is a NormalMultiPlaneShape shape consisting of the top, bottom, and
//   * endcaps. I replaced the previous createLine algorithm with this one on
//   * 2/14/2013. The previous one formed the union of two cylinders,
//   * representing the corner rounding, with one multiplane shape (8 planes)
//   * to form the line. This worked fine for lines that were wide enough, but
//   * when the corner rounding becomes large enough (e.g., in the extreme
//   * when the cylinder diameter is greater than the linewidth) the cylinder
//   * that rounds the right corner can protrude through the left sidewall.
//   * This new algorithm eliminates that issue at the cost of being slightly
//   * (less than 5%) slower. This is a price worth paying now that narrow
//   * lines are much more frequently done.
//   */
//
//   /* First, construct the enclosure */
//   NormalMultiPlaneShapeT enclosure;
//   // Add top plane
//   double signz = !topz ? 0 : (topz > 0 ? 1 : -1);
//   if (signz == 0.)
//      signz = 1.; // For rare case of 0-height specification
//   double n0[] = { 0., 0., signz }, p0[] = { 0., 0., topz };
//   PlaneT pl0(n0, 3, p0, 3);
//   enclosure.addPlane(pl0);
//   // Add bottom plane
//   double n1[] = { 0., 0., -signz }, p1[] = { 0., 0., 0. };
//   PlaneT pl1(n1, 3, p1, 3);
//   enclosure.addPlane(pl1);
//   // Add end caps
//   double n2[] = { 0., 1., 0. }, p2[] = { 0., length / 2., 0. }; // Right end 
//   PlaneT pl2(n2, 3, p2, 3);
//   enclosure.addPlane(pl2);
//   double n3[] = { 0., -1., 0. }, p3[] = { 0., -length / 2., 0. }; // Left end
//   PlaneT pl3(n3, 3, p3, 3);
//   enclosure.addPlane(pl3);
//
//   /* Now do the right side */
//
//   NormalMultiPlaneShapeT rightNMPS;
//   NormalShapeT* rightSide = NULL;
//
//   // Add right sidewall
//   double costhetar = ::cos(thetar);
//   double sinthetar = ::sin(thetar);
//
//   double n4[] = { costhetar, 0., signz * sinthetar }, p4[] = { width / 2, 0., 0. };
//   PlaneT pl4(n4, 3, p4, 3);
//   rightNMPS.addPlane(pl4);
//   // If radr>0 add a clipping plane and the cylinder
//   double root2 = ::sqrt(2.);
//   double absz = signz * topz;
//   if (radr > 0) {
//      double rad = ::sqrt(1 - sinthetar);
//      double n5[] = { rad / root2, 0., (signz * costhetar) / root2 / rad }, p5[] = { ((width / 2.) - (radr / costhetar)) + (((radr - absz) * sinthetar) / costhetar), 0., topz };
//      PlaneT pl5(n5, 3, p5, 3);
//      rightNMPS.addPlane(pl5);
//      // Construct cylinder for right corner
//      double xc = ((width / 2.) - (radr / ::cos(thetar))) + ((radr - absz) * ::tan(thetar));
//      double zc = topz - (signz * radr);
//      double n6[] = { xc, -length / 2., zc }, p6[] = { xc, length / 2., zc };
//      NormalCylindricalShapeT rcylinder(n6, 3, p6, 3, radr);
//      rightSide = new NormalUnionShape(rightNMPS, rcylinder);
//   }
//   else
//      rightSide = &rightNMPS;
//
//   /* Now do likewise for the left side */
//   NormalMultiPlaneShapeT leftNMPS;
//   NormalShapeT* leftSide;
//
//   // Add left sidewall
//   double costhetal = ::cos(thetal);
//   double sinthetal = ::sin(thetal);
//   double n7[] = { -costhetal, 0., signz * sinthetal }, p7[] = { -width / 2, 0., 0. };
//   PlaneT pl7(n7, 3, p7, 3);
//   leftNMPS.addPlane(pl7);
//   // If radl>0 add a clipping plane and the cylinder
//   if (radl > 0.) {
//      double rad = ::sqrt(1 - sinthetal);
//      double n8[] = { -rad / root2, 0., (signz * costhetal) / root2 / rad }, p8[] = { ((-width / 2.) + (radl / costhetal)) - (((radl - absz) * sinthetal) / costhetal), 0., topz };
//      PlaneT pl8(n8, 3, p8, 3);
//      leftNMPS.addPlane(pl8);
//      double xc = ((width / 2.) - (radl / ::cos(thetal))) + ((radl - absz) * ::tan(thetal));
//      double zc = topz - (signz * radl);
//      // Construct cylinder for left corner
//      double n9[] = { -xc, -length / 2., zc }, p9[] = { -xc, length / 2., zc };
//      PlaneT pl9(n9, 3, p9, 3);
//      NormalCylindricalShapeT lcylinder(pl9, radl);
//      leftSide = new NormalUnionShape(leftNMPS, lcylinder);
//   }
//   else
//      leftSide = leftNMPS;
//
//   /* The shape is the intersection of the 3 shapes just constructed */
//   auto nts = NormalIntersectionShapeT(*leftSide, *rightSide);
//   auto nis = NormalIntersectionShapeT(nts, enclosure);
//   return &nis;
//}
