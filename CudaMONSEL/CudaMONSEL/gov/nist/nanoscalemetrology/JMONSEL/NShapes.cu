#include "gov\nist\nanoscalemetrology\JMONSEL\NShapes.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalMultiPlaneShape.cuh"
#include "gov\nist\microanalysis\NISTMonte\MultiPlaneShape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalIntersectionShape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalCylindricalShape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalUnionShape.cuh"

namespace NShapes
{
   static void normalize(double vec[], double res[])
   {
      double norm = ::sqrt((vec[0] * vec[0]) + (vec[1] * vec[1]) + (vec[2] * vec[2]));
      res[0] = vec[0] / norm;
      res[1] = vec[1] / norm;
      res[2] = vec[2] / norm;
   }

   // returns a vector pointing in the opposite direction of v
   static void invert(double v[], double res[])
   {
      res[0] = -v[0];
      res[1] = -v[1];
      res[2] = -v[2];
   }

   // Add a plane offset by dist*normal from pt
   //static void addOffsetPlane(NormalMultiPlaneShapeT& shape, double normal[], double pt[], double dist)
   //{
   //   normal = normalize(normal).data();
   //   shape.addPlane(normal, new double[] {
   //      pt[0] + (normal[0] * dist),
   //         pt[1] + (normal[1] * dist),
   //         pt[2] + (normal[2] * dist)
   //   });
   //}

   static void addPlane(NormalMultiPlaneShapeT& shape, PlaneT& plane)
   {
      shape.addPlane(plane);
   }

   //NormalMultiPlaneShapeT createNormalFilm(double normal[], double pt1[], double thickness)
   //{
   //   NormalMultiPlaneShapeT mp;
   //   mp.addPlane(normal, pt1);
   //   addOffsetPlane(mp, invert(normal), pt1, thickness);
   //   return mp;
   //}

   // TODO: write a destructor
   NormalShapeT* createLine(
      double topz, // z of the top face
      double width, // line width
      double length, // length of line
      double thetal, // angle of right sidewall
      double thetar, // angle of left sidewall
      double radl, // radius of top right corner
      double radr // radius of top left corner
      ) {
      // Parameter checks
      if (radr < 0.)
         radr = 0.;
      if (radl < 0.)
         radl = 0.;

      /* First, construct the enclosure */
      NormalMultiPlaneShapeT* enclosure = new NormalMultiPlaneShapeT();
      // Add top plane
      double signz = !topz ? 0 : (topz > 0 ? 1 : -1);
      if (signz == 0.) signz = 1.; // For rare case of 0-height specification
      const double n0[] = { 0., 0., signz }, p0[] = { 0., 0., topz };
      PlaneT* pl0 = new PlaneT(n0, p0);
      enclosure->addPlane(*pl0);
      // Add bottom plane
      double n1[] = { 0., 0., -signz }, p1[] = { 0., 0., 0. };
      PlaneT* pl1 = new PlaneT(n1, p1);
      enclosure->addPlane(*pl1);
      // Add end caps
      const double n2[] = { 0., 1., 0. }, p2[] = { 0., length / 2., 0. }; // Right end 
      PlaneT* pl2 = new PlaneT(n2, p2);
      enclosure->addPlane(*pl2);
      const double n3[] = { 0., -1., 0. }, p3[] = { 0., -length / 2., 0. }; // Left end
      PlaneT* pl3 = new PlaneT(n3, p3);
      enclosure->addPlane(*pl3);

      NormalMultiPlaneShapeT* rightNMPS = new NormalMultiPlaneShapeT();
      NormalShapeT* rightSide = NULL;

      // Add right sidewall
      const double costhetar = ::cos(thetar);
      const double sinthetar = ::sin(thetar);

      const double n4[] = { costhetar, 0., signz * sinthetar }, p4[] = { width / 2, 0., 0. };
      PlaneT* pl4 = new PlaneT(n4, p4);
      rightNMPS->addPlane(*pl4);
      // If radr>0 add a clipping plane and the cylinder
      double root2 = ::sqrt(2.);
      double absz = signz * topz;
      if (radr > 0) {
         const double rad = ::sqrt(1 - sinthetar);
         const double nr[] = { rad / root2, 0., (signz * costhetar) / root2 / rad }, pr[] = { ((width / 2.) - (radr / costhetar)) + (((radr - absz) * sinthetar) / costhetar), 0., topz };
         PlaneT* plr = new PlaneT(nr, pr);
         rightNMPS->addPlane(*plr);
         // Construct cylinder for right corner
         const double xc = ((width / 2.) - (radr / ::cos(thetar))) + ((radr - absz) * ::tan(thetar));
         const double zc = topz - (signz * radr);
         const double end0r[] = { xc, -length / 2., zc }, end1r[] = { xc, length / 2., zc };
         NormalCylindricalShapeT* rcylinder = new NormalCylindricalShapeT(end0r, end1r, radr);
         rightSide = new NormalUnionShapeT(*rightNMPS, *rcylinder);
      }
      else
         rightSide = rightNMPS;

      NormalMultiPlaneShapeT* leftNMPS = new NormalMultiPlaneShapeT();
      NormalShapeT* leftSide = NULL;

      // Add left sidewall
      const double costhetal = ::cos(thetal);
      const double sinthetal = ::sin(thetal);
      const double n6[] = { -costhetal, 0., signz * sinthetal }, p6[] = { -width / 2, 0., 0. };
      PlaneT* pl6 = new PlaneT(n6, p6);
      leftNMPS->addPlane(*pl6);
      // If radl>0 add a clipping plane and the cylinder
      if (radl > 0.) {
         const double rad = ::sqrt(1 - sinthetal);
         const double n8[] = { -rad / root2, 0., (signz * costhetal) / root2 / rad }, p8[] = { ((-width / 2.) + (radl / costhetal)) - (((radl - absz) * sinthetal) / costhetal), 0., topz };
         PlaneT* pl = new PlaneT(n8, p8);
         leftNMPS->addPlane(*pl);
         const double xc = ((width / 2.) - (radl / ::cos(thetal))) + ((radl - absz) * ::tan(thetal));
         const double zc = topz - (signz * radl);
         // Construct cylinder for left corner
         const double end0[] = { -xc, -length / 2., zc }, end1[] = { -xc, length / 2., zc };
         NormalCylindricalShapeT* lcylinder = new NormalCylindricalShapeT(end0, end1, radl);
         leftSide = new NormalUnionShapeT(*leftNMPS, *lcylinder);
      }
      else
         leftSide = leftNMPS;

      NormalIntersectionShapeT* nts = new NormalIntersectionShapeT(*leftSide, *rightSide);
      NormalIntersectionShapeT* nis = new NormalIntersectionShapeT(*nts, *enclosure);
      return nis;
   }
}