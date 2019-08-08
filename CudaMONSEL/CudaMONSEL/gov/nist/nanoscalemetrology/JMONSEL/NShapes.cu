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

   //static void addPlane(NormalMultiPlaneShapeT& shape, PlaneT *plane)
   //{
   //   shape.addPlane(plane);
   //}

   //NormalMultiPlaneShapeT createNormalFilm(double normal[], double pt1[], double thickness)
   //{
   //   NormalMultiPlaneShapeT mp;
   //   mp.addPlane(normal, pt1);
   //   addOffsetPlane(mp, invert(normal), pt1, thickness);
   //   return mp;
   //}

   __host__ __device__ NormalShapeT* createLine(
      Line& line,
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
      line.enclosure = new NormalMultiPlaneShapeT();
      // Add top plane
      double signz = !topz ? 0 : (topz > 0 ? 1 : -1);
      if (signz == 0.) signz = 1.; // For rare case of 0-height specification
      const double n0[] = { 0., 0., signz }, p0[] = { 0., 0., topz };
      line.pl0 = new PlaneT(n0, p0);
      line.enclosure->addPlane(line.pl0);
      // Add bottom plane
      const double n1[] = { 0., 0., -signz }, p1[] = { 0., 0., 0. };
      line.pl1 = new PlaneT(n1, p1);
      line.enclosure->addPlane(line.pl1);
      // Add end caps
      const double n2[] = { 0., 1., 0. }, p2[] = { 0., length / 2., 0. }; // Right end
      line.pl2 = new PlaneT(n2, p2);
      line.enclosure->addPlane(line.pl2);
      const double n3[] = { 0., -1., 0. }, p3[] = { 0., -length / 2., 0. }; // Left end
      line.pl3 = new PlaneT(n3, p3);
      line.enclosure->addPlane(line.pl3);

      line.rightNMPS = new NormalMultiPlaneShapeT();
      line.rightSide = nullptr;

      // Add right sidewall
      const double costhetar = ::cos(thetar);
      const double sinthetar = ::sin(thetar);

      const double n4[] = { costhetar, 0., signz * sinthetar }, p4[] = { width / 2, 0., 0. };
      line.pl4 = new PlaneT(n4, p4);
      line.rightNMPS->addPlane(line.pl4);
      // If radr>0 add a clipping plane and the cylinder
      const double root2 = ::sqrt(2.);
      const double absz = signz * topz;
      if (radr > 0) {
         const double rad = ::sqrt(1 - sinthetar);
         const double nr[] = { rad / root2, 0., (signz * costhetar) / root2 / rad }, pr[] = { ((width / 2.) - (radr / costhetar)) + (((radr - absz) * sinthetar) / costhetar), 0., topz };
         line.plr = new PlaneT(nr, pr);
         line.rightNMPS->addPlane(line.plr);
         // Construct cylinder for right corner
         const double xc = ((width / 2.) - (radr / ::cos(thetar))) + ((radr - absz) * ::tan(thetar));
         const double zc = topz - (signz * radr);
         const double end0r[] = { xc, -length / 2., zc }, end1r[] = { xc, length / 2., zc };
         line.rcylinder = new NormalCylindricalShapeT(end0r, end1r, radr);
         line.rightSide = new NormalUnionShapeT(*line.rightNMPS, *line.rcylinder);
      }
      else
         line.rightSide = line.rightNMPS;

      line.leftNMPS = new NormalMultiPlaneShapeT();
      line.leftSide = nullptr;

      // Add left sidewall
      const double costhetal = ::cos(thetal);
      const double sinthetal = ::sin(thetal);
      const double n5[] = { -costhetal, 0., signz * sinthetal }, p5[] = { -width / 2, 0., 0. };
      line.pl5 = new PlaneT(n5, p5);
      line.leftNMPS->addPlane(line.pl5);
      // If radl>0 add a clipping plane and the cylinder
      if (radl > 0.) {
         const double rad = ::sqrt(1 - sinthetal);
         const double nl[] = { -rad / root2, 0., (signz * costhetal) / root2 / rad }, pl[] = { ((-width / 2.) + (radl / costhetal)) - (((radl - absz) * sinthetal) / costhetal), 0., topz };
         line.pll = new PlaneT(nl, pl);
         line.leftNMPS->addPlane(line.pll);
         const double xc = ((width / 2.) - (radl / ::cos(thetal))) + ((radl - absz) * ::tan(thetal));
         const double zc = topz - (signz * radl);
         // Construct cylinder for left corner
         const double end0[] = { -xc, -length / 2., zc }, end1[] = { -xc, length / 2., zc };
         line.lcylinder = new NormalCylindricalShapeT(end0, end1, radl);
         line.leftSide = new NormalUnionShapeT(*line.leftNMPS, *line.lcylinder);
      }
      else
         line.leftSide = line.leftNMPS;

      line.nts = new NormalIntersectionShapeT(*line.leftSide, *line.rightSide);
      line.nis = new NormalIntersectionShapeT(*line.nts, *line.enclosure);
      return line.nis;
   }

   __host__ __device__ extern void destroyLine(Line& line)
   {
      delete line.enclosure; line.enclosure = nullptr;
      delete line.pl0; line.pl0 = nullptr;
      delete line.pl1; line.pl1 = nullptr;
      delete line.pl2; line.pl2 = nullptr;
      delete line.pl3; line.pl3 = nullptr;

      delete line.rightNMPS; line.rightNMPS = nullptr;
      if (line.rightSide) {
         delete (NormalUnionShapeT*)line.rightSide;
         line.rightSide = nullptr;
      }
      delete line.pl4; line.pl4 = nullptr;
      delete line.plr; line.plr = nullptr;
      delete line.rcylinder; line.rcylinder = nullptr;

      delete line.leftNMPS; line.leftNMPS = nullptr;
      if (line.leftSide) {
         delete (NormalUnionShapeT*)line.leftSide;
         line.leftSide = nullptr;
      }
      delete line.pl5; line.pl5 = nullptr;
      delete line.pll; line.pll = nullptr;
      delete line.lcylinder; line.lcylinder = nullptr;

      delete line.nts; line.nts = nullptr;
      delete line.nis; line.nis = nullptr;
   }
}