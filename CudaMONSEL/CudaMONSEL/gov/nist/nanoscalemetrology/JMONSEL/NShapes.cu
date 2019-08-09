#include "gov\nist\nanoscalemetrology\JMONSEL\NShapes.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalMultiPlaneShape.cuh"
#include "gov\nist\microanalysis\NISTMonte\MultiPlaneShape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalIntersectionShape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalCylindricalShape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalUnionShape.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"

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

   __host__ __device__ static void getLineSegment(const PlaneT& bottom, const PlaneT& length, const PlaneT& rightEnd, const PlaneT& leftEnd, MultiPlaneShape::LineShape& segment)
   {
      if (MultiPlaneShape::intersect3D_2Planes(bottom, length, segment) != 2) printf("createLineProjection: error 1\n");

      double extensionVec[3];
      Math2::minus3d(segment.P1, segment.P0, extensionVec);
      Math2::minus3d(segment.P0, extensionVec, segment.P0);
      Math2::plus3d(segment.P1, extensionVec, segment.P1);

      double rightIntersect[3];
      if (MultiPlaneShape::intersect3D_SegmentPlane(segment, rightEnd, rightIntersect) != 1) printf("createLineProjection: error 2\n");
      double leftIntersect[3];
      if (MultiPlaneShape::intersect3D_SegmentPlane(segment, leftEnd, leftIntersect) != 1) printf("createLineProjection: error 3\n");

      memcpy(segment.P0, rightIntersect, sizeof(segment.P0[0]) * 3);
      memcpy(segment.P1, leftIntersect, sizeof(segment.P1[0]) * 3);
   }

   __host__ __device__ Line::Line(double topz, double width, double length, double thetal, double thetar, double radl, double radr) :
      topz(topz), width(width), length(length), thetal(thetal), thetar(thetar), radl(radl), radr(radr),
      enclosure(nullptr), pl0(nullptr), pl1(nullptr), pl2(nullptr), pl3(nullptr),
      rightNMPS(nullptr), rightSide(nullptr), pl4(nullptr), plr(nullptr), rcylinder(nullptr),
      leftNMPS(nullptr), leftSide(nullptr), pl5(nullptr), pll(nullptr), lcylinder(nullptr),
      nts(nullptr), nis(nullptr),
      segment1(nullptr), segment2(nullptr), segment3(nullptr), segment4(nullptr)
   {
      create();
      //project();
   }

   __host__ __device__ Line::~Line()
   {
      delete this->enclosure; this->enclosure = nullptr;
      delete this->pl0; this->pl0 = nullptr;
      delete this->pl1; this->pl1 = nullptr;
      delete this->pl2; this->pl2 = nullptr;
      delete this->pl3; this->pl3 = nullptr;

      delete this->rightNMPS; this->rightNMPS = nullptr;
      if (this->rightSide) {
         delete (NormalUnionShapeT*)this->rightSide;
         this->rightSide = nullptr;
      }
      if (this->pl4) {
         delete this->pl4;
         this->pl4 = nullptr;
      }
      if (this->plr) {
         delete this->plr;
         this->plr = nullptr;
      }
      if (this->rcylinder) {
         delete this->rcylinder;
         this->rcylinder = nullptr;
      }

      delete this->leftNMPS; this->leftNMPS = nullptr;
      if (this->leftSide) {
         delete (NormalUnionShapeT*)this->leftSide;
         this->leftSide = nullptr;
      }
      if (this->pl5) {
         delete this->pl5;
         this->pl5 = nullptr;
      }
      if (this->pll) {
         delete this->pll;
         this->pll = nullptr;
      }
      if (this->lcylinder) {
         delete this->lcylinder;
         this->lcylinder = nullptr;
      }

      delete this->nts; this->nts = nullptr;
      delete this->nis; this->nis = nullptr;

      if (segment1) delete segment1;
      if (segment2) delete segment2;
      if (segment3) delete segment3;
      if (segment4) delete segment4;
   }

   __host__ __device__ void Line::create()
   {
      // Parameter checks
      if (radr < 0.)
         radr = 0.;
      if (radl < 0.)
         radl = 0.;

      /* First, construct the enclosure */
      enclosure = new NormalMultiPlaneShapeT();
      double signz = !topz ? 0 : (topz > 0 ? 1 : -1);
      if (signz == 0.) signz = 1.; // For rare case of 0-height specification
      const double n0[] = { 0., 0., signz }, p0[] = { 0., 0., topz }; // Add top plane
      pl0 = new PlaneT(n0, p0);
      enclosure->addPlane(pl0);
      const double n1[] = { 0., 0., -signz }, p1[] = { 0., 0., 0. }; // Add bottom plane
      pl1 = new PlaneT(n1, p1);
      enclosure->addPlane(pl1);
      // Add end caps
      const double n2[] = { 0., 1., 0. }, p2[] = { 0., length / 2., 0. }; // Right end
      pl2 = new PlaneT(n2, p2);
      enclosure->addPlane(pl2);
      const double n3[] = { 0., -1., 0. }, p3[] = { 0., -length / 2., 0. }; // Left end
      pl3 = new PlaneT(n3, p3);
      enclosure->addPlane(pl3);

      rightNMPS = new NormalMultiPlaneShapeT();
      rightSide = nullptr;

      // Add right sidewall
      const double costhetar = ::cos(thetar);
      const double sinthetar = ::sin(thetar);
      const double n4[] = { costhetar, 0., signz * sinthetar }, p4[] = { width / 2, 0., 0. };
      pl4 = new PlaneT(n4, p4);
      rightNMPS->addPlane(pl4);
      // If radr>0 add a clipping plane and the cylinder
      const double root2 = ::sqrt(2.);
      const double absz = signz * topz;
      if (radr > 0) {
         const double rad = ::sqrt(1 - sinthetar);
         const double nr[] = { rad / root2, 0., (signz * costhetar) / root2 / rad }, pr[] = { ((width / 2.) - (radr / costhetar)) + (((radr - absz) * sinthetar) / costhetar), 0., topz };
         plr = new PlaneT(nr, pr);
         rightNMPS->addPlane(plr);
         // Construct cylinder for right corner
         const double xc = ((width / 2.) - (radr / ::cos(thetar))) + ((radr - absz) * ::tan(thetar));
         const double zc = topz - (signz * radr);
         const double end0r[] = { xc, -length / 2., zc }, end1r[] = { xc, length / 2., zc };
         rcylinder = new NormalCylindricalShapeT(end0r, end1r, radr);
         rightSide = new NormalUnionShapeT(*rightNMPS, *rcylinder);
      }
      else
         rightSide = rightNMPS;

      leftNMPS = new NormalMultiPlaneShapeT();
      leftSide = nullptr;

      // Add left sidewall
      const double costhetal = ::cos(thetal);
      const double sinthetal = ::sin(thetal);
      const double n5[] = { -costhetal, 0., signz * sinthetal }, p5[] = { -width / 2, 0., 0. };
      pl5 = new PlaneT(n5, p5);
      leftNMPS->addPlane(pl5);
      // If radl>0 add a clipping plane and the cylinder
      if (radl > 0.) {
         const double rad = ::sqrt(1 - sinthetal);
         const double nl[] = { -rad / root2, 0., (signz * costhetal) / root2 / rad }, pl[] = { ((-width / 2.) + (radl / costhetal)) - (((radl - absz) * sinthetal) / costhetal), 0., topz };
         pll = new PlaneT(nl, pl);
         leftNMPS->addPlane(pll);
         const double xc = ((width / 2.) - (radl / ::cos(thetal))) + ((radl - absz) * ::tan(thetal));
         const double zc = topz - (signz * radl);
         // Construct cylinder for left corner
         const double end0[] = { -xc, -length / 2., zc }, end1[] = { -xc, length / 2., zc };
         lcylinder = new NormalCylindricalShapeT(end0, end1, radl);
         leftSide = new NormalUnionShapeT(*leftNMPS, *lcylinder);
      }
      else
         leftSide = leftNMPS;

      nts = new NormalIntersectionShapeT(*leftSide, *rightSide);
      nis = new NormalIntersectionShapeT(*nts, *enclosure);
   }

   /*
   * Can only rotate along z-axis
   */
   __host__ __device__ void Line::createProjection()
   {
      segment1 = new LineShapeT();
      segment2 = new LineShapeT();
      segment3 = new LineShapeT();
      segment4 = new LineShapeT();

      getLineSegment(*pl1, *pl2, *pl4, *pl5, *segment1); // right cap
      getLineSegment(*pl1, *pl3, *pl4, *pl5, *segment2); // left cap
      getLineSegment(*pl1, *pl4, *pl2, *pl3, *segment3); // right length
      getLineSegment(*pl1, *pl5, *pl2, *pl3, *segment4); // left length

      printf("(%.5e, %.5e, %.5e) -> (%.5e, %.5e, %.5e)\n", segment1->P0[0], segment1->P0[1], segment1->P0[2], segment1->P1[0], segment1->P1[1], segment1->P1[2]);
      printf("(%.5e, %.5e, %.5e) -> (%.5e, %.5e, %.5e)\n", segment2->P0[0], segment2->P0[1], segment2->P0[2], segment2->P1[0], segment2->P1[1], segment2->P1[2]);
      printf("(%.5e, %.5e, %.5e) -> (%.5e, %.5e, %.5e)\n", segment3->P0[0], segment3->P0[1], segment3->P0[2], segment3->P1[0], segment3->P1[1], segment3->P1[2]);
      printf("(%.5e, %.5e, %.5e) -> (%.5e, %.5e, %.5e)\n", segment4->P0[0], segment4->P0[1], segment4->P0[2], segment4->P1[0], segment4->P1[1], segment4->P1[2]);
   }

   __host__ __device__ NormalIntersectionShapeT* Line::get()
   {
      return nis;
   }
}