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
      const double norm = ::sqrt((vec[0] * vec[0]) + (vec[1] * vec[1]) + (vec[2] * vec[2]));
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

   __host__ __device__ Line::Line(const double topz, const double width, const double length, const double thetal, const double thetar, const double radl, const double radr) :
      topz(topz), width(width), length(length), thetal(thetal), thetar(thetar), radl(radl < 0 ? 0 : radl), radr(radr < 0 ? 0 : radr),
      enclosure(nullptr), pl0(nullptr), pl1(nullptr), pl2(nullptr), pl3(nullptr),
      rightNMPS(nullptr), rightSide(nullptr), pl4(nullptr), plr(nullptr), rcylinder(nullptr),
      leftNMPS(nullptr), leftSide(nullptr), pl5(nullptr), pll(nullptr), lcylinder(nullptr),
      nts(nullptr), nis(nullptr),
      segment0(nullptr), segment1(nullptr), segment2(nullptr), segment3(nullptr)
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

      if (segment0) { delete segment0; segment0 = nullptr; }
      if (segment1) { delete segment1; segment1 = nullptr; }
      if (segment2) { delete segment2; segment2 = nullptr; }
      if (segment3) { delete segment3; segment3 = nullptr; }
   }

   __host__ __device__ void Line::create()
   {
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
      const double n4[] = { costhetar, 0., signz * sinthetar }, p4[] = { width / 2., 0., 0. };
      pl4 = new PlaneT(n4, p4);
      rightNMPS->addPlane(pl4);
      // If radr>0 add a clipping plane and the cylinder
      const double root2 = ::sqrt(2.);
      const double absz = signz * topz;
      if (radr > 0) {
         const double rad = ::sqrt(1. - sinthetar);
         const double nr[] = { rad / root2, 0., signz * costhetar / root2 / rad }, pr[] = { width / 2. - radr / costhetar + (radr - absz) * sinthetar / costhetar, 0., topz };
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
      const double n5[] = { -costhetal, 0., signz * sinthetal }, p5[] = { -width / 2., 0., 0. };
      pl5 = new PlaneT(n5, p5);
      leftNMPS->addPlane(pl5);
      // If radl>0 add a clipping plane and the cylinder
      if (radl > 0.) {
         const double rad = ::sqrt(1. - sinthetal);
         const double nl[] = { -rad / root2, 0., signz * costhetal / root2 / rad }, pl[] = { -width / 2. + radl / costhetal - (radl - absz) * sinthetal / costhetal, 0., topz };
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

   __host__ __device__ void Line::calcGroundtruth()
   {
      if (!segment0) segment0 = new LineShapeT();
      if (!segment1) segment1 = new LineShapeT();
      if (!segment2) segment2 = new LineShapeT();
      if (!segment3) segment3 = new LineShapeT();

      getLineSegment(*pl0, *pl2, *pl4, *pl5, *segment0); // right cap
      getLineSegment(*pl0, *pl3, *pl4, *pl5, *segment1); // left cap
      getLineSegment(*pl0, *pl4, *pl2, *pl3, *segment2); // right length
      getLineSegment(*pl0, *pl5, *pl2, *pl3, *segment3); // left length

      //getLineSegment(*pl1, *pl2, *pl4, *pl5, *segment0); // right cap
      //getLineSegment(*pl1, *pl3, *pl4, *pl5, *segment1); // left cap
      //getLineSegment(*pl1, *pl4, *pl2, *pl3, *segment2); // right length
      //getLineSegment(*pl1, *pl5, *pl2, *pl3, *segment3); // left length

      printf("(%.5e, %.5e, %.5e) -> (%.5e, %.5e, %.5e)\n", segment0->P0[0], segment0->P0[1], segment0->P0[2], segment0->P1[0], segment0->P1[1], segment0->P1[2]);
      printf("(%.5e, %.5e, %.5e) -> (%.5e, %.5e, %.5e)\n", segment1->P0[0], segment1->P0[1], segment1->P0[2], segment1->P1[0], segment1->P1[1], segment1->P1[2]);
      printf("(%.5e, %.5e, %.5e) -> (%.5e, %.5e, %.5e)\n", segment2->P0[0], segment2->P0[1], segment2->P0[2], segment2->P1[0], segment2->P1[1], segment2->P1[2]);
      printf("(%.5e, %.5e, %.5e) -> (%.5e, %.5e, %.5e)\n", segment3->P0[0], segment3->P0[1], segment3->P0[2], segment3->P1[0], segment3->P1[1], segment3->P1[2]);
   }

   __host__ __device__ static void getPerpendicularIntersection(
      const PlaneT& plane, // the projection plane
      const double* p, // absolute coordinate
      double* res // the absolute coordinate of the perpendicular intersection between p and plane
      )
   {
      // checks
      const double* n0 = plane.getNormal();
      const double* p0 = plane.getPoint();
      const float len = ::sqrtf(Math2::dot3d(n0, n0));
      if (len - 1.f > .000001f) printf("Line::calcProjection: plane normal length not unity %.5e\n", len);

      // calc
      double v[3];
      Math2::minus3d(p, p0, v);
      double s = Math2::dot3d(v, n0); // perpendicular distance between p and plane
      s = s == 0. ? 1. : s;
      MultiPlaneShape::LineShape l;
      static const float extensionFactor = 10.f;
      l.P0[0] = p[0] + extensionFactor * s * n0[0]; l.P0[1] = p[1] + extensionFactor * s * n0[1]; l.P0[2] = p[2] + extensionFactor * s * n0[2];
      l.P1[0] = p[0] - extensionFactor * s * n0[0]; l.P1[1] = p[1] - extensionFactor * s * n0[1]; l.P1[2] = p[2] - extensionFactor * s * n0[2];
      const int ind = MultiPlaneShape::intersect3D_SegmentPlane(l, plane, res);

      // checks
      if (ind == 0) printf("getPerpendicularIntersection: no intersect\n");
      else if (ind == 1) printf("getPerpendicularIntersection: intersect\n");
      else if (ind == 2) printf("getPerpendicularIntersection: on plane\n");
   }

   __host__ __device__ static void getLineProjection(
      const PlaneT& plane, // the projection plane
      const MultiPlaneShape::LineShape& line, // the absolute coordinate of the line to be projected
      MultiPlaneShape::LineShape& res // the absolute coordinate of the perpendicular intersection between line and plane
      )
   {
      getPerpendicularIntersection(plane, line.P0, res.P0);
      getPerpendicularIntersection(plane, line.P1, res.P1);
   }

   // get vectors relative to plane origin
   __host__ __device__ static void getRelativeCoord(
      const PlaneT& plane, // the axis-containing plane
      const double* axis0, // coordinate relative to plane origin (lies on plane)
      const double* axis1, // coordinate relative to plane origin orthogonal to axis0 and plane normal (lies on plane)
      const double* p, // absolute coordinate of a point on the plane
      double* res
      )
   {
      // check
      const float len = Math2::dot3d(plane.getNormal(), plane.getNormal());
      if (::fabsf(1 - len) > 0.00001f) printf("getRelativeCoord: plane normal not unity: %.5e\n", len);
      if (::fabsf(Math2::dot3d(plane.getNormal(), axis0)) > 0.000001f) printf("getRelativeCoord: axis0 not coplanar with plane\n");
      if (::fabsf(Math2::dot3d(plane.getNormal(), axis1)) > 0.000001f) printf("getRelativeCoord: axis1 not coplanar with plane\n");
      double checkNormal[3];
      Math2::cross3d(axis0, axis1, checkNormal);
      const float t = ::fabsf(Math2::dot3d(checkNormal, plane.getNormal()));
      if (::fabsf(1 - t) > 0.000001f) printf("getRelativeCoord(): t = %.5e\n", t);

      // check
      Math2::minus3d(p, plane.getPoint(), res);
      const double tmp[3] = { Math2::dot3d(res, axis0), Math2::dot3d(res, axis1), Math2::dot3d(res, plane.getNormal()) };
      memcpy(res, tmp, sizeof(res[0]) * 3);

      // check
      if (::fabsf(res[2]) > 0.00001f) printf("getRelativeCoord: res.P0 has out of plane component: %.5e\n", res[2]);
   }

   __host__ __device__ static void getRelativeProjection(
      const PlaneT& plane,
      const double* axis0, // vector from the plane origin
      const double* axis1, // vector from the plane origin
      const MultiPlaneShape::LineShape& line, // line in absolute coordinates to be projected onto plane
      MultiPlaneShape::LineShape& res // line with respect to origin of plane (ie. res + plane.origin = line)
      )
   {
      getRelativeCoord(plane, axis0, axis1, line.P0, res.P0);
      getRelativeCoord(plane, axis0, axis1, line.P1, res.P1);
   }

   __host__ __device__ static void calcLineProjection(
      const PlaneT& plane,
      const double* axis0, // vector from the plane origin
      const double* axis1, // vector from the plane origin
      const MultiPlaneShape::LineShape& line, // line in absolute coordinates to be projected onto plane
      MultiPlaneShape::LineShape& res // line with respect to origin of plane (ie. res + plane.origin = line)
      )
   {
      MultiPlaneShape::LineShape tmp;
      getLineProjection(plane, line, tmp);
      getRelativeProjection(plane, axis0, axis1, tmp, res);
   }

   struct Point
   {
      __host__ __device__ Point(int x, int y) : X(x), Y(y) {}

      int X, Y;
   };

   // Bresenham's algorithm
   __host__ __device__ void DrawLine(const Point& p0, const Point& p1, const char color, char* res, const unsigned int w, const unsigned int h)
   {
      int x = p0.X;
      int y = p0.Y;
      int dx = ::fabsf(p1.X - x);
      int dy = ::fabsf(p1.Y - y);
      const int s1 = (p1.X - x) < 0 ? -1 : ((p1.X - x) > 0 ? 1 : 0); // get the sign
      const int s2 = (p1.Y - y) < 0 ? -1 : ((p1.Y - y) > 0 ? 1 : 0); // get the sign
      bool interchange = dy > dx;
      if (interchange) {
         auto t = dy;
         dy = dx;
         dx = t;
      }
      int e = 2 * dy - dx;
      for (int i = 0; i < dx; i++) {
         if (y >= 0 && y < h && x >= 0 && x < w) res[y * w + x] = color;
         while (e >= 0) {
            if (interchange) {
               x += s1;
            }
            else {
               y += s2;
            }
            e += -2. * dx;
         }
         if (interchange) {
            y += s2;
         }
         else {
            x += s1;
         }
         e += 2. * dy;
      }
      if (p1.Y >= 0 && p1.Y < h && p1.X >= 0 && p1.X < w) res[p1.Y * w + p1.X] = color;
   }

   /*
   * Can only rotate along z-axis
   */
   // axis2 would be the normal of the plane, which creates no projection
   __host__ __device__ void Line::calcRasterization(
      const PlaneT& plane,
      const double* axis0, // vector from the plane origin
      const double* axis1, // vector from the plane origin
      const float xlenperpix,
      const float ylenperpix,
      char* res,
      const unsigned int w,
      const unsigned int h
      )
   {
      MultiPlaneShape::LineShape res0, res1, res2, res3; // with respect to origin of plane, along axis0 and axis1
      calcLineProjection(plane, axis0, axis1, *segment0, res0);
      calcLineProjection(plane, axis0, axis1, *segment1, res1);
      calcLineProjection(plane, axis0, axis1, *segment2, res2);
      calcLineProjection(plane, axis0, axis1, *segment3, res3);

      int start0[3] = { res0.P0[0] / xlenperpix, res0.P0[1] / ylenperpix, res0.P0[2] };
      int end0[3] = { res0.P1[0] / xlenperpix, res0.P1[1] / ylenperpix, res0.P1[2] };
      int start1[3] = { res1.P0[0] / xlenperpix, res1.P0[1] / ylenperpix, res1.P0[2] };
      int end1[3] = { res1.P1[0] / xlenperpix, res1.P1[1] / ylenperpix, res1.P1[2] };
      int start2[3] = { res2.P0[0] / xlenperpix, res2.P0[1] / ylenperpix, res2.P0[2] };
      int end2[3] = { res2.P1[0] / xlenperpix, res2.P1[1] / ylenperpix, res2.P1[2] };
      int start3[3] = { res3.P0[0] / xlenperpix, res3.P0[1] / ylenperpix, res3.P0[2] };
      int end3[3] = { res3.P1[0] / xlenperpix, res3.P1[1] / ylenperpix, res3.P1[2] };

      Point s0(start0[0], start0[1]);
      Point s1(start1[0], start1[1]);
      Point s2(start2[0], start2[1]);
      Point s3(start3[0], start3[1]);
      Point e0(end0[0], end0[1]);
      Point e1(end1[0], end1[1]);
      Point e2(end2[0], end2[1]);
      Point e3(end3[0], end3[1]);

      DrawLine(s0, e0, 255, res, w, h);
      DrawLine(s1, e1, 255, res, w, h);
      DrawLine(s2, e2, 255, res, w, h);
      DrawLine(s3, e3, 255, res, w, h);

      printf("(%.5e, %.5e, %.5e), (%.5e, %.5e, %.5e) \n", plane.getNormal()[0], plane.getNormal()[1], plane.getNormal()[2], plane.getPoint()[0], plane.getPoint()[1], plane.getPoint()[2]);
      printf("(%.5e, %.5e, %.5e) -> (%.5e, %.5e, %.5e)\n", res0.P0[0], res0.P0[1], res0.P0[2], res0.P1[0], res0.P1[1], res0.P1[2]);
      printf("(%.5e, %.5e, %.5e) -> (%.5e, %.5e, %.5e)\n", res1.P0[0], res1.P0[1], res1.P0[2], res1.P1[0], res1.P1[1], res1.P1[2]);
      printf("(%.5e, %.5e, %.5e) -> (%.5e, %.5e, %.5e)\n", res2.P0[0], res2.P0[1], res2.P0[2], res2.P1[0], res2.P1[1], res2.P1[2]);
      printf("(%.5e, %.5e, %.5e) -> (%.5e, %.5e, %.5e)\n", res3.P0[0], res3.P0[1], res3.P0[2], res3.P1[0], res3.P1[1], res3.P1[2]);

      printf("(%d, %d, %d) -> (%d, %d, %d)\n", start0[0], start0[1], start0[2], end0[0], end0[1], end0[2]);
      printf("(%d, %d, %d) -> (%d, %d, %d)\n", start1[0], start1[1], start1[2], end1[0], end1[1], end1[2]);
      printf("(%d, %d, %d) -> (%d, %d, %d)\n", start2[0], start2[1], start2[2], end2[0], end2[1], end2[2]);
      printf("(%d, %d, %d) -> (%d, %d, %d)\n", start3[0], start3[1], start3[2], end3[0], end3[1], end3[2]);
   }

   __host__ __device__ NormalIntersectionShapeT* Line::get()
   {
      return nis;
   }

   __host__ __device__ void TestProjection()
   {
      const double p[3] = { 
         -7., 
         -8., 
         -9. 
      };
      const double n[3] = {
         3. / ::sqrtf(50.f), 
         4. / ::sqrtf(50.f), 
         5. / ::sqrtf(50.f) 
      };
      const PlaneT plane(n, p); // projection plane
      const double axis0[3] = {
         -3. / ::sqrtf(50.f), 
         -4. / ::sqrtf(50.f), 
         5. / ::sqrtf(50.f) 
      }; // projection vector from projection plane origin
      const double axis1[3] = {
         4. / 5., 
         -3. / 5., 
         0. 
      }; // projection vector from projection plane origin

      const double p0[3] = {
         p[0] + axis0[0] * ::sqrtf(50.f) + n[0] * ::sqrtf(50.f), 
         p[1] + axis0[1] * ::sqrtf(50.f) + n[1] * ::sqrtf(50.f), 
         p[2] + axis0[2] * ::sqrtf(50.f) + n[2] * ::sqrtf(50.f) 
      };
      const double p1[3] = {
         p[0] + axis1[0] * 5. - n[0] * ::sqrtf(50.f), 
         p[1] + axis1[1] * 5. - n[1] * ::sqrtf(50.f), 
         p[2] + axis1[2] * 5. - n[2] * ::sqrtf(50.f) 
      };
      MultiPlaneShape::LineShape seg(p0, p1);
      MultiPlaneShape::LineShape res;
      calcLineProjection(plane, axis0, axis1, seg, res);

      static const float tol = 0.00001f;
      if (::fabsf(::sqrtf(50.f) - res.P0[0]) > tol) printf("res.P0[0] wrong: %.5e\n", res.P0[0]);
      if (::fabsf(res.P0[1]) > tol) printf("res.P0[1] wrong: %.5e\n", res.P0[1]);
      if (::fabsf(res.P0[2]) > tol) printf("res.P0[2] wrong: %.5e\n", res.P0[2]);
      if (::fabsf(res.P1[0]) > tol) printf("res.P1[0] wrong: %.5e\n", res.P1[0]);
      if (::fabsf(5.f - res.P1[1]) > tol) printf("res.P1[1] wrong: %.5e\n", res.P1[1]);
      if (::fabsf(res.P1[2]) > tol) printf("res.P1[2] wrong: %.5e\n", res.P1[2]);

      printf("(%.5e, %.5e, %.5e), (%.5e, %.5e, %.5e)", res.P0[0], res.P0[1], res.P0[2], res.P1[0], res.P1[1], res.P1[2]);
   }

   __host__ __device__ CrossSection::CrossSection(const double width)
   {
      enclosure = new NormalMultiPlaneShapeT();

      const double n0[] = { 0., 0., -1 }, p0[] = { 0., 0., 0 }; // Add top plane
      bnd = new PlaneT(n0, p0);
      enclosure->addPlane(bnd);

      const double n1[] = { 0., 1., 0. }, p1[] = { 0., width / 2., 0. }; // Right end
      pos = new PlaneT(n1, p1);
      enclosure->addPlane(pos);

      const double n2[] = { 0., -1., 0. }, p2[] = { 0., -width / 2., 0. }; // Left end
      double n2norm[3];
      Math2::normalize3d(n2, n2norm);
      neg = new PlaneT(n2norm, p2);
      enclosure->addPlane(neg);
   }

   __host__ __device__ CrossSection::~CrossSection()
   {
      delete enclosure;
      delete bnd;
      delete pos;
      delete neg;
   }

   __host__ __device__ NormalMultiPlaneShapeT* CrossSection::get()
   {
      return enclosure;
   }
}