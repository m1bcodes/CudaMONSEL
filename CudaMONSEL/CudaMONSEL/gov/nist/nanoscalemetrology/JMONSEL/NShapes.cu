#include "gov\nist\nanoscalemetrology\JMONSEL\NShapes.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalMultiPlaneShape.cuh"
#include "gov\nist\microanalysis\NISTMonte\MultiPlaneShape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalIntersectionShape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalCylindricalShape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalUnionShape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalComplementShape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalDifferenceShape.cuh"
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
      extensionVec[0] *= 10.;
      extensionVec[1] *= 10.;
      extensionVec[2] *= 10.;
      Math2::minus3d(segment.P0, extensionVec, segment.P0);
      Math2::plus3d(segment.P1, extensionVec, segment.P1);

      double rightIntersect[3];
      if (MultiPlaneShape::intersect3D_SegmentPlane(segment, rightEnd, rightIntersect) != 1) printf("createLineProjection: error 2\n");
      double leftIntersect[3];
      if (MultiPlaneShape::intersect3D_SegmentPlane(segment, leftEnd, leftIntersect) != 1) printf("createLineProjection: error 3\n");

      memcpy(segment.P0, rightIntersect, sizeof(segment.P0[0]) * 3);
      memcpy(segment.P1, leftIntersect, sizeof(segment.P1[0]) * 3);
   }

   __host__ __device__ Line::Line(const double topz, const double width, const double length, const double thetal, const double thetar, const double radtl, const double radtr) :
      topz(topz), width(width), length(length), thetal(thetal), thetar(thetar), radtl(radtl < 0 ? 0 : radtl), radbl(this->radtl), radtr(radtr < 0 ? 0 : radtr), radbr(this->radtr),
      enclosure(nullptr), pl0(nullptr), pl1(nullptr), pl2(nullptr), pl3(nullptr),
      // right
      rightSide(nullptr), pl4(nullptr),
      toprightNMPS(nullptr), toprightSide(nullptr), pltr(nullptr), trcylinder(nullptr),
      botrightNMPS(nullptr), plbr(nullptr), brcylinder(nullptr), plbr_0(nullptr), plbr_1(nullptr), botrightSide(nullptr),
      // left
      leftSide(nullptr), pl5(nullptr),
      topleftNMPS(nullptr), topleftSide(nullptr), pltl(nullptr), tlcylinder(nullptr),
      botleftNMPS(nullptr), plbl(nullptr), blcylinder(nullptr), plbl_0(nullptr), plbl_1(nullptr), botleftSide(nullptr),
      // encl
      nts(nullptr), nis(nullptr),
      // gt
      gt0(nullptr), gt1(nullptr), gt2(nullptr), gt3(nullptr),
      // others
      plrestr0(nullptr)//, pl7(nullptr)
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

      delete this->plrestr0; this->plrestr0 = nullptr;
      //delete this->pl7; this->pl7 = nullptr;

      // top
      // right
      delete this->toprightNMPS; this->toprightNMPS = nullptr;
      if (this->pltr) {
         delete this->pltr;
         this->pltr = nullptr;
      }
      if (this->trcylinder) {
         delete this->trcylinder;
         this->trcylinder = nullptr;
      }
      if (this->toprightSide) {
         delete this->toprightSide;
         this->toprightSide = nullptr;
      }
      // bot
      delete this->botrightNMPS; this->botrightNMPS = nullptr;
      if (this->plbr) {
         delete this->plbr;
         this->plbr = nullptr;
      }
      if (this->plbr_0) {
         delete this->plbr_0;
         this->plbr_0 = nullptr;
      }
      if (this->plbr_1) {
         delete this->plbr_1;
         this->plbr_1 = nullptr;
      }
      if (this->brcylinder) {
         delete this->brcylinder;
         this->brcylinder = nullptr;
      }
      if (this->botrightSide) {
         delete this->botrightSide;
         this->botrightSide = nullptr;
      }
      // others
      if (this->pl4) {
         delete this->pl4;
         this->pl4 = nullptr;
      }
      if (this->rightSide) {
         delete (NormalUnionShapeT*)this->rightSide;
         this->rightSide = nullptr;
      }

      // left
      // top
      delete this->topleftNMPS; this->topleftNMPS = nullptr;
      if (this->pltl) {
         delete this->pltl;
         this->pltl = nullptr;
      }
      if (this->tlcylinder) {
         delete this->tlcylinder;
         this->tlcylinder = nullptr;
      }
      if (this->topleftSide) {
         delete this->topleftSide;
         this->topleftSide = nullptr;
      }
      // bot
      delete this->botleftNMPS; this->botleftNMPS = nullptr;
      if (this->plbl) {
         delete this->plbl;
         this->plbl = nullptr;
      }
      if (this->plbl_0) {
         delete this->plbl_0;
         this->plbl_0 = nullptr;
      }
      if (this->plbl_1) {
         delete this->plbl_1;
         this->plbl_1 = nullptr;
      }
      if (this->blcylinder) {
         delete this->blcylinder;
         this->blcylinder = nullptr;
      }
      if (this->botleftSide) {
         delete this->botleftSide;
         this->botleftSide = nullptr;
      }
      // others
      if (this->pl5) {
         delete this->pl5;
         this->pl5 = nullptr;
      }
      if (this->leftSide) {
         delete (NormalUnionShapeT*)this->leftSide;
         this->leftSide = nullptr;
      }

      // encl
      delete this->nts; this->nts = nullptr;
      delete this->nis; this->nis = nullptr;

      // gt
      if (gt0) { delete gt0; gt0 = nullptr; }
      if (gt1) { delete gt1; gt1 = nullptr; }
      if (gt2) { delete gt2; gt2 = nullptr; }
      if (gt3) { delete gt3; gt3 = nullptr; }
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

      toprightNMPS = new NormalMultiPlaneShapeT();
      rightSide = nullptr;

      // Add right sidewall
      const double costhetar = ::cos(thetar);
      const double sinthetar = ::sin(thetar);
      const double tanthetar = sinthetar / costhetar;
      const double n4[] = { costhetar, 0., signz * sinthetar }, p4[] = { width / 2., 0., 0. };
      pl4 = new PlaneT(n4, p4);
      toprightNMPS->addPlane(pl4);
      // If radtr>0 add a clipping plane and the cylinder
      const double root2 = ::sqrt(2.);
      const double absz = signz * topz;
      if (radtr > 0) {
         const double rad = ::sqrt(1. - sinthetar);

         // top curve
         //const double xc = ((width / 2.) - (radtr / costhetar)) + ((radtr - absz) * ::tan(thetar));
         const double xtc = (width / 2.) - (radtr * costhetar) + (radtr - absz) * tanthetar;
         //const double ntr[] = { rad / root2, 0., signz * costhetar / root2 / rad }, ptr[] = { width / 2. - radtr / costhetar + (radtr - absz) * sinthetar / costhetar, 0., topz };
         const double ntr[] = { rad / root2, 0., signz * costhetar / root2 / rad }, ptr[] = { xtc, 0., topz };
         pltr = new PlaneT(ntr, ptr); // slicer
         toprightNMPS->addPlane(pltr);
         // Construct cylinder for right corner
         const double ztc = topz - (signz * radtr);
         const double end0tr[] = { xtc, -length / 2., ztc }, end1tr[] = { xtc, length / 2., ztc };
         trcylinder = new NormalCylindricalShapeT(end0tr, end1tr, radtr);
         toprightSide = new NormalUnionShapeT(*toprightNMPS, *trcylinder);
         //rightSide = new NormalUnionShapeT(*toprightNMPS, *trcylinder);

         // bottom curve
         botrightNMPS = new NormalMultiPlaneShapeT();
         const double nbr_0[] = { 0., 0., -signz }, pbr_0[] = { 0., 0., 0. };
         plbr_0 = new PlaneT(nbr_0, pbr_0); // z = 0 plane
         botrightNMPS->addPlane(plbr_0);

         //const double nbr_1[] = { -n4[0], -n4[1], -n4[2] }, pbr_1[] = { p4[0], p4[1], p4[2] };
         //plbr_1 = new PlaneT(nbr_1, pbr_1);
         //botrightNMPS->addPlane(plbr_1);

         const double xbc = width / 2. + radbr * (1. - tanthetar);
         //const double nbr[] = { rad / root2, 0., signz * costhetar / root2 / rad }, pbr[] = { xbc, 0., 0. };
         const double nbr[] = { ntr[0], ntr[1], ntr[2] }, pbr[] = { xbc, 0., 0. };
         plbr = new PlaneT(nbr, pbr); // shifted pltr (slicer)
         botrightNMPS->addPlane(plbr);

         const double zbc = signz * radbr;
         const double end0br[] = { xbc, -length / 2., zbc }, end1br[] = { xbc, length / 2., zbc };
         brcylinder = new NormalCylindricalShapeT(end0br, end1br, radbr);

         botrightSide = new NormalDifferenceShapeT(*botrightNMPS, *brcylinder);

         // union
         rightSide = new NormalUnionShapeT(*toprightSide, *botrightSide);
      }
      else
         rightSide = toprightNMPS;

      topleftNMPS = new NormalMultiPlaneShapeT();
      leftSide = nullptr;

      // Add left sidewall
      const double costhetal = ::cos(thetal);
      const double sinthetal = ::sin(thetal);
      const double tanthetal = sinthetal / costhetal;
      const double n5[] = { -costhetal, 0., signz * sinthetal }, p5[] = { -width / 2., 0., 0. };
      pl5 = new PlaneT(n5, p5);
      topleftNMPS->addPlane(pl5);
      // If radtl>0 add a clipping plane and the cylinder
      if (radtl > 0.) {
         const double rad = ::sqrt(1. - sinthetal);

         // top curve
         const double xtc = -(width / 2. - radtl * costhetal + (radtl - absz) * tanthetal);
         //const double nl[] = { -rad / root2, 0., signz * costhetal / root2 / rad }, pl[] = { -width / 2. + radtl / costhetal - (radtl - absz) * sinthetal / costhetal, 0., topz };
         const double ntl[] = { -rad / root2, 0., signz * costhetal / root2 / rad }, ptl[] = { xtc, 0., topz };
         pltl = new PlaneT(ntl, ptl);
         topleftNMPS->addPlane(pltl);
         //const double xc = width / 2. - radtl / ::cos(thetal) + (radtl - absz) * ::tan(thetal);
         const double ztc = topz - (signz * radtl);
         // Construct cylinder for left corner
         const double end0[] = { xtc, -length / 2., ztc }, end1[] = { xtc, length / 2., ztc };
         tlcylinder = new NormalCylindricalShapeT(end0, end1, radtl);
         topleftSide = new NormalUnionShapeT(*topleftNMPS, *tlcylinder);
         //leftSide = new NormalUnionShapeT(*topleftNMPS, *tlcylinder);

         // bottom curve
         botleftNMPS = new NormalMultiPlaneShapeT();
         const double nbl_0[] = { 0., 0., -signz }, pbl_0[] = { 0., 0., 0. };
         plbl_0 = new PlaneT(nbl_0, pbl_0); // z = 0 plane
         botleftNMPS->addPlane(plbl_0);

         //const double nbr_1[] = { -n4[0], -n4[1], -n4[2] }, pbr_1[] = { p4[0], p4[1], p4[2] };
         //plbr_1 = new PlaneT(nbr_1, pbr_1);
         //botrightNMPS->addPlane(plbr_1);

         const double xbc = -(width / 2. + radbl * (1. - tanthetal));
         //const double nbr[] = { rad / root2, 0., signz * costhetar / root2 / rad }, pbr[] = { xbc, 0., 0. };
         const double nbl[] = { ntl[0], ntl[1], ntl[2] }, pbl[] = { xbc, 0., 0. };
         plbl = new PlaneT(nbl, pbl); // shifted pltl (slicer)
         botleftNMPS->addPlane(plbl);

         const double zbc = signz * radbl;
         const double end0bl[] = { xbc, -length / 2., zbc }, end1bl[] = { xbc, length / 2., zbc };
         blcylinder = new NormalCylindricalShapeT(end0bl, end1bl, radbl);

         botleftSide = new NormalDifferenceShapeT(*botleftNMPS, *blcylinder);

         // union
         leftSide = new NormalUnionShapeT(*topleftSide, *botleftSide);
      }
      else
         leftSide = topleftNMPS;

      nts = new NormalIntersectionShapeT(*leftSide, *rightSide);
      nis = new NormalIntersectionShapeT(*nts, *enclosure);
   }

   __host__ __device__ void Line::calcGroundtruth()
   {
      if (!gt0) gt0 = new LineShapeT();
      if (!gt1) gt1 = new LineShapeT();
      if (!gt2) gt2 = new LineShapeT();
      if (!gt3) gt3 = new LineShapeT();

      // top view
      //getLineSegment(*pl0, *pl2, *pl4, *pl5top, *segment0); // right cap
      //getLineSegment(*pl0, *pl3, *pl4, *pl5top, *gt1); // left cap
      //getLineSegment(*pl0, *pl4, *pl2, *pl3, *gt2); // right length
      //getLineSegment(*pl0, *pl5top, *pl2, *pl3, *gt3); // left length

      // side view (back)
      //{
      //   if (!pltl) {
      //      getLineSegment(*pl3, *pl5top, *pl0, *pl1, *gt1); // left cap
      //      getLineSegment(*pl3, *pl0, *pl4, *pl5top, *gt3); // top length
      //   }
      //   else {
      //      getLineSegment(*pl3, *pl5top, *pltl, *pl1, *gt1); // left cap
      //      getLineSegment(*pl3, *pl0, *pltr, *pltl, *gt3); // top length
      //   }
      //   if (!pltr) {
      //      getLineSegment(*pl3, *pl4, *pl0, *pl1, *gt0); // right cap
      //      getLineSegment(*pl3, *pl1, *pl4, *pl5top, *gt2); // bottom length
      //   }
      //   else {
      //      getLineSegment(*pl3, *pl4, *pltr, *pl1, *gt0); // right cap
      //      getLineSegment(*pl3, *pl1, *pl4, *pl5top, *gt2); // bottom length
      //   }
      //}

      // side view (front)
      {
         // without bottom curve
         //if (!pltl) {
         //   getLineSegment(*pl2, *pl5, *pl0, *pl1, *gt1); // left cap
         //   getLineSegment(*pl2, *pl0, *pl4, *pl5, *gt3); // top length
         //}
         //else {
         //   getLineSegment(*pl2, *pl5, *pltl, *pl1, *gt1); // left cap
         //   getLineSegment(*pl2, *pl0, *pltr, *pltl, *gt3); // top length
         //}
         //if (!pltr) {
         //   getLineSegment(*pl2, *pl4, *pl0, *pl1, *gt0); // right cap
         //   getLineSegment(*pl2, *pl1, *pl4, *pl5, *gt2); // bottom length
         //}
         //else {
         //   getLineSegment(*pl2, *pl4, *pltr, *pl1, *gt0); // right cap
         //   getLineSegment(*pl2, *pl1, *pl4, *pl5, *gt2); // bottom length
         //}

         // with bottom curve
         // 0 <-> top
         // 1 <-> bottom
         // 4 <-> right
         // 5 <-> left
         getLineSegment(*pl2, *pl5, pltl ? *pltl : *pl0, plbl ? *plbl : *pl1, *gt1); // left cap (second arg = 5 <-> l)
         getLineSegment(*pl2, *pl0, pltr ? *pltr : *pl4, pltl ? *pltl : *pl5, *gt3); // top length (second arg = 0 <-> t)
         getLineSegment(*pl2, *pl4, pltr ? *pltr : *pl0, plbr ? *plbr : *pl1, *gt0); // right cap (second arg = 4 <-> r)
         getLineSegment(*pl2, *pl1, plbr ? *plbr : *pl4, plbl ? *plbl : *pl5, *gt2); // bottom length (second arg = 1 <-> b)
      }

      //printf("(%.5e, %.5e, %.5e) -> (%.5e, %.5e, %.5e)\n", gt0->P0[0], gt0->P0[1], gt0->P0[2], gt0->P1[0], gt0->P1[1], gt0->P1[2]);
      //printf("(%.5e, %.5e, %.5e) -> (%.5e, %.5e, %.5e)\n", gt1->P0[0], gt1->P0[1], gt1->P0[2], gt1->P1[0], gt1->P1[1], gt1->P1[2]);
      //printf("(%.5e, %.5e, %.5e) -> (%.5e, %.5e, %.5e)\n", gt2->P0[0], gt2->P0[1], gt2->P0[2], gt2->P1[0], gt2->P1[1], gt2->P1[2]);
      //printf("(%.5e, %.5e, %.5e) -> (%.5e, %.5e, %.5e)\n", gt3->P0[0], gt3->P0[1], gt3->P0[2], gt3->P1[0], gt3->P1[1], gt3->P1[2]);
   }

   /*
   * assumes perpendicular projection onto the plane
   */
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
      if (len - 1.f > .000001f) printf("Line::getPerpendicularIntersection: plane normal length not unity %.5e\n", len);

      // calc
      double v[3];
      Math2::minus3d(p, p0, v);
      double s = Math2::dot3d(v, n0); // perpendicular distance between p and plane
      s = s == 0. ? 1. : s;
      MultiPlaneShape::LineShape l;
      static const float extensionFactor = 20.f;
      l.P0[0] = p[0] + extensionFactor * s * n0[0]; l.P0[1] = p[1] + extensionFactor * s * n0[1]; l.P0[2] = p[2] + extensionFactor * s * n0[2];
      l.P1[0] = p[0] - extensionFactor * s * n0[0]; l.P1[1] = p[1] - extensionFactor * s * n0[1]; l.P1[2] = p[2] - extensionFactor * s * n0[2];
      const int ind = MultiPlaneShape::intersect3D_SegmentPlane(l, plane, res);

      // checks
      if (ind == 0) printf("getPerpendicularIntersection: no intersect (check if decreasing MultiPlaneShape::SMALL_NUM helps)\n");
      else if (ind == 1) printf("getPerpendicularIntersection: intersect\n");
      else if (ind == 2) printf("getPerpendicularIntersection: on plane\n");
   }

   /*
   * note that the projdir and the plane normal is assumed to point in 'opposite' directions (cross a plane in opposite directions)
   */
   __host__ __device__ static void getIntersection(
       const PlaneT& plane, // the projection plane
       const double* p, // absolute coordinate to be projected
       const double* projdir, // absolute projection direction from p
       double* res // the absolute intersection coordinate between p and plane (lies on the plane)
       )
   {
       // checks
       const double* n0 = plane.getNormal();
       const float len = ::sqrtf(Math2::dot3d(n0, n0));
       if (Math2::fabs(len - 1.f) > .000001f) printf("Line::getIntersection: plane normal length not unity %.5e\n", len);

       const float projdirlen = ::sqrtf(Math2::dot3d(projdir, projdir));
       if (Math2::fabs(projdirlen - 1.f) > .000001f) printf("Line::getIntersection: direction vector length not unity %.5e\n", projdirlen);

       // calc
       double v[3];
       const double* p0 = plane.getPoint();
       Math2::minus3d(p, p0, v);
       double s = Math2::dot3d(v, n0); // perpendicular distance between p and plane
       s = s == 0. ? 1. : s;
       //const double d = s / -Math2::dot3d(v, n0); // length of the hypothenuse
       MultiPlaneShape::LineShape l;
       static const float extensionFactor = 50.f;
       l.P0[0] = p[0] + extensionFactor * s * projdir[0]; l.P0[1] = p[1] + extensionFactor * s * projdir[1]; l.P0[2] = p[2] + extensionFactor * s * projdir[2]; // WARNING: using s (ie the perpendicular ditance) is wrong, need hypothenuse
       l.P1[0] = p[0] - extensionFactor * s * projdir[0]; l.P1[1] = p[1] - extensionFactor * s * projdir[1]; l.P1[2] = p[2] - extensionFactor * s * projdir[2]; // WARNING: using s (ie the perpendicular ditance) is wrong, need hypothenuse
       const int ind = MultiPlaneShape::intersect3D_SegmentPlane(l, plane, res);

       // checks
       if (ind == 0) printf("getPerpendicularIntersection: no intersect (check if decreasing MultiPlaneShape::SMALL_NUM helps)\n");
       else if (ind == 1) printf("getPerpendicularIntersection: intersect\n");
       else if (ind == 2) printf("getPerpendicularIntersection: on plane\n");
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

   __host__ __device__ static void calcPointProjection(
      const PlaneT& plane, // plane on which point is projected
      const double* axis0, // vector wrt the plane origin (plane coordinate)
      const double* axis1, // vector wrt the plane origin (plane coordinate)
      const double* lightdir, // abosolute vector of where the light source is pointing
      const double line[], // line in absolute coordinates to be projected onto plane
      double res[] // output: line with respect to origin of plane (ie. res + plane.origin = line)
      )
   {
      double tmp[3]; // required
      //getPerpendicularIntersection(plane, line, tmp);
      const double reverseddir[3] = { -lightdir[0], -lightdir[1], -lightdir[2] };
      getIntersection(plane, line, reverseddir, tmp);
      getRelativeCoord(plane, axis0, axis1, tmp, res);
   }

   __host__ __device__ static void calcLineProjection(
      const PlaneT& plane,
      const double* axis0, // vector from the plane origin
      const double* axis1, // vector from the plane origin
      const double* lightdir,
      const MultiPlaneShape::LineShape& line, // line in absolute coordinates to be projected onto plane
      MultiPlaneShape::LineShape& res // line with respect to origin of plane (ie. res + plane.origin = line)
      )
   {
      calcPointProjection(plane, axis0, axis1, lightdir, line.P0, res.P0);
      calcPointProjection(plane, axis0, axis1, lightdir, line.P1, res.P1);
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

   __host__ __device__ void drawRaster(
      double start[3],
      double end[3],
      const PlaneT& plane,
      const double* axis0, // vector from the plane origin
      const double* axis1, // vector from the plane origin
      const double* lightdir,
      const float xlenperpix,
      const float ylenperpix,
      char* res,
      const unsigned int w,
      const unsigned int h
      )
   {
      LineShapeT l(start, end);
      MultiPlaneShape::LineShape lr;
      calcLineProjection(plane, axis0, axis1, lightdir, l, lr); // TODO: 
      const int s4[3] = { lr.P0[0] / xlenperpix, lr.P0[1] / ylenperpix, lr.P0[2] };
      const int e4[3] = { lr.P1[0] / xlenperpix, lr.P1[1] / ylenperpix, lr.P1[2] };
      const Point c0(s4[0], s4[1]);
      const Point c1(e4[0], e4[1]);
      DrawLine(c0, c1, 255, res, w, h);
   }

   /*
   * project from point to light source in reverse lightdir
   */
   // axis2 would be the normal of the plane, which creates no projection
   __host__ __device__ void Line::calcRasterization(
      const PlaneT& plane, // plane where the light source is on/will be moving on
      const double* axis0, // vector from the plane origin
      const double* axis1, // vector from the plane origin
      const double* lightdir,
      const float xlenperpix,
      const float ylenperpix,
      char* res,
      const unsigned int w,
      const unsigned int h
      ) const
   {
      MultiPlaneShape::LineShape res0, res1, res2, res3; // with respect to origin of plane, along axis0 and axis1
      calcLineProjection(plane, axis0, axis1, lightdir, *gt0, res0); // TODO: 
      calcLineProjection(plane, axis0, axis1, lightdir, *gt1, res1); // TODO: 
      calcLineProjection(plane, axis0, axis1, lightdir, *gt2, res2); // TODO: 
      calcLineProjection(plane, axis0, axis1, lightdir, *gt3, res3); // TODO: 

      const int start0[3] = { res0.P0[0] / xlenperpix, res0.P0[1] / ylenperpix, res0.P0[2] };
      const int end0[3] = { res0.P1[0] / xlenperpix, res0.P1[1] / ylenperpix, res0.P1[2] };
      const int start1[3] = { res1.P0[0] / xlenperpix, res1.P0[1] / ylenperpix, res1.P0[2] };
      const int end1[3] = { res1.P1[0] / xlenperpix, res1.P1[1] / ylenperpix, res1.P1[2] };
      const int start2[3] = { res2.P0[0] / xlenperpix, res2.P0[1] / ylenperpix, res2.P0[2] };
      const int end2[3] = { res2.P1[0] / xlenperpix, res2.P1[1] / ylenperpix, res2.P1[2] };
      const int start3[3] = { res3.P0[0] / xlenperpix, res3.P0[1] / ylenperpix, res3.P0[2] };
      const int end3[3] = { res3.P1[0] / xlenperpix, res3.P1[1] / ylenperpix, res3.P1[2] };

      const Point s0(start0[0], start0[1]);
      const Point s1(start1[0], start1[1]);
      const Point s2(start2[0], start2[1]);
      const Point s3(start3[0], start3[1]);
      const Point e0(end0[0], end0[1]);
      const Point e1(end1[0], end1[1]);
      const Point e2(end2[0], end2[1]);
      const Point e3(end3[0], end3[1]);

      DrawLine(s0, e0, 255, res, w, h);
      DrawLine(s1, e1, 255, res, w, h);
      DrawLine(s2, e2, 255, res, w, h); // assuming same material as top layer
      DrawLine(s3, e3, 255, res, w, h);

      // top left
      if (tlcylinder) {
         // crosshair
         //double cstart0[3] = { lcylinder->getEnd0()[0], lcylinder->getEnd0()[1], lcylinder->getEnd0()[2] };
         //double cend0[3] = { lcylinder->getEnd0()[0] + lcylinder->getRadius(), lcylinder->getEnd0()[1], lcylinder->getEnd0()[2] };
         //drawRaster(cstart0, cend0, plane, axis0, axis1, xlenperpix, ylenperpix, res, w, h);

         //double cstart1[3] = { lcylinder->getEnd0()[0], lcylinder->getEnd0()[1], lcylinder->getEnd0()[2] };
         //double cend1[3] = { lcylinder->getEnd0()[0], lcylinder->getEnd0()[1] + lcylinder->getRadius(), lcylinder->getEnd0()[2] };
         //drawRaster(cstart1, cend1, plane, axis0, axis1, xlenperpix, ylenperpix, res, w, h);

         //double cstart2[3] = { lcylinder->getEnd0()[0], lcylinder->getEnd0()[1], lcylinder->getEnd0()[2] };
         //double cend2[3] = { lcylinder->getEnd0()[0] - lcylinder->getRadius(), lcylinder->getEnd0()[1], lcylinder->getEnd0()[2] };
         //drawRaster(cstart2, cend2, plane, axis0, axis1, xlenperpix, ylenperpix, res, w, h);

         //double cstart3[3] = { lcylinder->getEnd0()[0], lcylinder->getEnd0()[1], lcylinder->getEnd0()[2] };
         //double cend3[3] = { lcylinder->getEnd0()[0], lcylinder->getEnd0()[1] - lcylinder->getRadius(), lcylinder->getEnd0()[2] };
         //drawRaster(cstart3, cend3, plane, axis0, axis1, xlenperpix, ylenperpix, res, w, h);

         //LineShapeT clipper0, clipper0Proj;
         //getLineSegment(*pl3, *pltl, *pl0, *pl5top, clipper0);
         //calcLineProjection(plane, axis0, axis1, clipper0, clipper0Proj);
         //const int clstart0[3] = { clipper0Proj.P0[0] / xlenperpix, clipper0Proj.P0[1] / ylenperpix, clipper0Proj.P0[2] };
         //const int clend0[3] = { clipper0Proj.P1[0] / xlenperpix, clipper0Proj.P1[1] / ylenperpix, clipper0Proj.P1[2] };
         //const Point cls0(clstart0[0], clstart0[1]);
         //const Point cle0(clend0[0], clend0[1]);
         //DrawLine(cls0, cle0, 255, res, w, h);

         // points on the quarter circle for top left corner curve, uncomment the code above to see why this is valid
         for (int i = 180; i < 270; ++i) {
            const double theta = i / 180. * Math2::PI;
            double orig[3] = {
               tlcylinder->getEnd1()[0] + tlcylinder->getRadius() * ::cosf(theta),
               tlcylinder->getEnd1()[1] + tlcylinder->getRadius() * ::sinf(theta),
               tlcylinder->getEnd1()[2]
            }; // eg. point on the circle
            double proj[3];
            calcPointProjection(plane, axis0, axis1, lightdir, orig, proj);
            const int s[3] = { proj[0] / xlenperpix, proj[1] / ylenperpix, proj[2] };
            const unsigned int idx = s[1] * w + s[0];
            if (idx >= 0 && idx < h * w) {
               if (s[1] >= 0 && s[1] < h && s[0] >= 0 && s[0] < w)
                  res[idx] = 255;
            }
            else printf("Line::calcRasterization: out of bound index top left: %d, %d\n", s[0], s[1]);
         }
      }
      // top right
      if (trcylinder) {
         // crosshair
         //double cstart0[3] = { rcylinder->getEnd0()[0], rcylinder->getEnd0()[1], rcylinder->getEnd0()[2] };
         //double cend0[3] = { rcylinder->getEnd0()[0] + rcylinder->getRadius(), rcylinder->getEnd0()[1], rcylinder->getEnd0()[2] };
         //drawRaster(cstart0, cend0, plane, axis0, axis1, xlenperpix, ylenperpix, res, w, h);

         //double cstart1[3] = { rcylinder->getEnd0()[0], rcylinder->getEnd0()[1], rcylinder->getEnd0()[2] };
         //double cend1[3] = { rcylinder->getEnd0()[0], rcylinder->getEnd0()[1] + rcylinder->getRadius(), rcylinder->getEnd0()[2] };
         //drawRaster(cstart1, cend1, plane, axis0, axis1, xlenperpix, ylenperpix, res, w, h);

         //double cstart2[3] = { rcylinder->getEnd0()[0], rcylinder->getEnd0()[1], rcylinder->getEnd0()[2] };
         //double cend2[3] = { rcylinder->getEnd0()[0] - rcylinder->getRadius(), rcylinder->getEnd0()[1], rcylinder->getEnd0()[2] };
         //drawRaster(cstart2, cend2, plane, axis0, axis1, xlenperpix, ylenperpix, res, w, h);

         //double cstart3[3] = { rcylinder->getEnd0()[0], rcylinder->getEnd0()[1], rcylinder->getEnd0()[2] };
         //double cend3[3] = { rcylinder->getEnd0()[0], rcylinder->getEnd0()[1] - rcylinder->getRadius(), rcylinder->getEnd0()[2] };
         //drawRaster(cstart3, cend3, plane, axis0, axis1, xlenperpix, ylenperpix, res, w, h);

         //LineShapeT clipper0, clipper0Proj;
         //getLineSegment(*pl3, *pltr, *pl0, *pl4, clipper0);
         //calcLineProjection(plane, axis0, axis1, clipper0, clipper0Proj);
         //const int clstart0[3] = { clipper0Proj.P0[0] / xlenperpix, clipper0Proj.P0[1] / ylenperpix, clipper0Proj.P0[2] };
         //const int clend0[3] = { clipper0Proj.P1[0] / xlenperpix, clipper0Proj.P1[1] / ylenperpix, clipper0Proj.P1[2] };
         //const Point cls0(clstart0[0], clstart0[1]);
         //const Point cle1(clend0[0], clend0[1]);
         //DrawLine(cls0, cle1, 255, res, w, h);

         // points on the quarter circle for top right corner curve, uncomment the code above to see why this is valid
         for (int i = 270; i < 360; ++i) {
            const double theta = i / 180. * Math2::PI;
            const double orig[3] = {
               trcylinder->getEnd1()[0] + trcylinder->getRadius() * ::cosf(theta),
               trcylinder->getEnd1()[1] + trcylinder->getRadius() * ::sinf(theta),
               trcylinder->getEnd1()[2]
            }; // eg. point on the circle
            double proj[3];
            calcPointProjection(plane, axis0, axis1, lightdir, orig, proj);
            const int s[3] = { proj[0] / xlenperpix, proj[1] / ylenperpix, proj[2] };
            const unsigned int idx = s[1] * w + s[0];
            if (idx >= 0 && idx < h * w) {
               if (s[1] >= 0 && s[1] < h && s[0] >= 0 && s[0] < w)
                  res[idx] = 255;
            }
            else printf("Line::calcRasterization: out of bound index top right: %d, %d\n", s[0], s[1]);
         }
      }
      // bottom left
      if (blcylinder) {
         for (int i = 0; i < 90; ++i) {
            const double theta = i / 180. * Math2::PI;
            double orig[3] = {
               blcylinder->getEnd1()[0] + blcylinder->getRadius() * ::cosf(theta),
               blcylinder->getEnd1()[1] + blcylinder->getRadius() * ::sinf(theta),
               blcylinder->getEnd1()[2]
            }; // eg. point on the circle
            double proj[3];
            calcPointProjection(plane, axis0, axis1, lightdir, orig, proj);
            const int s[3] = { proj[0] / xlenperpix, proj[1] / ylenperpix, proj[2] };
            const unsigned int idx = s[1] * w + s[0];
            if (idx >= 0 && idx < h * w) {
               if (s[1] >= 0 && s[1] < h && s[0] >= 0 && s[0] < w)
                  res[idx] = 255;
            }
            else printf("Line::calcRasterization: out of bound index bot left: %d, %d\n", s[0], s[1]);
         }
      }
      // bottom right
      if (brcylinder) {
         for (int i = 90; i < 180; ++i) {
            const double theta = i / 180. * Math2::PI;
            const double orig[3] = {
               brcylinder->getEnd1()[0] + brcylinder->getRadius() * ::cosf(theta),
               brcylinder->getEnd1()[1] + brcylinder->getRadius() * ::sinf(theta),
               brcylinder->getEnd1()[2]
            }; // eg. point on the circle
            double proj[3];
            calcPointProjection(plane, axis0, axis1, lightdir, orig, proj);
            const int s[3] = { proj[0] / xlenperpix, proj[1] / ylenperpix, proj[2] };
            const unsigned int idx = s[1] * w + s[0];
            if (idx >= 0 && idx < h * w) {
               if (s[1] >= 0 && s[1] < h && s[0] >= 0 && s[0] < w)
                  res[idx] = 255;
            }
            else printf("Line::calcRasterization: out of bound index bot right: %d, %d\n", s[0], s[1]);
         }
      }

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

   __host__ __device__ unsigned int Line::calcRasterPoint(
      const PlaneT& plane,
      const double* axis0, // vector from the plane origin
      const double* axis1, // vector from the plane origin
      const double* lightdir, 
      const float xlenperpix,
      const float ylenperpix,
      const unsigned int w,
      const double pt[],
      const double theta) const
   {
      double proj[3];
      calcPointProjection(plane, axis0, axis1, lightdir, pt, proj);
      const int s[3] = { proj[0] / xlenperpix, proj[1] / ylenperpix, proj[2] };
      return s[1] * w + s[0];
   }

   // axis2 would be the normal of the plane, which creates no projection
   __host__ __device__ void Line::calcRasterizationCorrection(
      const PlaneT& plane, // plane where the light source is on
      const double* axis0, // vector from the plane origin
      const double* axis1, // vector from the plane origin
      const double* lightdir,
      const float xlenperpix,
      const float ylenperpix,
      char* res,
      const unsigned int w,
      const unsigned int h
      ) const
   {
      MultiPlaneShape::LineShape res2;
      calcLineProjection(plane, axis0, axis1, lightdir, *gt2, res2); // TODO: 

      const int start2[3] = { res2.P0[0] / xlenperpix, res2.P0[1] / ylenperpix, res2.P0[2] };
      const int end2[3] = { res2.P1[0] / xlenperpix, res2.P1[1] / ylenperpix, res2.P1[2] };

      const Point s2(start2[0], start2[1]);
      const Point e2(end2[0], end2[1]);

      DrawLine(s2, e2, 0, res, w, h);

      // bottom left
      if (blcylinder) {
         for (int i = 0; i <= 90; ++i) {
            const double theta = i / 180. * Math2::PI;

            double pt0[3] = {
               blcylinder->getEnd1()[0] + blcylinder->getRadius() * ::cosf(theta),
               blcylinder->getEnd1()[1] + blcylinder->getRadius() * ::sinf(theta),
               blcylinder->getEnd1()[2]
            }; // eg. point on the circle
            const unsigned int idx = calcRasterPoint(plane, axis0, axis1, lightdir, xlenperpix, ylenperpix, w, pt0, theta);
            
            double pt1[3] = {
               blcylinder->getEnd1()[0] + blcylinder->getRadius() * ::cosf(theta),
               blcylinder->getEnd1()[1] + blcylinder->getRadius(),
               blcylinder->getEnd1()[2]
            }; // eg. point on the boundary
            const unsigned int idx1 = calcRasterPoint(plane, axis0, axis1, lightdir, xlenperpix, ylenperpix, w, pt1, theta);
            
            if (idx >= 0 && idx < h * w) { // valid index
               if (idx == idx1) { // on the boundary
                  res[idx1] = 255;
               }
            }
         }
      }
      // bottom right
      if (brcylinder) {
         for (int i = 90; i <= 180; ++i) {
            const double theta = i / 180. * Math2::PI;
            const double pt0[3] = {
               brcylinder->getEnd1()[0] + brcylinder->getRadius() * ::cosf(theta),
               brcylinder->getEnd1()[1] + brcylinder->getRadius() * ::sinf(theta),
               brcylinder->getEnd1()[2]
            }; // eg. point on the circle
            const unsigned int idx = calcRasterPoint(plane, axis0, axis1, lightdir, xlenperpix, ylenperpix, w, pt0, theta);
            
            double pt1[3] = {
               brcylinder->getEnd1()[0] + brcylinder->getRadius() * ::cosf(theta),
               brcylinder->getEnd1()[1] + brcylinder->getRadius(),
               brcylinder->getEnd1()[2]
            }; // eg. point on the boundary
            const unsigned int idx1 = calcRasterPoint(plane, axis0, axis1, lightdir, xlenperpix, ylenperpix, w, pt1, theta);
            
            if (idx >= 0 && idx < h * w) { // valid index
               if (idx == idx1) { // on the boundary
                  res[idx1] = 255;
               }
            }
         }
      }
   }

   __host__ __device__ NormalIntersectionShapeT* Line::get()
   {
      return nis;
   }

   __host__ __device__ void Line::addRestrainingPlanes()
   {
      const double nrestr0[] = { 0., 0., -1. }, prestr0[] = { 0., 0., 0. }; // Add top plane
      plrestr0 = new PlaneT(nrestr0, prestr0);
      enclosure->addPlane(plrestr0);

      //const double nrestr1[] = { 0., 1., 0. }, prestr1[] = { 0., 0., 0 }; // Add top plane
      //plrestr1 = new PlaneT(nrestr1, prestr1);
      //enclosure->addPlane(plrestr1);
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

      const double st = ::sinf(Math2::toRadians(0.));
      const double beamdir[3] = {
         ::cosf(Math2::toRadians(0.)) * st,
         ::sinf(Math2::toRadians(0.)) * st,
         ::cosf(Math2::toRadians(1.))
      };

      calcLineProjection(plane, axis0, axis1, beamdir, seg, res); // TODO: 

      static const float tol = 0.00001f;
      if (::fabsf(::sqrtf(50.f) - res.P0[0]) > tol) printf("res.P0[0] wrong: %.5e\n", res.P0[0]);
      if (::fabsf(res.P0[1]) > tol) printf("res.P0[1] wrong: %.5e\n", res.P0[1]);
      if (::fabsf(res.P0[2]) > tol) printf("res.P0[2] wrong: %.5e\n", res.P0[2]);
      if (::fabsf(res.P1[0]) > tol) printf("res.P1[0] wrong: %.5e\n", res.P1[0]);
      if (::fabsf(5.f - res.P1[1]) > tol) printf("res.P1[1] wrong: %.5e\n", res.P1[1]);
      if (::fabsf(res.P1[2]) > tol) printf("res.P1[2] wrong: %.5e\n", res.P1[2]);

      printf("(%.5e, %.5e, %.5e), (%.5e, %.5e, %.5e)", res.P0[0], res.P0[1], res.P0[2], res.P1[0], res.P1[1], res.P1[2]);
   }

   __host__ __device__ HorizontalStrip::HorizontalStrip(const double width, const bool fadeBottom) :
      enclosure(nullptr), bnd(nullptr), pos(nullptr), neg(nullptr),
      gt0(nullptr), gt1(nullptr)
   {
      enclosure = new NormalMultiPlaneShapeT();

      const double n0[] = { 0., 0., -1 }, p0[] = { 0., 0., 0 }; // Add top plane
      bnd = new PlaneT(n0, p0);
      enclosure->addPlane(bnd);

      if (fadeBottom) {
         double n1[] = { 0., 0.1, 0.2 }, p1[] = { 0., width / 2., 0. }; // Right end
         Math2::normalize3d(n1, n1);
         pos = new PlaneT(n1, p1);
      }
      else {
         const double n1[] = { 0., 1., 0. }, p1[] = { 0., width / 2., 0. }; // Right end
         pos = new PlaneT(n1, p1);
      }
      enclosure->addPlane(pos);

      const double n2[] = { 0., -1., 0. }, p2[] = { 0., -width / 2., 0. }; // Left end
      double n2norm[3];
      Math2::normalize3d(n2, n2norm);
      neg = new PlaneT(n2norm, p2);
      enclosure->addPlane(neg);
   }

   __host__ __device__ HorizontalStrip::~HorizontalStrip()
   {
      delete enclosure;
      delete bnd;
      delete pos;
      delete neg;

      delete gt0;
      delete gt1;
   }

   __host__ __device__ NormalMultiPlaneShapeT* HorizontalStrip::get()
   {
      return enclosure;
   }

   __host__ __device__ void HorizontalStrip::calcGroundtruth()
   {
      if (!gt0) gt0 = new LineShapeT();
      if (!gt1) gt1 = new LineShapeT();

      const double rightNormal[] = { 1., 0., 0. };
      const double rightPoint[] = { 1.e-5, 0., 0. };
      PlaneT rightBound(rightNormal, rightPoint);
      const double leftNormal[] = { -1., 0., 0. };
      const double leftPoint[] = { -1.e-5, 0., 0. };
      PlaneT leftBound(leftNormal, leftPoint);

      getLineSegment(*bnd, *pos, leftBound, rightBound, *gt0); // right cap
      getLineSegment(*bnd, *neg, leftBound, rightBound, *gt1); // left cap

      //printf("(%.5e, %.5e, %.5e) -> (%.5e, %.5e, %.5e)\n", gt0->P0[0], gt0->P0[1], gt0->P0[2], gt0->P1[0], gt0->P1[1], gt0->P1[2]);
      //printf("(%.5e, %.5e, %.5e) -> (%.5e, %.5e, %.5e)\n", gt1->P0[0], gt1->P0[1], gt1->P0[2], gt1->P1[0], gt1->P1[1], gt1->P1[2]);
   }

   __host__ __device__ void HorizontalStrip::calcRasterization(
      const PlaneT& plane, // light source plane
      const double* axis0, // vector from the plane origin
      const double* axis1, // vector from the plane origin
      const double * lightdir,
      const float xlenperpix,
      const float ylenperpix,
      char* res,
      const unsigned int w,
      const unsigned int h
      ) const
   {
      MultiPlaneShape::LineShape res0, res1; // with respect to origin of plane, along axis0 and axis1
      calcLineProjection(plane, axis0, axis1, lightdir, *gt0, res0); // TODO: 
      calcLineProjection(plane, axis0, axis1, lightdir, *gt1, res1); // TODO: 

      int start0[3] = { res0.P0[0] / xlenperpix, res0.P0[1] / ylenperpix, res0.P0[2] };
      int end0[3] = { res0.P1[0] / xlenperpix, res0.P1[1] / ylenperpix, res0.P1[2] };
      int start1[3] = { res1.P0[0] / xlenperpix, res1.P0[1] / ylenperpix, res1.P0[2] };
      int end1[3] = { res1.P1[0] / xlenperpix, res1.P1[1] / ylenperpix, res1.P1[2] };

      Point s0(start0[0], start0[1]);
      Point s1(start1[0], start1[1]);
      Point e0(end0[0], end0[1]);
      Point e1(end1[0], end1[1]);

      DrawLine(s0, e0, 255, res, w, h);
      DrawLine(s1, e1, 255, res, w, h);

      printf("(%.5e, %.5e, %.5e), (%.5e, %.5e, %.5e) \n", plane.getNormal()[0], plane.getNormal()[1], plane.getNormal()[2], plane.getPoint()[0], plane.getPoint()[1], plane.getPoint()[2]);
      printf("(%.5e, %.5e, %.5e) -> (%.5e, %.5e, %.5e)\n", res0.P0[0], res0.P0[1], res0.P0[2], res0.P1[0], res0.P1[1], res0.P1[2]);
      printf("(%.5e, %.5e, %.5e) -> (%.5e, %.5e, %.5e)\n", res1.P0[0], res1.P0[1], res1.P0[2], res1.P1[0], res1.P1[1], res1.P1[2]);

      printf("(%d, %d, %d) -> (%d, %d, %d)\n", start0[0], start0[1], start0[2], end0[0], end0[1], end0[2]);
      printf("(%d, %d, %d) -> (%d, %d, %d)\n", start1[0], start1[1], start1[2], end1[0], end1[1], end1[2]);
   }

   __host__ __device__ void HorizontalStrip::calcRasterizationCorrection(
      const PlaneT& plane,
      const double* axis0, // vector from the plane origin
      const double* axis1, // vector from the plane origin
      const double* lightdir,
      const float xlenperpix,
      const float ylenperpix,
      char* res,
      const unsigned int w,
      const unsigned int h
      ) const
   {
      MultiPlaneShape::LineShape res1; // with respect to origin of plane, along axis0 and axis1
      calcLineProjection(plane, axis0, axis1, lightdir, *gt1, res1); // TODO: 

      int start1[3] = { res1.P0[0] / xlenperpix, res1.P0[1] / ylenperpix, res1.P0[2] };
      int end1[3] = { res1.P1[0] / xlenperpix, res1.P1[1] / ylenperpix, res1.P1[2] };

      Point s1(start1[0], start1[1]);
      Point e1(end1[0], end1[1]);

      DrawLine(s1, e1, 0, res, w, h);
   }

   __host__ __device__ Washer::Washer(const double innerRadius, const double outerRadius)
   {
      if (innerRadius >= outerRadius) printf("Washer::Washer: innerRadius >= outerRadius (%.5e >= %.5e)\n", innerRadius, outerRadius);

      const double end0Inner[] = { 0., 0., 0. }, end1Inner[] = { 0., 0., -10. };
      inner = new NormalCylindricalShapeT(end0Inner, end1Inner, innerRadius);

      const double end0Outer[] = { 0., 0., 0. }, end1Outer[] = { 0., 0., -10. };
      outer = new NormalCylindricalShapeT(end0Outer, end1Outer, outerRadius);

      diff = new NormalDifferenceShapeT(*outer, *inner);
   }

   __host__ __device__ Washer::~Washer()
   {
      delete diff;
      delete inner;
      delete outer;
   }

   __host__ __device__ void Washer::calcGroundtruth()
   {
   }

   __host__ __device__ void Washer::calcRasterization(const PlaneT&, const double*, const double*, const float, const float, char*, const unsigned int, const unsigned int) const
   {
   }

   __host__ __device__ NormalDifferenceShapeT* Washer::get()
   {
      return diff;
   }
}