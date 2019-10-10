#include "gov\nist\microanalysis\NISTMonte\MultiPlaneShape.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\Utility\Transform3D.cuh"

namespace MultiPlaneShape
{
   // construct a Plane object with the specified normal containing the specified point
   __host__ __device__ Plane::Plane(const double normal[], const double point[])
   {
      memcpy(mNormal, normal, sizeof(mNormal[0]) * 3);
      memcpy(mPoint, point, sizeof(mPoint[0]) * 3);
   }

   Plane::Plane(const Plane& other)
   {
      memcpy(mNormal, other.mNormal, sizeof(mNormal[0]) * 3);
      memcpy(mPoint, other.mPoint, sizeof(mPoint[0]) * 3);
   }

   // contains - Is the point p on the inside side of the plane? (Side opposite to the direction of the normal)
   __host__ __device__ bool Plane::contains(const double p[]) const
   {
      //if (!(len == 3)) printf("Plane::contains: len != 3 (%d)", len);
      return ((p[0] - mPoint[0]) * mNormal[0] + (p[1] - mPoint[1]) * mNormal[1] + (p[2] - mPoint[2]) * mNormal[2]) <= 0.0;
   }

   // Is the point close to being contained by this plane?
   bool Plane::almostContains(const double p[])
   {
      //assert(p.length == 3);
      return ((p[0] - mPoint[0]) * mNormal[0] + (p[1] - mPoint[1]) * mNormal[1] + (p[2] - mPoint[2]) * mNormal[2]) <= 0.0;
   }

   // intersection - Where does a line through p1 and p2 intersect the plane?
   // (Double.MAX_VALUE if the line and the plane are parallel or the line intersection occurs before p1)
   __host__ __device__ double Plane::getFirstIntersection(const double p1[], const double p2[])
   {
      //assert(p1.length == 3);
      //assert(p2.length == 3);
      double den = (p2[0] - p1[0]) * mNormal[0] + (p2[1] - p1[1]) * mNormal[1] + (p2[2] - p1[2]) * mNormal[2];
      if (den != 0.0) {
         double res = ((mPoint[0] - p1[0]) * mNormal[0] + (mPoint[1] - p1[1]) * mNormal[1] + (mPoint[2] - p1[2]) * mNormal[2]) / den;
         return res < 0.0 ? INFINITY : res;
      }
      else
         return INFINITY;
   }

   __host__ __device__ void Plane::rotate(const double pivot[], double phi, double theta, double psi)
   {
      Transform3D::rotate3d(mNormal, phi, theta, psi, mNormal);
      Transform3D::rotate3d(mPoint, pivot, phi, theta, psi, mPoint);
   }

   __host__ __device__ void Plane::translate(const double distance[])
   {
      mPoint[0] += distance[0];
      mPoint[1] += distance[1];
      mPoint[2] += distance[2];
   }

   __host__ __device__ const double* Plane::getNormal() const
   {
      return mNormal;
   }

   __host__ __device__ const double* Plane::getPoint() const
   {
      return mPoint;
   }

   static void intersection(const Plane planes[], int len, double res[])
   {
      if (!(len == 3)) printf("MultiPlaneShape::intersection: len == 3 (%d)", len);
      const double* n0 = planes[0].getNormal();
      const double* n1 = planes[1].getNormal();
      const double* n2 = planes[2].getNormal();
      // The determinant of the matrix made from the vector normals
      double det = n0[0] * (n1[1] * n2[2] - n2[1] * n1[2]) // 
         - n1[0] * (n0[1] * n2[2] - n2[1] * n0[2]) // 
         + n2[0] * (n0[1] * n1[2] - n1[1] * n0[2]);
      if (::abs(det) > 1.0e-10) {
         double cross12[3];
         Math2::cross3d(n1, n2, cross12);
         double cross20[3];
         Math2::cross3d(n2, n0, cross20);
         double cross01[3];
         Math2::cross3d(n0, n1, cross01);

         double m0[3];
         Math2::multiply3d(Math2::dot3d(planes[0].getPoint(), n0), cross12, m0);
         double m1[3];
         Math2::multiply3d(Math2::dot3d(planes[1].getPoint(), n1), cross20, m1);
         double m2[3];
         Math2::multiply3d(Math2::dot3d(planes[2].getPoint(), n2), cross01, m2);

         double p0[3];
         Math2::plus3d(m0, m1, p0);
         double p1[3];
         Math2::plus3d(p0, m2, p1);

         Math2::divide3d(p1, det, res);
      }
   }

   static void normalize(const double vec[], double res[])
   {
      double norm = ::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
      res[0] = vec[0] / norm;
      res[1] = vec[1] / norm;
      res[2] = vec[2] / norm;
   }

   static void invert(const double v[], double res[])
   {
      res[0] = -v[0];
      res[1] = -v[1];
      res[2] = -v[2];
   }

   //void MultiPlaneShape::addOffsetPlane(const double normal[], const double pt[], double dist)
   //{
   //   normal = normalize(normal).data();
   //   VectorXd ptVec ({ pt[0] + normal[0] * dist, pt[1] + normal[1] * dist, pt[2] + normal[2] * dist });
   //   addPlane(normal, ptVec.data());
   //}

   __host__ __device__ MultiPlaneShape::MultiPlaneShape()
   {
   }

   __host__ __device__ MultiPlaneShape::MultiPlaneShape(Plane* const planes[], int len) : mPlanes(planes, planes + len)
   {
   }

   //MultiPlaneShape createFilm(const double normal[], const double pt1[], double thickness)
   //{
   //   MultiPlaneShape mp;
   //   mp.addPlane(normal, pt1);
   //   mp.addOffsetPlane(invert(normal).data(), pt1, thickness);
   //   return mp;
   //}

   //MultiPlaneShape createSubstrate(const double normal[], const double pt[])
   //{
   //   MultiPlaneShape mp;
   //   mp.addPlane(normal, pt);
   //   return mp;
   //}

   //MultiPlaneShape createBlock(const double dims[], const double point[], double phi, double theta, double psi)
   //{
   //   double cphi = ::cos(phi), sphi = ::sin(phi);
   //   double cpsi = ::cos(psi), spsi = ::sin(psi);
   //   double cth = ::cos(theta), sth = ::sin(theta);
   //   // Rotated x, y and z-axis normals
   //   MatrixXd normals = {
   //      {
   //         cphi * cth * cpsi - sphi * spsi,
   //         sphi * cpsi + cphi * cth * spsi,
   //         -cphi * sth
   //      },
   //      {
   //         -sphi * cth * cpsi - cphi * spsi,
   //         -sphi * cth * spsi + cphi * cpsi,
   //         sth * sphi
   //      },
   //      {
   //         sth * cpsi,
   //         sth * spsi,
   //         cth
   //      }
   //   };

   //   MultiPlaneShape mp;

   //   for (int i = 0; i < 3; ++i) {
   //      double pt0[] = { point[0] + dims[i] * normals[i][0] / 2.0, point[1] + dims[i] * normals[i][1] / 2.0, point[2] + dims[i] * normals[i][2] / 2.0 };
   //      mp.addPlane(normals[i].data(), pt0);
   //      double pt1[] = { point[0] - dims[i] * normals[i][0] / 2.0, point[1] - dims[i] * normals[i][1] / 2.0, point[2] - dims[i] * normals[i][2] / 2.0 };
   //      mp.addPlane(invert(normals[i].data()).data(), pt1);
   //   }
   //   return mp;
   //}

   //MultiPlaneShape createNamid(double center[], int n, double height, double base)
   //{
   //   if (!(height > 0)) printf("Height must be greater than zero.");
   //   if (!(base > 0)) printf("Base must be greater than zero.");
   //   MultiPlaneShape mp;
   //   mp.addPlane(Math2::Z_AXIS, Math2::ORIGIN_3D);
   //   double theta = -::atan2(base / 2.0, height);
   //   for (int i = 0; i < n; ++i) {
   //      double perp[] = Transform3D.rotate(Math2::X_AXIS, ((double)i / (double)n) * (2.0 * Math2::PI), 0.0, 0.0);
   //      double nn[] = Transform3D.rotate(Math2::X_AXIS, 0.0, -theta, ((double)i / (double)n) * (2.0 * Math2::PI));
   //      mp.addPlane(nn, Math2::multiply(base / 2.0, VectorXd(perp)).data());
   //   }
   //   mp.translate(Math2::multiply(-height, VectorXd(Math2::Z_AXIS, Math2::Z_AXIS+3)).data());
   //   mp.translate(center);
   //   return mp;
   //}

   //void MultiPlaneShape::addPlane(const double normal[], const double point[])
   //{
   //   Plane newPlane(normalize(normal).data(), 3, point, 3); // TODO: will not work
   //   mPlanes.push_back(&newPlane);
   //}

   __host__ __device__ void MultiPlaneShape::addPlane(Plane *plane)
   {
      mPlanes.push_back(plane);
   }

   __host__ __device__ bool MultiPlaneShape::contains(const double pos[]) const
   {
      if (mPlanes.size() == 1)
         return (mPlanes.at(0))->contains(pos);
      else {
         for (auto &pl : mPlanes)
            if (!pl->contains(pos))
               return false;
         return true;
      }
   }

   __host__ __device__ static bool IsSamePosition(const double a[], const double b[])
   {
      return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
   }

   __host__ __device__ double MultiPlaneShape::getFirstIntersection(const double pos0[], const double pos1[])
   {
      if (mPlanes.size() == 1)
         return (mPlanes.at(0))->getFirstIntersection(pos0, pos1);
      else {
         double minU = INFINITY;
         if (IsSamePosition(pos0, mInsidePos) || contains(pos0)) { // Easy part...
            // If we start inside then any plane we strike will take us outside
            for (auto &pl : mPlanes) {
               const double u = pl->getFirstIntersection(pos0, pos1);
               // assert (u >= 0.0);
               if (u < minU)
                  minU = u;
            }
            // Can optimize 'contains(pos0)' by storing pos1 when u>1.0. Since
            // in all likelyhood during the next
            // call to MCSS_MultiPlaneShape.getFirstIntersection, pos0 will
            // equal the stored pos1.
            if ((minU > 1.0) && (minU != INFINITY)) {
               mInsidePos[0] = pos1[0];
               mInsidePos[1] = pos1[1];
               mInsidePos[2] = pos1[2];
            }
         }
         else { // A little more difficult...
            // If we start outside then we may intersect a plane and yet remain
            // outside
            double testPt[3];
            for (auto &pli : mPlanes) {
               const double u = pli->getFirstIntersection(pos0, pos1);
               if (!(u >= 0.0)) printf("MultiPlaneShape::getFirstIntersection: u < 0.0 (%.10e)", u);
               if (u < minU) {
                  // Don't bother to test the intersection unless it
                  // is potentially closer...
                  testPt[0] = pos0[0] + u * (pos1[0] - pos0[0]);
                  testPt[1] = pos0[1] + u * (pos1[1] - pos0[1]);
                  testPt[2] = pos0[2] + u * (pos1[2] - pos0[2]);
                  bool inside = true;
                  for (auto &plj : mPlanes)
                     if (plj != pli)
                        if (!plj->contains(testPt)) {
                           inside = false;
                           break;
                        }
                  if (inside)
                     minU = u;
               }
            }
            if (minU < 1.0) {
               mInsidePos[0] = pos1[0];
               mInsidePos[1] = pos1[1];
               mInsidePos[2] = pos1[2];
            }
         }
         return minU;
      }
   }

   __host__ __device__ void MultiPlaneShape::rotate(const double pivot[], double phi, double theta, double psi)
   {
      for (auto t : mPlanes)
         t->rotate(pivot, phi, theta, psi);
   }

   __host__ __device__ void MultiPlaneShape::translate(const double distance[])
   {
      for (auto t : mPlanes)
         t->translate(distance);
   }

   //public void render(TrajectoryVRML.RenderContext rc, Writer wr)
   //   throws IOException{
   //   // Configure the number format
   //   final NumberFormat nf = NumberFormat.getInstance(Locale.US);
   //   nf.setMaximumFractionDigits(3);
   //   nf.setGroupingUsed(false);
   //   // Configure the rendering color
   //   final Color color = rc.getCurrentColor();
   //   final String trStr = nf.format(rc.getTransparency());
   //   final String colorStr = nf.format(color.getRed() / 255.0) + " " + nf.format(color.getGreen() / 255.0) + " "
   //      + nf.format(color.getBlue() / 255.0);
   //   // Points are the intersection of three planes
   //   final Plane[] ps = new Plane[3];
   //   // For each plane find the points located on it...
   //   for (int i = 0; i < mPlanes.size(); ++i) {
   //      ps[0] = mPlanes.get(i);
   //      // An array of plane indexes int[2]
   //      final ArrayList<int[]> idxList = new ArrayList<int[]>();
   //      final HashMap<int[], double[]> ptMap = new HashMap<int[], double[]>();
   //      // For each plane find the other 2*n planes that define the
   //      for (int j = 0; j < mPlanes.size(); ++j)
   //         if (i != j) {
   //            ps[1] = mPlanes.get(j);
   //            for (int k = j + 1; k < mPlanes.size(); ++k)
   //               if (i != k) {
   //                  ps[2] = mPlanes.get(k);
   //                  final double[] pt = intersection(ps);
   //                  if (pt != null) {
   //                     // Check to ensure that this point is inside
   //                     // all the other planes...
   //                     boolean inside = true;
   //                     for (int m = 0; inside && (m < mPlanes.size()); m++)
   //                        if ((m != i) && (m != j) && (m != k)) {
   //                           final Plane p = mPlanes.get(m);
   //                           inside = p.almostContains(pt);
   //                        }
   //                     if (inside) {
   //                        boolean add = true;
   //                        for (final Iterator<double[]> m = ptMap.values().iterator(); add && m.hasNext();) {
   //                           final double[] inclPt = m.next();
   //                           final double d = Math2.distance(inclPt, pt);
   //                           add = (d > 1.0e-12);
   //                        }
   //                        if (add) {
   //                           final int[] idxItem = new int[] {
   //                              j,
   //                                 k
   //                           };
   //                           // This point forms one of the corners for this
   //                           // plane
   //                           idxList.add(idxItem);
   //                           ptMap.put(idxItem, pt);
   //                        }
   //                     }
   //                  }
   //               }
   //         }
   //      if (ptMap.size() > 2) {
   //         // Put idxList in order such that we can step through them by
   //         // permuting only one index per step.
   //         for (int j = 0; j < idxList.size(); ++j) {
   //            final int[] iJ = idxList.get(j);
   //            for (int k = j + 1; k < idxList.size(); ++k) {
   //               final int[] iK = idxList.get(k);
   //               // Is the first or second index duplicated?
   //               final boolean b0 = ((iJ[0] == iK[0]) || (iJ[0] == iK[1]));
   //               final boolean b1 = ((iJ[1] == iK[0]) || (iJ[1] == iK[1]));
   //               // Only permute one index
   //               if (b0 ^ b1) {
   //                  // Swap iK into j+1
   //                  final int[] iJp1 = idxList.get(j + 1);
   //                  idxList.set(j + 1, iK);
   //                  idxList.set(k, iJp1);
   //                  break;
   //               }
   //            }
   //         }
   //         wr.append("\n#Rendering side " + Integer.toString(i + 1) + " of " + Integer.toString(mPlanes.size()));
   //         // Write each face separately (not solid)
   //         wr.append("\nShape {\n");
   //         wr.append(" geometry IndexedFaceSet {\n");
   //         wr.append("  coord Coordinate {\n");
   //         wr.append("   point [ ");
   //         for (int j = 0; j < idxList.size(); ++j) {
   //            if (j != 0)
   //               wr.append(", ");
   //            final double[] pt = ptMap.get(idxList.get(j));
   //            wr.append(nf.format(pt[0] / TrajectoryVRML.SCALE));
   //            wr.append(" ");
   //            wr.append(nf.format(pt[1] / TrajectoryVRML.SCALE));
   //            wr.append(" ");
   //            wr.append(nf.format(pt[2] / TrajectoryVRML.SCALE));
   //         }
   //         wr.append(" ]\n");
   //         wr.append("  }\n");
   //         wr.append("  coordIndex [ ");
   //         for (int j = 0; j < idxList.size(); ++j) {
   //            if (j != 0)
   //               wr.append(" ");
   //            wr.append(Integer.toString(j));
   //         }
   //         wr.append(" ]\n");
   //         wr.append("  solid FALSE\n");
   //         wr.append(" }\n");
   //         wr.append(" appearance Appearance {\n");
   //         wr.append("  material Material {\n");
   //         wr.append("   emissiveColor " + colorStr + "\n");
   //         wr.append("   transparency " + trStr + "\n");
   //         wr.append("  }\n");
   //         wr.append(" }\n");
   //         wr.append("}");
   //      }
   //   }
   //   wr.flush();
   //}

   MultiPlaneShape::Planes MultiPlaneShape::getPlanes() const
   {
      return mPlanes;
   }

   __host__ __device__ StringT Plane::toString() const
   {
      return "Plane";
   }

   __host__ __device__ StringT MultiPlaneShape::toString() const
   {
      return "MultiPlaneShape";
   }

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ > 0))
   __constant__ const double SMALL_NUM = 1.e-25f; // anything that avoids division overflow
#else
   const double SMALL_NUM = 1.e-25f; // anything that avoids division overflow
#endif

   __host__ __device__ LineShape::LineShape()
   {
      memset(P0, 0, sizeof(P0[0]) * 3);
      memset(P1, 0, sizeof(P1[0]) * 3);
   }

   __host__ __device__ LineShape::LineShape(const double* p0, const double* p1)
   {
      memcpy(P0, p0, sizeof(P0[0]) * 3);
      memcpy(P1, p1, sizeof(P1[0]) * 3);
   }

   // intersect3D_SegmentPlane(): find the 3D intersection of a segment and a plane
   //    Input:  S = a segment, and Pn = a plane = {Point V0;  Vector n;}
   //    Output: *I0 = the intersect point (when it exists)
   //    Return: 0 = disjoint (no intersection)
   //            1 = intersection in the unique point *I0
   //            2 = the segment lies in the plane
   __host__ __device__ int intersect3D_SegmentPlane(const LineShape& S, const Plane& Pn, double* I)
   {
      double u[3];
      Math2::minus3d(S.P1, S.P0, u);
      double w[3];
      Math2::minus3d(S.P0, Pn.getPoint(), w);

      const double D = Math2::dot3d(Pn.getNormal(), u);
      const double N = -Math2::dot3d(Pn.getNormal(), w);

      if (fabs(D) < SMALL_NUM) { // segment is parallel to plane
         if (N == 0) // segment lies in plane
            return 2;
         else
            return 0; // no intersection
      }
      // they are not parallel
      // compute intersect param
      const double sI = N / D;
      if (sI < 0 || sI > 1)
         return 0; // no intersection
      Math2::multiply3d(sI, u, I);
      Math2::plus3d(S.P0, I, I); // compute segment intersect point

      return 1;
   }

   // http://geomalgorithms.com/a05-_intersect-1.html
   // intersect3D_2Planes(): find the 3D intersection of two planes
   //    Input: two planes Pn1 and Pn2
   //    Output: L = the intersection line (when it exists)
   //    Return: 0 = disjoint (no intersection)
   //            1 = the two planes coincide
   //            2 = intersection in the unique line L
   __host__ __device__ int intersect3D_2Planes(const Plane& Pn1, const Plane& Pn2, LineShape& L)
   {
      double u[3];
      Math2::cross3d(Pn1.getNormal(), Pn2.getNormal(), u); // cross product

      double ax = (u[0] >= 0 ? u[0] : -u[0]);
      double ay = (u[1] >= 0 ? u[1] : -u[1]);
      double az = (u[2] >= 0 ? u[2] : -u[2]);

      // test if the two planes are parallel
      if ((ax + ay + az) < SMALL_NUM) { // Pn1 and Pn2 are near parallel
         // test if disjoint or coincide
         double v[3];
         Math2::minus3d(Pn2.getPoint(), Pn1.getPoint(), v);
         
         if (Math2::dot3d(Pn1.getNormal(), v) == 0) // Pn2.V0 lies in Pn1
            return 1; // Pn1 and Pn2 coincide
         else
            return 0; // Pn1 and Pn2 are disjoint
      }

      // Pn1 and Pn2 intersect in a line
      // first determine max abs coordinate of cross product
      int maxc; // max coordinate
      if (ax > ay) {
         if (ax > az)
            maxc = 1;
         else
            maxc = 3;
      }
      else {
         if (ay > az)
            maxc = 2;
         else maxc = 3;
      }

      // next, to get a point on the intersect line
      // zero the max coord, and solve for the other two
      double iP[3]; // intersect point
      double d1, d2; // the constants in the 2 plane equations
      d1 = -Math2::dot3d(Pn1.getNormal(), Pn1.getPoint()); // note: could be pre-stored  with plane
      d2 = -Math2::dot3d(Pn2.getNormal(), Pn2.getPoint()); // ditto

      switch (maxc)
      { // select max coordinate
      case 1: // intersect with x=0
         iP[0] = 0;
         iP[1] = (d2*Pn1.getNormal()[2] - d1*Pn2.getNormal()[2]) / u[0];
         iP[2] = (d1*Pn2.getNormal()[1] - d2*Pn1.getNormal()[1]) / u[0];
         break;
      case 2: // intersect with y=0
         iP[0] = (d1*Pn2.getNormal()[2] - d2*Pn1.getNormal()[2]) / u[1];
         iP[1] = 0;
         iP[2] = (d2*Pn1.getNormal()[0] - d1*Pn2.getNormal()[0]) / u[1];
         break;
      case 3: // intersect with z=0
         iP[0] = (d2*Pn1.getNormal()[1] - d1*Pn2.getNormal()[1]) / u[2];
         iP[1] = (d1*Pn2.getNormal()[0] - d2*Pn1.getNormal()[0]) / u[2];
         iP[2] = 0;
      }
      memcpy(L.P0, iP, sizeof(L.P0[0]) * 3);
      Math2::plus3d(iP, u, iP);
      memcpy(L.P1, iP, sizeof(L.P1[0]) * 3);

      return 2;
   }

   void intersect3D_2PlanesTest()
   {
      const double n0[3] = { 0.f, 0.f, 1.f };
      const double s0[3] = { 0.f, 0.f, 1.f };

      const double n1[3] = { 0.f, 1.f, 0.f };
      const double s1[3] = { 0.f, 1.f, 0.f };

      PlaneT pl0(n0, s0);
      PlaneT pl1(n1, s1);
      LineShape l;
      printf("%d\n", intersect3D_2Planes(pl0, pl1, l));
      printf("%lf, %lf, %lf\n", l.P0[0], l.P0[1], l.P0[2]);
      printf("%lf, %lf, %lf\n", l.P1[0], l.P1[1], l.P1[2]);
   }
}