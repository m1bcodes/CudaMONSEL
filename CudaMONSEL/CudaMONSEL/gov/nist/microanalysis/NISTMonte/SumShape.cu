#include "gov\nist\microanalysis\NISTMonte\SumShape.cuh"

//#include "gov.nist.microanalysis.EPQLibrary.EPQFatalException.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ITransform.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"
#include "gov\nist\microanalysis\NISTMonte\CylindricalShape.cuh"
#include "gov\nist\microanalysis\NISTMonte\MultiPlaneShape.cuh"
#include "gov\nist\microanalysis\NISTMonte\Sphere.cuh"

#include "gov\nist\nanoscalemetrology\JMONSEL\NormalShapeTransformer.cuh"

namespace SumShape
{
   SumShape::SumShape()
   {
   }

   SumShape::SumShape(ShapeT* const shapes[], int len) : mShapes(shapes, shapes + len)
   {
   }

   __host__ __device__ SumShape::SumShape(ShapeT* const a, ShapeT* const b)
   {
      mShapes.push_back(a);
      mShapes.push_back(b);
   }

   void SumShape::addShape(ShapeT* s)
   {
      mShapes.push_back(s);
   }

   __host__ __device__ bool SumShape::contains(const double pos[]) const
   {
      for (auto &shape : mShapes)
         if (shape->contains(pos))
            return true;
      return false;
   }

   __host__ __device__ double SumShape::getFirstIntersection(const double pos0[], const double pos1[])
   {
      double u = INFINITY;
      if (contains(pos0)) {
         // Starting inside...
         double start[3];
         memcpy(start, pos0, sizeof(double) * 3);
         do {
            double uInc = INFINITY;
            double end[3] = { INFINITY, INFINITY, INFINITY };
            for (auto &shape : mShapes)
               if (shape->contains(start)) {
                  const double ui = shape->getFirstIntersection(start, pos1);
                  if ((ui != INFINITY) && ((uInc == INFINITY) || (ui > uInc))) {
                     double ptbtw0[3];
                     Math2::pointBetween3d(start, pos1, 0.99 * ui, ptbtw0);
                     if (!shape->contains(ptbtw0)) printf("SumShape::getFirstIntersection1 %s%s%s%lf%s%lf%s%lf%s%s%lf\n", shape->toString().c_str(), ", ", "(", pos0[0], ", ", start[1], ", ", pos1[2], ")", ", ", ui);
                     double ptbtw1[3];
                     Math2::pointBetween3d(start, pos1, 1.01 * (ui + 1.0e-3), ptbtw1);
                     if (shape->contains(ptbtw1)) printf("SumShape::getFirstIntersection2: %s%s%s%lf%s%lf%s%lf%s%s%lf\n", shape->toString().c_str(), ", ", "(", pos0[0], ", ", start[1], ", ", pos1[2], ")", ", ", ui);
                     Math2::pointBetween3d(start, pos1, ui, end);
                     uInc = ui;
                     u = Math2::distance3d(end, pos0) / Math2::distance3d(pos1, pos0);
                     if (u > 1.0)
                        break;
                  }
               }
            if (end[0] == INFINITY && end[1] == INFINITY && end[2] == INFINITY)
               break;
            double m0[3];
            Math2::minus3d(pos1, pos0, m0);
            double norm[3];
            Math2::normalize3d(m0, norm);
            double extra[3];
            Math2::multiply3d(1.0e-14, norm, extra);
            // Bump the start point into the next Shape...
            Math2::plus3d(end, extra, start);
            // Repeat until we can take a full step or
            // the step can't be enlarged...
         } while (u < 1.0);
      }
      else
         // Starting outside so get the shortest distance to a boundary
         for (auto &shape : mShapes) {
            const double ui = shape->getFirstIntersection(pos0, pos1);
            if (ui < u)
               u = ui;
         }
      return u;
   }

   __host__ __device__ void rotateShape(const double pivot[], double phi, double theta, double psi, ShapeT* shape)
   {
      StringT& name = shape->toString();
      if (name.starts_with("Normal")) {
         //NormalIntersectionShape::rotateNormalShape(pivot, phi, theta, psi, (NormalShapeT&)(*shape));
         NormalShapeTransformer::rotate(pivot, phi, theta, psi, (NormalShapeT&)(*shape));
      }
      else if (name.starts_with("SumShape")) {
         ((SumShapeT*)shape)->rotate(pivot, phi, theta, psi);
      }
      else if (name.starts_with("CylindricalShape")) {
         ((CylindricalShapeT*)shape)->rotate(pivot, phi, theta, psi);
      }
      else if (name.starts_with("MultiPlaneShape")) {
         ((MultiPlaneShapeT*)shape)->rotate(pivot, phi, theta, psi);
      }
      else if (name.starts_with("Plane")) {
         ((PlaneT*)shape)->rotate(pivot, phi, theta, psi);
      }
      else if (name.starts_with("Sphere")) {
         ((SphereT*)shape)->rotate(pivot, phi, theta, psi);
      }
      else {
         printf("rotateShape: shape not supported %s\n", name.c_str());
      }
   }

   __host__ __device__ void SumShape::rotate(const double pivot[], double phi, double theta, double psi)
   {
      for (auto shape : mShapes) {
         rotateShape(pivot, phi, theta, psi, shape);
      }
   }

   __host__ __device__ void translateShape(const double distance[], ShapeT* shape)
   {
      StringT& name = shape->toString();
      if (name.starts_with("Normal")) {
         NormalShapeTransformer::translate(distance, (NormalShapeT&)(*shape));
      }
      else if (name.starts_with("SumShape")) {
         ((SumShapeT*)shape)->translate(distance);
      }
      else if (name.starts_with("CylindricalShape")) {
         ((CylindricalShapeT*)shape)->translate(distance);
      }
      else if (name.starts_with("MultiPlaneShape")) {
         ((MultiPlaneShapeT*)shape)->translate(distance);
      }
      else if (name.starts_with("Plane")) {
         ((PlaneT*)shape)->translate(distance);
      }
      else if (name.starts_with("Sphere")) {
         ((SphereT*)shape)->translate(distance);
      }
      else {
         printf("translateShape: shape not supported %s\n", name.c_str());
      }
   }

   __host__ __device__ void SumShape::translate(const double distance[])
   {
      for (auto shape : mShapes) {
         translateShape(distance, shape);
      }
   }

   // void render(TrajectoryVRML.RenderContext rc, Writer wr)
   //   throws IOException{
   //   for (const MonteCarloSS.Shape shape : mShapes)
   //   if (shape instanceof TrajectoryVRML.IRender)
   //      ((TrajectoryVRML.IRender) shape).render(rc, wr);
   //}

   __host__ __device__ const amp::vector<ShapeT*>& SumShape::getShapes() const
   {
      return mShapes;
   }

   __host__ __device__ StringT SumShape::toString() const
   {
      StringT res = "SumShape[";
      bool first = true;
      for (auto &shape : mShapes) {
         if (!first)
            res += ", ";
         res += shape->toString();
         first = false;
      }
      res += "]";
      return res;
   }
};
