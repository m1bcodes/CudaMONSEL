#include "gov\nist\nanoscalemetrology\JMONSEL\NormalShapeTransformer.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalIntersectionShape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalCylindricalShape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalMultiPlaneShape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalDifferenceShape.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalUnionShape.cuh"

namespace NormalShapeTransformer
{
   __host__ __device__ void translate(const double distance[], NormalShapeT& shape)
   {
      const StringT& name = shape.toString();
      if (name.starts_with("NormalIntersectionShape")) {
         ((NormalIntersectionShapeT&)shape).translate(distance);
      }
      else if (name.starts_with("NormalCylindricalShape")) {
         ((NormalCylindricalShapeT&)shape).translate(distance);
      }
      else if (name.starts_with("NormalMultiPlaneShape")) {
         ((NormalMultiPlaneShapeT&)shape).translate(distance);
      }
      else if (name.starts_with("NormalDifferenceShape")) {
         ((NormalDifferenceShapeT&)shape).translate(distance);
      }
      else if (name.starts_with("NormalUnionShape")) {
         ((NormalUnionShapeT&)shape).translate(distance);
      }
      else {
         printf("translateNormal: shape not supported: %s\n", shape.toString().c_str());
      }
   }

   __host__ __device__ void rotate(const double pivot[], double phi, double theta, double psi, NormalShapeT& shape)
   {
      const StringT& name = shape.toString();
      if (name.starts_with("NormalIntersectionShape")) {
         ((NormalIntersectionShapeT&)shape).rotate(pivot, phi, theta, psi);
      }
      else if (name.starts_with("NormalCylindricalShape")) {
         ((NormalCylindricalShapeT&)shape).rotate(pivot, phi, theta, psi);
      }
      else if (name.starts_with("NormalMultiPlaneShape")) {
         ((NormalMultiPlaneShapeT&)shape).rotate(pivot, phi, theta, psi);
      }
      else if (name.starts_with("NormalDifferenceShape")) {
         ((NormalDifferenceShapeT&)shape).rotate(pivot, phi, theta, psi);
      }
      else if (name.starts_with("NormalUnionShape")) {
         ((NormalUnionShapeT&)shape).rotate(pivot, phi, theta, psi);
      }
      else {
         printf("rotateNormalShape: shape not supported: %s\n", shape.toString().c_str());
      }
   }
};