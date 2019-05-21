#include "gov\nist\microanalysis\NISTMonte\SumShape.cuh"

//#include "gov.nist.microanalysis.EPQLibrary.EPQFatalException.cuh"
#include "gov\nist\microanalysis\EPQLibrary\ITransform.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"

namespace SumShape
{
   /**
   * Creates a sum shape that represents the sum of an array of Shapes.
   *
   * @param shapes Shape[]
   */
   SumShape::SumShape(ShapeT* const shapes[], int len) : mShapes(shapes, shapes + len)
   {
   }

   /**
   * Create a sum shape that represents the sum of two shapes.
   *
   * @param a Shape
   * @param b Shape
   */
   SumShape::SumShape(ShapeT* const a, ShapeT* const b)
   {
      mShapes.push_back(a);
      mShapes.push_back(b);
   }

   bool SumShape::isNormalShape() const
   {
      return false;
   }

   /**
   * @see gov.nist.microanalysis.NISTMonte.MonteCarloSS.Shape#contains(double[])
   */
   bool SumShape::contains(const double pos[]) const
   {
      for (auto shape : mShapes)
         if (shape->contains(pos))
            return true;
      return false;
   }

   /**
   * @see gov.nist.microanalysis.NISTMonte.MonteCarloSS.Shape#getFirstIntersection(double[],
   *      double[])
   */
   double SumShape::getFirstIntersection(const double pos0[], const double pos1[])
   {
      double u = INFINITY;
      if (contains(pos0)) {
         // Starting inside...
         VectorXd start(pos0, pos0 + 3);
         do {
            double uInc = INFINITY;
            VectorXd end;
            for (auto shape : mShapes)
               if (shape->contains(start.data())) {
                  const double ui = shape->getFirstIntersection(start.data(), pos1);
                  if ((ui != INFINITY) && ((uInc == INFINITY) || (ui > uInc))) {
                     if (!shape->contains(Math2::pointBetween(start, VectorXd(pos1, pos1 + 3), 0.99 * ui).data())) printf("%s", (shape->toString() + ", "
                        + "(" + std::to_string(pos0[0]) + ", " + std::to_string(pos0[1]) + ", " + std::to_string(pos0[2]) + ")" + ", "
                        + "(" + std::to_string(start[0]) + ", " + std::to_string(start[1]) + ", " + std::to_string(start[2]) + ")" + ", "
                        + "(" + std::to_string(pos1[0]) + ", " + std::to_string(pos1[1]) + ", " + std::to_string(pos1[2]) + ")" + ", "
                        + std::to_string(ui)).c_str());
                     if (shape->contains(Math2::pointBetween(start, VectorXd(pos1, pos1 + 3), 1.01 * (ui + 1.0e-3)).data())) printf("%s\n", (shape->toString() + ", "
                        + "(" + std::to_string(pos0[0]) + ", " + std::to_string(pos0[1]) + ", " + std::to_string(pos0[2]) + ")" + ", "
                        + "(" + std::to_string(start[0]) + ", " + std::to_string(start[1]) + ", " + std::to_string(start[2]) + ")" + ", "
                        + "(" + std::to_string(pos1[0]) + ", " + std::to_string(pos1[1]) + ", " + std::to_string(pos1[2]) + ")" + ", "
                        + std::to_string(ui)).c_str());
                     end = Math2::pointBetween(start, VectorXd(pos1, pos1+3), ui);
                     uInc = ui;
                     u = Math2::distance(end, VectorXd(pos0, pos0 + 3)) / Math2::distance(VectorXd(pos1, pos1 + 3), VectorXd(pos0, pos0+3));
                     if (u > 1.0)
                        break;
                  }
               }
            if (end.empty())
               break;
            const VectorXd extra = Math2::multiply(1.0e-14, Math2::normalize(Math2::minus(VectorXd(pos1, pos1 + 3), VectorXd(pos0, pos0 + 3))));
            // Bump the start point into the next Shape...
            start = Math2::plus(end, extra);
            // Repeat until we can take a full step or
            // the step can't be enlarged...
         } while (u < 1.0);
      }
      else
         // Starting outside so get the shortest distance to a boundary
         for (auto shape : mShapes) {
            const double ui = shape->getFirstIntersection(pos0, pos1);
            if (ui < u)
               u = ui;
         }
      return u;
   }

   /**
   * @see gov.nist.microanalysis.EPQLibrary.ITransform#rotate(double[], double,
   *      double, double)
   */
   void SumShape::rotate(const double pivot[], double phi, double theta, double psi)
   {
      for (auto shape : mShapes) {
         //if (!(shape instanceof ITransform))
         //   throw new EPQFatalException(shape.toString() + " does not support transformation.");
         ((ITransformT*)shape)->rotate(pivot, phi, theta, psi);
      }
   }

   /**
   * @see gov.nist.microanalysis.EPQLibrary.ITransform#translate(double[])
   */
   void SumShape::translate(const double distance[])
   {
      for (auto shape : mShapes) {
         //if (!(shape instanceof ITransform))
         //   throw new EPQFatalException(shape.toString() + " does not support transformation.");
         ((ITransformT*)shape)->translate(distance);
      }
   }

   /**
   * Render the SumShape by rendering each of the sub-Shapes. If a sub-Shape
   * does not implement the interface TrajectoryVRML.IRender then it will be
   * missing from the rendered VRML world.
   *
   * @param rc
   * @param wr
   * @throws IOException
   * @see gov.nist.microanalysis.NISTMonte.TrajectoryVRML.IRender#render(gov.nist.microanalysis.NISTMonte.TrajectoryVRML.RenderContext,
   *      java.io.Writer)
   */
   // void render(TrajectoryVRML.RenderContext rc, Writer wr)
   //   throws IOException{
   //   for (const MonteCarloSS.Shape shape : mShapes)
   //   if (shape instanceof TrajectoryVRML.IRender)
   //      ((TrajectoryVRML.IRender) shape).render(rc, wr);
   //}

   /**
   * Returns an immutable list of the Shapes which define this SumShape object.
   *
   * @return Returns the shapes.
   */
   std::vector<ShapeT*> SumShape::getShapes() const
   {
      return mShapes;
   }

   StringT SumShape::toString() const
   {
      StringT res = "Sum[";
      bool first = true;
      for (auto shape : mShapes) {
         if (!first)
            res.append(", ");
         res.append(shape->toString());
         first = false;
      }
      res.append("]");
      return res;
   }
};
