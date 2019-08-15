#include "gov\nist\microanalysis\NISTMonte\SimpleBlock.cuh"

namespace SimpleBlock
{
   __host__ __device__ static bool between(double x, double b0, double b1)
   {
      if (!(b0 <= b1)) printf("SimpleBlock::between error: %.10e, %.10e\n", b0, b1);
      return (x >= b0) && (x <= b1);
   }

   SimpleBlock::SimpleBlock(const double corner0[], const double corner1[])
   {
      //assert(corner0.length == 3);
      //assert(corner1.length == 3);
      // Normalize coordinates so that mCorner0[i]<=mCorner1[i]

      memcpy(mCorner0, corner0, sizeof(double) * 3);
      memcpy(mCorner1, corner1, sizeof(double) * 3);

      for (int i = 0; i < 3; ++i)
         if (mCorner0[i] > mCorner1[i]) {
            const double tmp = mCorner0[i];
            mCorner0[i] = mCorner1[i];
            mCorner1[i] = tmp;
         }
   }

   __host__ __device__ bool SimpleBlock::contains(const double pos[]) const
   {
      //assert(pos.length == 3);
      return between(pos[0], mCorner0[0], mCorner1[0]) && between(pos[1], mCorner0[1], mCorner1[1]) && between(pos[2], mCorner0[2], mCorner1[2]);
   }

   __host__ __device__ double SimpleBlock::getFirstIntersection(const double pos0[], const double pos1[])
   {
      //assert(pos0.length == 3);
      //assert(pos1.length == 3);
      double t = INFINITY;
      for (int i = 2; i >= 0; --i) {
         const int j = (i + 1) % 3, k = (i + 2) % 3;
         if (pos1[i] != pos0[i]) {
            double u = (mCorner0[i] - pos0[i]) / (pos1[i] - pos0[i]);
            if ((u >= 0.0) && (u <= t) && between(pos0[j] + u * (pos1[j] - pos0[j]), mCorner0[j], mCorner1[j])
               && between(pos0[k] + u * (pos1[k] - pos0[k]), mCorner0[k], mCorner1[k]))
               t = u;
            // Bottom of block
            u = (mCorner1[i] - pos0[i]) / (pos1[i] - pos0[i]);
            if ((u >= 0.0) && (u <= t) && between(pos0[j] + u * (pos1[j] - pos0[j]), mCorner0[j], mCorner1[j])
               && between(pos0[k] + u * (pos1[k] - pos0[k]), mCorner0[k], mCorner1[k]))
               t = u;
         }
      }
      return t >= 0 ? t : INFINITY;
   }

   //public void render(TrajectoryVRML.RenderContext rc, Writer wr)
   //   throws IOException{
   //   final NumberFormat nf = NumberFormat.getNumberInstance(Locale.US);
   //   nf.setMaximumFractionDigits(3);
   //   nf.setGroupingUsed(false);
   //   final Color color = rc.getCurrentColor();
   //   final String trStr = nf.format(rc.getTransparency());
   //   final String colorStr = nf.format(color.getRed() / 255.0) + " " + nf.format(color.getGreen() / 255.0) + " "
   //      + nf.format(color.getBlue() / 255.0);
   //   wr.append("\nTransform {\n");
   //   wr.append(" translation " + nf.format((mCorner0[0] + mCorner1[0]) / (2.0 * TrajectoryVRML.SCALE)) + " ");
   //   wr.append(nf.format((mCorner0[1] + mCorner1[1]) / (2.0 * TrajectoryVRML.SCALE)) + " ");
   //   wr.append(nf.format((mCorner0[2] + mCorner1[2]) / (2.0 * TrajectoryVRML.SCALE)) + "\n");
   //   wr.append(" children [\n");
   //   wr.append("  Shape {\n");
   //   wr.append("   geometry Box {\n");
   //   wr.append("    size " + nf.format(Math.abs(mCorner1[0] - mCorner0[0]) / TrajectoryVRML.SCALE) + " ");
   //   wr.append(nf.format(Math.abs(mCorner1[1] - mCorner0[1]) / TrajectoryVRML.SCALE) + " ");
   //   wr.append(nf.format(Math.abs(mCorner1[2] - mCorner0[2]) / TrajectoryVRML.SCALE) + "\n");
   //   wr.append("   }\n");
   //   wr.append("   appearance Appearance {\n");
   //   wr.append("    material Material {\n");
   //   wr.append("     emissiveColor " + colorStr + "\n");
   //   wr.append("     transparency " + trStr + "\n");
   //   wr.append("    }\n");
   //   wr.append("   }\n");
   //   wr.append("  }\n");
   //   wr.append(" ]\n");
   //   wr.append("}");
   //}

   const double*  SimpleBlock::getCorner0() const
   {
      return mCorner0;
   }

   const double*  SimpleBlock::getCorner1() const
   {
      return mCorner1;
   }

   
   __host__ __device__ StringT SimpleBlock::toString() const
   {
      return "SimpleBlock(" + amp::to_string(mCorner1[0]) + " " + amp::to_string(mCorner1[1]) + " " + amp::to_string(mCorner1[2]) + "," +
         amp::to_string(mCorner1[0]) + " " + amp::to_string(mCorner1[1]) + " " + amp::to_string(mCorner1[2]) + ")";
   }
}