#include "gov\nist\microanalysis\EPQLibrary\PhysicalConstants.cuh"
#include "gov\nist\microanalysis\NISTMonte\RegionBase.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
#include "gov\nist\microanalysis\NISTMonte\Shape.cuh"
#include "gov\nist\microanalysis\NISTMonte\MonteCarloSS.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\ExpQMBarrierSM.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\SEMaterial.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NormalShape.cuh"

namespace ExpQMBarrierSM
{
   ExpQMBarrierSM::ExpQMBarrierSM(const MaterialT* mat) :
      u0(mat->isSEmaterial() ? ((SEmaterialT*)mat)->getEnergyCBbottom() : 0),
      classical(true)
   {
   }

   ExpQMBarrierSM::ExpQMBarrierSM(const MaterialT* mat, double lambda) : 
      u0(mat->isSEmaterial() ? ((SEmaterialT*)mat)->getEnergyCBbottom() : 0),
      classical(false),
      lambda(lambda),
      lambdaFactor((PhysicalConstants::PI * lambda * ::sqrt(PhysicalConstants::ElectronMass / 2.)) / PhysicalConstants::PlanckReduced)
   {
   }

   static double sharpBarrierT(double rootPerpE, double rootDiff)
   {
      double ksum = rootPerpE + rootDiff;
      return (4. * rootPerpE * rootDiff) / (ksum * ksum);
   }

   double ExpQMBarrierSM::generalBarrierT(double rootPerpE, double rootDiff) const
   {
      double k1 = lambdaFactor * rootPerpE;
      if (k1 > 50.)
         return 1.;
      double k2 = lambdaFactor * rootDiff;
      double kplus = k1 + k2;
      double kminus = k1 - k2;
      double sinhPlus = ::sinh(kplus);
      double sinhMinus = ::sinh(kminus);
      double ratio = sinhMinus / sinhPlus;
      return 1. - (ratio * ratio);
   }

   ElectronT* ExpQMBarrierSM::barrierScatter(ElectronT* pe, const RegionBaseT* nextRegion) const
   {
      const MaterialT& nextmaterial = nextRegion->getMaterial();
      const RegionBaseT* currentRegion = pe->getCurrentRegion();

      if (!(currentRegion != nextRegion)) {
         printf("ExpQMBarrierSM::barrierScatter: currentRegion == nextRegion");
         return NULL;
      }

      double deltaU;
      const MaterialT& currentMaterial = currentRegion->getMaterial();
      if (currentMaterial.isSEmaterial())
         deltaU = -((SEmaterialT)currentMaterial).getEnergyCBbottom();
      else
         deltaU = 0.;

      if (!(deltaU == -u0)) {
         printf("ExpQMBarrierSM::barrierScatter: deltaU != -u0 (%.10e, %.10e)", deltaU, -u0);
         return NULL;
      }
      if (nextmaterial.isSEmaterial())
         deltaU += ((SEmaterialT)nextmaterial).getEnergyCBbottom();

      /* FIND THE OUTWARD POINTING NORMAL AT THE BOUNDARY */
      PositionVecT nb; // We'll store it here

      if (currentRegion->isContainingRegion(*nextRegion)) {
         const RegionBaseT* struckRegion = nextRegion; // usually this is true
         /*
         * Sometimes we cross multiple boundaries at once. The while loop
         * checks and corrects for this.
         */
         while (struckRegion->getParent() != currentRegion)
            struckRegion = struckRegion->getParent();
         const ShapeT* intersectedshape = struckRegion->getShape();
         if (intersectedshape->isNormalShape()) {
            nb = ((NormalShapeT*)intersectedshape)->getPreviousNormal();
            for (int i = 0; i < 3; i++)
               nb[i] *= -1;
         }
      }
      else {
         const ShapeT* intersectedshape = currentRegion->getShape();
         if (intersectedshape->isNormalShape())
            nb = ((NormalShapeT*)intersectedshape)->getPreviousNormal();
      }

      // GET THE VECTOR IN THE ELECTRON'S DIRECTION OF MOTION
      double theta0 = pe->getTheta();
      double phi0 = pe->getPhi();
      double sintheta0 = ::sin(theta0);
      PositionVecT n0({ sintheta0 * ::cos(phi0), sintheta0 * ::sin(phi0), ::cos(theta0) });

      /*
      * If the intersected shape is not a NormalShape, we still haven't
      * initialized nb. We have no data in this case, so we must make do with
      * an arbitrary assignment: Let nb be the same as the electron direction.
      * This choice gives maximum transmission probability and no deflection of
      * the electron's path.
      */
      if (!nb.empty())
         nb = n0;

      /*
      * Let the angle of incidence be called alpha. Cos(alpha) is given by the
      * dot product
      */
      double cosalpha = (n0[0] * nb[0]) + (n0[1] * nb[1]) + (n0[2] * nb[2]);

      if (cosalpha <= 0.) {
         /*
         * This case corresponds to the electron "hitting" the barrier while
         * moving away from it. I.e., it didn't really hit the barrier. This
         * can happen, e.g., if electric field alters the electron's direction
         * of motion. We give it a nudge away from the barrier towards the
         * inside
         */
         PositionVecT pos0 = pe->getPosition();
         double tmppos[] = { pos0[0] - (MonteCarloSS::SMALL_DISP * nb[0]), pos0[1] - (MonteCarloSS::SMALL_DISP * nb[1]), pos0[2] - (MonteCarloSS::SMALL_DISP * nb[2]) };
         pe->setPosition(tmppos);

         return NULL;
      }

      if (deltaU == 0.) {
         /*
         * This corresponds to no barrier. This is usually due to a
         * mathematical boundary with the same material on both sides. It
         * transmits, so we give it a nudge off of the barrier toward the
         * outside, update the electron's region, and return.
         */
         auto pos0 = pe->getPosition();
         double tmppos[] = { pos0[0] + (MonteCarloSS::SMALL_DISP * nb[0]), pos0[1] + (MonteCarloSS::SMALL_DISP * nb[1]), pos0[2] + (MonteCarloSS::SMALL_DISP * nb[2]) };
         pe->setPosition(tmppos);
         pe->setCurrentRegion(nextRegion);

         return NULL;
      }

      double kE0 = pe->getEnergy();
      double perpE;
      if (kE0 <= 0.)
         perpE = 0.;
      else
         perpE = cosalpha * cosalpha * kE0;
      double rootPerpE = 0.;
      double rootDiff = 0.;

      /* DECIDE WHETHER IT TRANSMITS OR NOT */
      bool transmits;
      if ((perpE == 0.) || (perpE <= deltaU))
         /*
         * Even if deltaU<0 (the electron is stepping downhill) the quantum
         * mechanical formula gives transmission = 0 when perpE = 0.
         */
         transmits = false;
      else {
         rootPerpE = ::sqrt(perpE);
         rootDiff = ::sqrt(perpE - deltaU);
         if (classical)
            transmits = true; // Since we already know perpE>deltaU
         else {
            double transmissionProb;
            if (lambda == 0.)
               transmissionProb = sharpBarrierT(rootPerpE, rootDiff);
            else
               transmissionProb = generalBarrierT(rootPerpE, rootDiff);

            double r = Math2::random();
            transmits = r < transmissionProb;
         }
      }

      /*
      * COMPUTE DIRECTION AND ENERGY FOR EACH OF THE CASES: TRANSMISSION
      * OR REFLECTION
      */
      PositionVecT nf(3, 0); // Direction vector after scattering
      if (transmits) { // Transmission
         double factor = cosalpha * ((rootDiff / rootPerpE) - 1.);
         for (int i = 0; i < 3; i++)
            nf[i] = n0[i] + (factor * nb[i]);
         /* Normalize the z component to use later computing theta. */
         // nf[2] /= Math.sqrt(1. + (2. * cosalpha + factor) * factor);
         nf[2] /= ::sqrt((nf[0] * nf[0]) + (nf[1] * nf[1]) + (nf[2] * nf[2]));

         pe->setEnergy(kE0 - deltaU);
         pe->setCurrentRegion(nextRegion);
         auto pos0 = pe->getPosition();

         double tmppos[] = { pos0[0] + (MonteCarloSS::SMALL_DISP * nb[0]), pos0[1] + (MonteCarloSS::SMALL_DISP * nb[1]), pos0[2] + (MonteCarloSS::SMALL_DISP * nb[2]) };
         pe->setPosition(tmppos);
      }
      else { // Total internal reflection
         double twocosalpha = 2. * cosalpha;
         for (int i = 0; i < 3; i++)
            nf[i] = n0[i] - (nb[i] * twocosalpha);

         auto pos0 = pe->getPosition();
         double tmppos[] = { pos0[0] - (MonteCarloSS::SMALL_DISP * nb[0]), pos0[1] - (MonteCarloSS::SMALL_DISP * nb[1]), pos0[2] - (MonteCarloSS::SMALL_DISP * nb[2]) };
         pe->setPosition(tmppos);
      }

      double thetaf = ::acos(nf[2]);
      double phif = ::atan2(nf[1], nf[0]);

      pe->setDirection(thetaf, phif);
      return NULL;
   }

   StringT ExpQMBarrierSM::toString() const
   {
      return "ExpQMBarrierSM(" + mat->toString() + "," + std::to_string(u0) + "," + std::to_string(lambda) + ")";
   }
}