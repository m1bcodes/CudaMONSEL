#include "gov\nist\nanoscalemetrology\JMONSEL\TabulatedInelasticSM.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NUTableInterpolation.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\SEmaterial.cuh"
#include "gov\nist\microanalysis\NISTMonte\Electron.cuh"
#include "gov\nist\microanalysis\Utility\Math2.cuh"

#include "Amphibian\Algorithm.cuh"

namespace TabulatedInelasticSM
{
   TabulatedInelasticSM::TabulatedInelasticSM(const SEmaterialT& mat, int methodSE, const char* tables[]) :
      methodSE((methodSE != 2) && (methodSE != 3) ? 1 : methodSE),
      energyOffset(0.),
      minEgenSE(0.),
      tableIIMFP(NUTableInterpolation::getInstance(tables[0])),
      tableReducedDeltaE(NUTableInterpolation::getInstance(tables[1])),
      tableTheta(NUTableInterpolation::getInstance(tables[2])),
      tableSEE0((methodSE == 2) || (methodSE == 3) ? NUTableInterpolation::getInstance(tables[3]) : NULL),
      rateMult(1.),
      E0fromDispersion(false),
      kEa(VectorXd(1)),
      interpInput(VectorXd(3)),
      defaultRatios(true)
   {
      setMaterial(&mat);
   }

   TabulatedInelasticSM::TabulatedInelasticSM(const SEmaterialT& mat, int methodSE, const char* tables[], double energyOffset) :
      methodSE((methodSE != 2) && (methodSE != 3) ? 1 : methodSE),
      energyOffset(energyOffset),
      minEgenSE(0.),
      tableIIMFP(NUTableInterpolation::getInstance(tables[0])),
      tableReducedDeltaE(NUTableInterpolation::getInstance(tables[1])),
      tableTheta(NUTableInterpolation::getInstance(tables[2])),
      tableSEE0((methodSE == 2) || (methodSE == 3) ? NUTableInterpolation::getInstance(tables[3]) : NULL),
      energyRangeSE0((methodSE == 2) || (methodSE == 3) ? tableSEE0->getRange() : VectorXd()),
      rateMult(1.), 
      E0fromDispersion(false),
      kEa(VectorXd(1)), 
      interpInput(VectorXd(3)),
      defaultRatios(true)
   {
      ///* Read interpolation tables into memory */
      //try {
      //   tableIIMFP = NUTableInterpolation.getInstance(tables[0]);
      //}
      //catch (final FileNotFoundException e1) {
      //   throw new EPQFatalException("File " + tables[0] + " not found.");
      //}
      //try {
      //   tableReducedDeltaE = NUTableInterpolation.getInstance(tables[1]);
      //}
      //catch (final FileNotFoundException e1) {
      //   throw new EPQFatalException("File " + tables[1] + " not found.");
      //}
      //try {
      //   tableTheta = NUTableInterpolation.getInstance(tables[2]);
      //}
      //catch (final FileNotFoundException e1) {
      //   throw new EPQFatalException("File " + tables[2] + " not found.");
      //}
      //if ((methodSE == 2) || (methodSE == 3))
      //   try {
      //   tableSEE0 = NUTableInterpolation.getInstance(tables[3]);
      //   energyRangeSE0 = tableSEE0.getRange();
      //}
      //catch (final FileNotFoundException e1) {
      //   throw new EPQFatalException("File " + tables[3] + " not found.");
      //}
      setMaterial(&mat);
   }

   VectorXd updateDirection(double theta, double phi, double dTheta, double dPhi)
   {
      double ct = ::cos(theta), st = ::sin(theta);
      double cp = ::cos(phi), sp = ::sin(phi);
      double ca = ::cos(dTheta), sa = ::sin(dTheta);
      double cb = ::cos(dPhi);

      double xx = (cb * ct * sa) + (ca * st);
      double yy = sa * ::sin(dPhi);
      double dx = (cp * xx) - (sp * yy);
      double dy = (cp * yy) + (sp * xx);
      double dz = (ca * ct) - (cb * sa * st);

      theta = ::atan2(::sqrt((dx * dx) + (dy * dy)), dz);
      phi = ::atan2(dy, dx);
      VectorXd res(2);
      res[0] = theta;
      res[1] = phi;
      return res;
   }

   VectorXd TabulatedInelasticSM::simESEf(double Eq, double deltaE, double r)
   {
      double q = ::sqrt(Eq);
      double kz = (deltaE - Eq) / 2. / q;
      double kzf = kz + q;
      double Ezq = kzf * kzf;
      double minE = (offsetFermiEnergy + bandgap) - Ezq;
      if (minE < 0.) minE = 0.;
      double maxE = offsetFermiEnergy - (kz * kz);
      if (!(minE <= maxE)) printf("TabulatedInelasticSM::simESEf: minE (%.10e) <= maxE (%.10e)\n", minE, maxE);
      double Exy = (minE * (1. - r)) + (maxE * r);
      double ESEf = Exy + Ezq;
      double theta = ::acos(kzf / ::sqrt(ESEf));
      VectorXd res(2);
      res[0] = ESEf;
      res[1] = theta;
      return res;
   }

   static double computeE0fromDispersion(double Eq, double deltaE)
   {
      /*
      * Note to self: The derivation of this algorithm is in
      * Y:\proj\linewidth\jvillar
      * \develop\NewMONSEL\Physics\DielectricDevelopment\
      * SimulatingALADingShimizu.nb
      */
      if (Eq == 0.)
         return deltaE;
      /* Precompute quantities we're going to need more than once. */
      double x = Eq / deltaE;
      /*
      * The numeric constant on the next line is a constant that appears in the
      * solution of Penn's dispersion equation. It has units of 1/Joule.
      */
      double y = 3.77300614251479e17 * Eq;
      double x13 = ::pow(x, 1 / 3.);
      double x2 = x * x;
      double c1 = x2 * x2 * (27. + (18. * y) + (2 * y * y));
      double c2 = x2 * (27. + (4. * y));
      double c3 = c2 - 27.;
      double c4 = x2 * (3. + y);
      double c6 = 1 - x2;
      if (c3 > 0.) {
         double tanPart = ::atan2(3. * c6 * ::sqrt(3. * c3 * c6), (18. * c4) - c1 - 27.);
         double trigPart = ::cos(tanPart / 3.);
         double prefactor = 2. * x * ::sqrt(y * (y - (c6 * (y + 6.))));
         trigPart *= prefactor;
         double result = deltaE * ::sqrt(((3. - c4) + trigPart) / 3.);
         return result;
      }
      else {
         double c5 = ::pow(-c1 + (18. * c4) + (3. * (-9. + (c6 * ::sqrt(-(3. * c3 * c6))))), 1. / 3.);
         double term1 = -2. * c4;
         double x23 = x13 * x13;
         double x43 = x * x13;
         double y13 = ::pow(y, 1. / 3.);
         double y23 = y13 * y13;
         double x43y23 = x43 * y23;
         double two13 = ::pow(2., 1. / 3.);
         double term2 = -12. * two13 * x43y23 * c6;
         double term3 = 2. * two13 * x43y23 * x2 * y;
         double term23 = (term2 + term3) / c5;
         double term4 = two13 * two13 * x23 * y13 * c5;
         double result = deltaE * ::sqrt((6 + term1 + term23 + term4) / 6.);
         return result;
      }
   }

   double TabulatedInelasticSM::pickBE(double Eq, double deltaE)
   {
      int i;
      /*
      * Detect and return immediately in the most common case (deltaE too small
      * for inner shell excitation)
      */
      if ((coreEnergies.size() > 0) && (deltaE <= coreEnergies[0]))
         return 0.;
      /*
      * Arrive here if there is enough energy in principle to free an inner
      * shell electron. In this case we must compute E0 from the dispersion
      * equation, given Eq and deltaE.
      */
      double energyForIonization;
      if (E0fromDispersion)
         energyForIonization = computeE0fromDispersion(Eq, deltaE);
      else
         energyForIonization = deltaE;

      for (i = 0; (i < coreEnergies.size()) && (coreEnergies[i] <= energyForIonization); i++)
         ;
      if (i == 0)
         return 0.;
      if (defaultRatios)
         /*
         * The advertised default behavior, as described by Ding & Shimizu
         * (Scanning).
         */
         return coreEnergies[i - 1];
      else {
         auto cprob = cumulativeBranchingProbabilities[i - 1];
         double r = Math2::random();
         int index = Algorithm::binarySearch(cprob.data(), 0, cprob.size()-1, r);
         /*
         * Above binary search returns a positive index in the rare case when r
         * is exactly in the table. In such a case the energy we want is
         * coreEnergies[index]. Usually, though, there is no exact match for r.
         * In such case binarySearch returns index=-insertionPoint-1, so
         * insertionPoint (the index of the first table value that is larger
         * than r) is -(index+1). The energy we want though, is the next
         * smaller one than this, i.e., the one at -(index+1)-1 = -index-2. If
         * this remains less than 0, it means ALL table entries were greater
         * than r, in which case the energy we want is 0.
         */

         if (index < 0) {
            index = -2 - index;
            if (index < 0)
               return 0.;
         }

         return coreEnergies[index];
      }
   }

   ElectronT* TabulatedInelasticSM::scatter(ElectronT* pe)
   {
      double kE0 = pe->getEnergy(); // PE initial energy rel to CB bottom
      double kE = kE0 + energyOffset; // PE initial energy rel to
      // scattering band bottom
      if (kE < tableEiDomain[0])
         /*
         * This might happen if something, e.g., electrostatic potential
         * difference, reduces electron energy between the time we determine
         * that scattering happens (at the beginning of a step) and the time it
         * actually occurs (at the end of the step). In this case, we simply
         * don't scatter.
         */
         return NULL;
      if (kE > tableEiDomain[1])
         printf("PE energy %.10e is outside the interpolation table interval of [%.10e, %.10e]\n", kE, tableEiDomain[0], tableEiDomain[1]);
      double theta = 0.;
      double phi = 0.; // PE trajectory parameters
      double energySE, thetaSE, phiSE; // SE trajectory parameters
      // TODO Do I need to check that kE>offsetFermiEnergy?
      double randoms[] = { Math2::random(), Math2::random(), Math2::random(), Math2::random() };
      interpInput[0] = kE;
      interpInput[1] = randoms[0];
      // Energy loss by PE
      double deltaE = kE * tableReducedDeltaE->interpolate(interpInput.data(), interpInput.size(), 3);
      /*
      * Cubic interpolation of the table can undershoot. Treat deltaE close to
      * but below the energyGap as such undershoot and correct it.
      */
      if ((deltaE < energyGap) && (deltaE >(0.95 * energyGap)))
         deltaE = energyGap;
      /*
      * Larger discrepancies are most likely because we've been supplied an
      * empirical table that includes non-electronic energy losses (e.g.,
      * scattering from phonons). These should really be handled separately
      * because our model is only valid for electrons & plasmons. (E.g., phonon
      * of energy deltaE carries very different momentum from electron of
      * energy deltaE, so scattering angles can't be determined in the present
      * model.) We skip the angular scattering part for such events. Any
      * generated SE will be in the bandgap, so most likely dropped anyway. We
      * return after we deal with the PE energy loss.
      */

      double theta0PE = pe->getTheta(); // Remember original direction;
      double phi0PE = pe->getPhi(); // to use for SE
      if (deltaE >= bandgap) {
         // Determine theta and phi here
         /*
         * First, the reduced energy. This parameter ranges from 0 to 1 as
         * deltaE ranges from its minimium to maximum value.
         */
         interpInput[1] = (deltaE - energyGap) / (kE - offsetFermiEnergy - (2. * energyGap));
         /*
         * The reduced energy can on rare occasions, as a result of
         * interpolation error, lie slightly outside its physically determined
         * interval of [0,1]. If it does, clip it to the boundary.
         */
         if (interpInput[1] > 1.)
            interpInput[1] = 1.;
         else if (interpInput[1] < 0.)
            interpInput[1] = 0.;
         interpInput[2] = randoms[1];
         theta = tableTheta->interpolate(interpInput.data(), interpInput.size(), 3);
         phi = 2. * Math2::PI * randoms[2];
         /*
         * Update PE trajectory. Note that the energy of the PE is decremented
         * by deltaE. Any continuous energy loss formula should only account
         * for losses not including this one.
         */
         pe->updateDirection(theta, phi);
      }

      pe->setEnergy(kE0 - deltaE);

      /*
      * I originally reset the previous energy to kE0 (next line), but I'm now
      * commenting it though I keep it here as a place-marker. My thinking is
      * that pe.getPreviousEnergy() should return the energy at the end of the
      * last step. In takeStep(), that energy has already been stored as the
      * previous energy, when takeStep() took account of any CSD energy loss.
      * Any new energy loss in the current routine should simply reduce the
      * current energy, so that pe.getEnergy()-pe.getPreviousEnergy() should
      * reflect the total energy loss this step, CSD + inelastic mechanisms
      * combined.
      */

      // Determine SE final energy and trajectory
      ElectronT* se = NULL;
      double be = 0.;

      /*
      * Some measured ELF data have nonzero values for deltaE less than the
      * bandgap, and it is consequently possible that some scattering tables
      * that retain this part of the data will include such loss events. These
      * may, for example, correspond to phonon losses. They presumably do not
      * correspond to generation of mobile SE, since there are no empty mobile
      * states in the gap. We therefore return no SE for such events.
      */
      if (deltaE < bandgap)
         return NULL;

      double Eq = (2. * kE) - deltaE - (2. * ::sqrt(kE * (kE - deltaE)) * ::cos(theta));

      switch (methodSE) {
      case 1:
         /*
         * In the following formula, offsetFermiEnergy - energyOffset is the
         * Fermi energy re-referenced to the bottom of the conduction band.
         * If b (the nearest lower core level binding energy) is zero, then
         * this is the electron's initial energy. (mode 1 assumes target
         * electrons come from the Fermi level.) If b>0 then the Fermi
         * energy - b is still the electron's initial energy, since the core
         * level energies are referenced to the Fermi level. Either way,
         * adding deltaE gives the SE's final energy.
         */
         energySE = (deltaE + bEref) - pickBE(Eq, deltaE);
         if ((energySE + energyCBbottom) < minEgenSE)
            return NULL;
         thetaSE = Math2::PI / 2. - theta;
         phiSE = phi + Math2::PI;
         // Generate SE, apply energy loss and trajectory change to SE here
         se = new ElectronT(*pe, theta0PE, phi0PE, energySE);
         se->updateDirection(thetaSE, phiSE);
         break;
      case 2:
         be = pickBE(Eq, deltaE);
         if (be > 0.)
            energySE = (deltaE + bEref) - be;
         else {
            interpInput[0] = deltaE;
            interpInput[1] = randoms[3];
            double energy0SE = tableSEE0->interpolate(interpInput.data(), interpInput.size(), 3);
            /*
            * The values in the SEE0 table should range from 0 to EFermi,
            * which represents the range of allowed values. If the
            * interpolated value overshoots, clip it.
            */
            if (energy0SE < energyRangeSE0[0])
               energy0SE = energyRangeSE0[0];
            else if (energy0SE > energyRangeSE0[1])
               energy0SE = energyRangeSE0[1];
            energySE = (deltaE + energy0SE) - energyOffset;
         }
         if ((energySE + energyCBbottom) < minEgenSE)
            return NULL;
         thetaSE = Math2::PI / 2. - theta;
         phiSE = phi + Math2::PI;
         // Generate SE, apply energy loss and trajectory change to SE here
         se = new ElectronT(*pe, theta0PE, phi0PE, energySE);
         se->updateDirection(thetaSE, phiSE);
         break;
      case 3:
         be = pickBE(Eq, deltaE);
         if (be > 0.) { // core level excitation
            energySE = (deltaE + bEref) - be;
            if ((energySE + energyCBbottom) < minEgenSE)
               return NULL;
            /*
            * I'm going to approximate the angle distribution as isotropic
            * for now.
            */
            thetaSE = ::acos(1. - (2. * Math2::random()));
            phiSE = 2. * Math2::PI * Math2::random();
            // Generate SE, apply energy loss and trajectory change to SE
            // here
            se = new ElectronT(*pe, thetaSE, phiSE, energySE);
         }
         else { // SE generation from extended band
            double root = 2. * ::sqrt(offsetFermiEnergy * (offsetFermiEnergy + deltaE));
            double sum = (2. * offsetFermiEnergy) + deltaE;
            double Eqmin = sum - root;
            double Eqmax = sum + root;
            if ((Eqmin <= Eq) && (Eq <= Eqmax)) { // single-electron
               // scattering
               auto energytheta = simESEf(Eq, deltaE, randoms[3]);
               energySE = energytheta[0] - energyOffset;
               if ((energySE + energyCBbottom) < minEgenSE)
                  return NULL;
               // Generate SE in PE direction with correct energy
               se = new ElectronT(*pe, theta0PE, phi0PE, energySE);
               // Determine angles of SE q vector relative to PE original direction
               thetaSE = Math2::PI / 2. - theta;
               phiSE = phi + Math2::PI;
               // Combine with adjustment for additional simESEf deflection
               auto newdir = updateDirection(thetaSE, phiSE, energytheta[1], 2. * Math2::PI * Math2::random());
               //Update SE direction by this combined amount
               se->updateDirection(newdir[0], newdir[1]);

            }
            else { // plasmon scattering
               interpInput[0] = deltaE;
               interpInput[1] = randoms[3];
               double energy0SE = tableSEE0->interpolate(interpInput.data(), interpInput.size(), 3);
               /*
               * The values in the SEE0 table should range from 0 to EFermi,
               * which represents the range of allowed values. If the
               * interpolated value overshoots, clip it.
               */
               if (energy0SE < energyRangeSE0[0])
                  energy0SE = energyRangeSE0[0];
               else if (energy0SE > energyRangeSE0[1])
                  energy0SE = energyRangeSE0[1];
               energySE = (deltaE + energy0SE) - energyOffset;
               if ((energySE + energyCBbottom) < minEgenSE)
                  return NULL;
               /*
               * For plasmon scattering, mode 3 assumes the plasmon
               * "forgets" the momentum of the event that created it before
               * it decays into an electron-hole pair. The angular
               * distribution is therefore isotropic.
               */
               thetaSE = ::acos(1. - (2. * Math2::random()));
               phiSE = 2 * Math2::PI * Math2::random();
               // Generate SE, apply energy loss and trajectory change to SE
               // here
               se = new ElectronT(*pe, thetaSE, phiSE, energySE);
            }
         }
         break;
      default:
         se = NULL;
         break;
      }

      return se;
   }

   double TabulatedInelasticSM::scatterRate(const ElectronT* pe)
   {
      kEa[0] = pe->getEnergy() + energyOffset; // The PE kinetic energy
      /*
      * The PE kinetic energy can fall below the minimum in the table for
      * materials with a energyGap. In this case the actual scatter rate is 0.
      */
      if (kEa[0] < tableIIMFPEiDomain[0])
         return 0.;
      if (kEa[0] > tableIIMFPEiDomain[1])
         printf("TabulatedInelasticSM::scatterRate: PE energy %.10e exceeds interpolation table maximum energy of %.10e\n", kEa[0], tableIIMFPEiDomain[1]);

      /*
      * I do only first order interpolation below because I noticed for some
      * tables that I get negative interpolated values. This happens despite
      * having all positive values in the table, because the scatter rate
      * approaches 0 at the Fermi level, leaving open the possibility of
      * overshoot. Possible approaches to avoid overshoot are to use linear
      * interpolation or to clip the result (as I do below for other tables).
      * Clipping to 0 seems a bad choice here, because it results in an
      * infinite inelastic free path.
      */
      return rateMult * tableIIMFP->interpolate(kEa.data(), kEa.size(), 1);
   }

   void TabulatedInelasticSM::setMaterial(const MaterialT* mat)
   {
      if (!(mat->isSEmaterial()))
         printf("Material %s is not an SEmaterial as required for TabulatedInelasticSM.", mat->toString());

      const SEmaterialT* semat = (SEmaterialT*)mat;

      energyCBbottom = semat->getEnergyCBbottom();
      workfunction = semat->getWorkfunction();
      bandgap = semat->getBandgap();
      energyGap = bandgap;
      offsetFermiEnergy = semat->getEFermi() + energyOffset;
      auto coreEnergiesSet = semat->getCoreEnergyArray();
      for (auto ce : coreEnergiesSet) {
         coreEnergies.push_back(ce);
      }
      //coreEnergies.insert(coreEnergies.begin(), coreEnergiesSet.begin(), coreEnergiesSet.end());
      /*
      * My source for binding energies provides them relative to vacuum for a
      * few nobel and binary gases that we're unlikely to use, relative to the
      * Fermi level for metals, and relative to the top of the valence band for
      * insulators and semiconductors. bEref gives us our offsets to bottom of
      * conduction band.
      */
      if (bandgap > 0.)
         bEref = -bandgap;
      else
         bEref = semat->getEFermi();

      /*
      * tableEiDomain must be the energy range that is valid for *all* required
      * tables that take the PE initial energy as an input parameter.
      */
      tableEiDomain = tableReducedDeltaE->getDomain()[0];
      auto thetaTableEiDomain = tableTheta->getDomain()[0];
      if (thetaTableEiDomain[0] > tableEiDomain[0])
         tableEiDomain[0] = thetaTableEiDomain[0];
      if (thetaTableEiDomain[1] < tableEiDomain[1])
         tableEiDomain[1] = thetaTableEiDomain[1];
      tableIIMFPEiDomain = tableIIMFP->getDomain()[0];
   }

   void TabulatedInelasticSM::setMinEgenSE(double minEgenSE)
   {
      if (minEgenSE > -workfunction)
         this->minEgenSE = minEgenSE;
      else
         printf("TabulatedInelasticSM::setMinEgenSE: Illegal minEgenSE %.10e.", minEgenSE);
   }

   /**
   * @return Returns the minEgenSE.
   */
   double TabulatedInelasticSM::getMinEgenSE() const
   {
      return minEgenSE;
   }

   /**
   * Returns the value of the SE method parameter that is being used. This is
   * normally the same as the value supplied to the constructor.
   *
   * @return the methodSE
   */
   int TabulatedInelasticSM::getMethodSE() const
   {
      return methodSE;
   }

   /**
   * Branching ratios control how this class associates a core (binding) energy
   * with an excitation. If deltaE is the energy lost by the primary electron
   * in a scattering event, the secondary electron's final energy is equal to
   * its initial energy plus deltaE, but what is its initial energy? Was it
   * originally in the valence band, or in a deeper bound state? Only initial
   * states with binding energies less than deltaE contribute a channel for
   * energy loss. Typically the energy loss function (ELF) has a discontinuous
   * increase when deltE is equal to its energy. TabulatedInelasticSM assumes
   * the excitation channel that does not involve one of these bound electrons
   * is unaffected by the presence of the new excitation channel. Thus, the
   * ratio, r, of the ELF value from the left to the ELF value from the right
   * represents the probability that the new channel is not involved, and 1-r
   * is the probability that it is. These r ratios can easily be determined
   * from the ELF vs. energy curve. ratios (the argument of this method) is an
   * array of these ratios. The first ratio in the array is associated with the
   * lowest nonzero core energy (i.e., the first non-valence band bound state).
   * The length of the array must be equal to the length of the material's
   * coreEnergies array.
   * </p>
   * <p>
   * The default behavior, if this method is not called or if it is called with
   * no argument, is to assume all entries are 0. That is, the largest eligible
   * binding energy state is assumed to be the one associated with the
   * excitation channel. This is the method described by Ding & Shimizu in
   * SCANNING.
   */
   void TabulatedInelasticSM::setBranchingRatios()
   {
      defaultRatios = true;
   }

   void TabulatedInelasticSM::setBranchingRatios(double ratios[], int len)
   {
      defaultRatios = false;
      int cElen = coreEnergies.size();
      if (len != cElen)
         printf("The number of branching ratios must be equal to the number of core energies, %d, in this case.", cElen);
      MatrixXd probabilities(cElen);
      cumulativeBranchingProbabilities.resize(cElen);
      probabilities[0].push_back(ratios[0]);
      cumulativeBranchingProbabilities[0].push_back(ratios[0]);
      for (int i = 1; i < cElen; i++) {
         probabilities[i].resize(i + 1);
         cumulativeBranchingProbabilities[i].resize(i + 1);
         for (int j = 0; j < i; j++)
            probabilities[i][j] = probabilities[i - 1][j] * ratios[i];
         probabilities[i][i] = (1. - ratios[i - 1]) * ratios[i];
         cumulativeBranchingProbabilities[i][0] = probabilities[i][0];
         for (int j = 1; j <= i; j++)
            cumulativeBranchingProbabilities[i][j] = cumulativeBranchingProbabilities[i][j - 1] + probabilities[i][j];
      }
   }

   void TabulatedInelasticSM::setEnergyGap(double energyGap)
   {
      this->energyGap = energyGap;
   }

   /**
   * The scatterRate or IIMFP (inverse inelastic mean free path) derived from
   * the associated table is multiplied by the factor provided via this routine
   * (default = 1). I'm not sure if I'm going to keep this. I added it as an ad
   * hoc way to simulate the effect of porosity in a material. The multiplier
   * can be set equal to the fill fraction, i.e., the fraction of the space
   * actually occupied by the material, in order to simulate the effect of many
   * tiny pores in the material.
   *
   * @param rateMult
   */
   void TabulatedInelasticSM::setRateMult(double rateMult)
   {
      this->rateMult = rateMult;
   }

   /**
   * Gets the current value assigned to E0fromDispersion
   *
   * @return Returns the E0fromDispersion.
   */
   bool TabulatedInelasticSM::isE0fromDispersion() const
   {
      return E0fromDispersion;
   }

   /**
   * Sets the value assigned to E0fromDispersion. If E0fromDispersion = false
   * (the default) JMONSEL continutes to use its original method (Ding &
   * Shimizu's SCANNING method) for determining the energy available to ionize
   * an inner shell. This method assumes the shell may be ionized whenever
   * deltaE (the energy transferred to the SE in an inelastic event) is greater
   * than the ionization energy. If E0fromDispersion = true, it uses Ding et
   * al's later method, in which 0-momentum part (E0) of the energy is computed
   * from the plasmon dispersion for an event which transfers deltaE. Inner
   * shell ionization can only happen if E0 > ionization energy. This is more
   * restrictive than deltaE > ionization energy.
   *
   * @param e0fromDispersion The value to which to set E0fromDispersion.
   */
   void TabulatedInelasticSM::setE0fromDispersion(bool e0fromDispersion)
   {
      E0fromDispersion = e0fromDispersion;
   }
}