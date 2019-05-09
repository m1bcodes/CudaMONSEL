#include "gov\nist\nanoscalemetrology\JMONSEL\TabulatedInelasticSM.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\NUTableInterpolation.cuh"
#include "gov\nist\nanoscalemetrology\JMONSEL\SEmaterial.cuh"

namespace TabulatedInelasticSM
{
   //TabulatedInelasticSM::TabulatedInelasticSM(const SEmaterialT* mat, int methodSE, StringT tables[]) : ScatterMechanismT(*mat, methodSE, tables, 0.)
   //{
   //}

   TabulatedInelasticSM::TabulatedInelasticSM(const SEmaterialT& mat, int methodSE, StringT tables[], int tableslen, double energyOffset) :
      methodSE((methodSE != 2) && (methodSE != 3) ? 1 : methodSE),
      energyOffset(energyOffset),
      tableIIMFP(NUTableInterpolation::getInstance(tables[0].c_str())),
      tableReducedDeltaE(NUTableInterpolation::getInstance(tables[1].c_str())),
      tableTheta(NUTableInterpolation::getInstance(tables[2].c_str())),
      tableSEE0((methodSE == 2) || (methodSE == 3) ? NUTableInterpolation::getInstance(tables[3].c_str()) : NULL),
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
}