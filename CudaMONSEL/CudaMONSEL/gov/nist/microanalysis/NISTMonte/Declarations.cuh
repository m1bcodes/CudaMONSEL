#ifndef _DECLARATIONS_CUH_
#define _DECLARATIONS_CUH_

//#include <string>
#include <set>
#include <map>
#include <unordered_map>
#include <stack>

#include "Amphibian\String.cuh"
#include "Amphibian\vector.cuh"

typedef amp::string StringT;
typedef amp::vector<float> VectorXf;
typedef amp::vector<VectorXf> MatrixXf;
typedef amp::vector<MatrixXf> Matrix3DXf;
typedef amp::vector<Matrix3DXf> Matrix4DXf;
typedef amp::vector<double> VectorXd;
typedef amp::vector<VectorXd> MatrixXd;
typedef amp::vector<MatrixXd> Matrix3DXd;
typedef amp::vector<Matrix3DXd> Matrix4DXd;
typedef amp::vector<int> VectorXi;
typedef amp::vector<VectorXi> MatrixXi;

namespace Element
{
   class Element;
}

typedef Element::Element ElementT;

namespace Composition
{
   class Composition;
}

typedef Composition::Composition CompositionT;

namespace Material
{
   class Material;
}

typedef Material::Material MaterialT;

namespace Electron
{
   class Electron;
}

typedef Electron::Electron ElectronT;

namespace RegionBase
{
   class RegionBase;
}

typedef RegionBase::RegionBase RegionBaseT;

namespace RegionBase
{
   class TransformableRegion;
}

typedef RegionBase::TransformableRegion TransformableRegionT;

namespace IMaterialScatterModel
{
   class IMaterialScatterModel;
}

typedef IMaterialScatterModel::IMaterialScatterModel IMaterialScatterModelT;

namespace Shape
{
   class Shape;
}

typedef Shape::Shape ShapeT;

namespace ITransform
{
   class ITransform;
}

typedef ITransform::ITransform ITransformT;

namespace RegionBase
{
   class Region;
}

typedef RegionBase::Region RegionT;

namespace ElectronGun
{
   class ElectronGun;
}

typedef ElectronGun::ElectronGun ElectronGunT;

namespace NullMaterialScatterModel
{
   class NullMaterialScatterModel;
}

typedef NullMaterialScatterModel::NullMaterialScatterModel NullMaterialScatterModelT;

namespace GaussianBeam
{
   class GaussianBeam;
}

typedef GaussianBeam::GaussianBeam GaussianBeamT;

namespace Sphere
{
   class Sphere;
}

typedef Sphere::Sphere SphereT;

namespace Reference
{
   class Reference;
}

typedef Reference::Reference ReferenceT;

namespace AlgorithmUser
{
   class AlgorithmUser;
}

typedef AlgorithmUser::AlgorithmUser AlgorithmUserT;

namespace AlgorithmClass
{
   class AlgorithmClass;
}

typedef AlgorithmClass::AlgorithmClass AlgorithmClassT;

namespace RandomizedScatter
{
   class RandomizedScatter;
}

typedef RandomizedScatter::RandomizedScatter RandomizedScatterT;

namespace RandomizedScatterFactory
{
   class RandomizedScatterFactory;
}

typedef RandomizedScatterFactory::RandomizedScatterFactory RandomizedScatterFactoryT;

namespace ScreenedRutherfordRandomizedScatterFactory
{
   class ScreenedRutherfordRandomizedScatterFactory;
}

typedef ScreenedRutherfordRandomizedScatterFactory::ScreenedRutherfordRandomizedScatterFactory ScreenedRutherfordRandomizedScatterFactoryT;

namespace CzyzewskiMottCrossSection
{
   class CzyzewskiMottCrossSection;
}

typedef CzyzewskiMottCrossSection::CzyzewskiMottCrossSection CzyzewskiMottCrossSectionT;

namespace ScreenedRutherfordScatteringAngle
{
   class ScreenedRutherfordScatteringAngle;
}

typedef ScreenedRutherfordScatteringAngle::ScreenedRutherfordScatteringAngle ScreenedRutherfordScatteringAngleT;

namespace Strategy
{
   class Strategy;
}

typedef Strategy::Strategy StrategyT;

namespace AtomicShell
{
   class AtomicShell;
}

typedef AtomicShell::AtomicShell AtomicShellT;

namespace EdgeEnergy
{
   class EdgeEnergy;
}

typedef EdgeEnergy::EdgeEnergy EdgeEnergyT;

namespace BarrierScatterMechanism
{
   class BarrierScatterMechanism;
}

typedef BarrierScatterMechanism::BarrierScatterMechanism BarrierScatterMechanismT;

namespace BrowningEmpiricalCrossSection
{
   class BrowningEmpiricalCrossSection;
}
typedef BrowningEmpiricalCrossSection::BrowningEmpiricalCrossSection BrowningEmpiricalCrossSectionT;

namespace SimpleBlock
{
   class SimpleBlock;
}

typedef SimpleBlock::SimpleBlock SimpleBlockT;

namespace ScatterMechanism
{
   class ScatterMechanism;
}

typedef ScatterMechanism::ScatterMechanism ScatterMechanismT;

namespace SlowingDownAlg
{
   class SlowingDownAlg;
}

typedef SlowingDownAlg::SlowingDownAlg SlowingDownAlgT;

namespace MeanIonizationPotential
{
   class MeanIonizationPotential;
}

typedef MeanIonizationPotential::MeanIonizationPotential MeanIonizationPotentialT;

namespace NUTableInterpolation
{
   class NUTableInterpolation;
}

typedef NUTableInterpolation::NUTableInterpolation NUTableInterpolationT;

namespace MultiPlaneShape
{
   class MultiPlaneShape;
}

typedef MultiPlaneShape::MultiPlaneShape MultiPlaneShapeT;

namespace MultiPlaneShape
{
   class Plane;
}

typedef MultiPlaneShape::Plane PlaneT;

namespace CylindricalShape
{
   class CylindricalShape;
}

typedef CylindricalShape::CylindricalShape CylindricalShapeT;

namespace NormalCylindricalShape
{
   class NormalCylindricalShape;
}

typedef NormalCylindricalShape::NormalCylindricalShape NormalCylindricalShapeT;

namespace SumShape
{
   class SumShape;
}

typedef SumShape::SumShape SumShapeT;

namespace NormalUnionShape
{
   class NormalUnionShape;
}

typedef NormalUnionShape::NormalUnionShape NormalUnionShapeT;

namespace NormalComplementShape
{
   class NormalComplementShape;
}

typedef NormalComplementShape::NormalComplementShape NormalComplementShapeT;

namespace NormalDifferenceShape
{
   class NormalDifferenceShape;
}

typedef NormalDifferenceShape::NormalDifferenceShape NormalDifferenceShapeT;

namespace BetheElectronEnergyLoss
{
   class BetheElectronEnergyLoss;
}

typedef BetheElectronEnergyLoss::BetheElectronEnergyLoss BetheElectronEnergyLossT;

namespace BasicMaterialModel
{
   class BasicMaterialModel;
}

typedef BasicMaterialModel::BasicMaterialModel BasicMaterialModelT;

namespace Histogram
{
   class Histogram;
}

typedef Histogram::Histogram HistogramT;

namespace MonteCarloSS
{
   class MonteCarloSS;
}

typedef MonteCarloSS::MonteCarloSS MonteCarloSST;

namespace Histogram
{
   class Histogram;
}

typedef Histogram::Histogram HistogramT;

namespace ActionListener
{
   class ActionListener;
}

typedef ActionListener::ActionListener ActionListenerT;

namespace BackscatterStats
{
   class BackscatterStats;
}

typedef BackscatterStats::BackscatterStats BackscatterStatsT;

namespace MultiPlaneShape
{
   struct LineShape;
}

typedef MultiPlaneShape::LineShape LineShapeT;

// nano

namespace SEmaterial
{
   class SEmaterial;
}

typedef SEmaterial::SEmaterial SEmaterialT;

namespace ExpQMBarrierSM
{
   class ExpQMBarrierSM;
}

typedef ExpQMBarrierSM::ExpQMBarrierSM ExpQMBarrierSMT;

namespace MONSEL_MaterialScatterModel
{
   class MONSEL_MaterialScatterModel;
}

typedef MONSEL_MaterialScatterModel::MONSEL_MaterialScatterModel MONSEL_MaterialScatterModelT;

namespace SelectableElasticSM
{
   class SelectableElasticSM;
}

typedef SelectableElasticSM::SelectableElasticSM SelectableElasticSMT;

namespace JoyLuoNieminenCSD
{
   class JoyLuoNieminenCSD;
}

typedef JoyLuoNieminenCSD::JoyLuoNieminenCSD JoyLuoNieminenCSDT;

namespace FittedInelSM
{
   class FittedInelSM;
}

typedef FittedInelSM::FittedInelSM FittedInelSMT;

namespace GanachaudMokraniPolaronTrapSM
{
   class GanachaudMokraniPolaronTrapSM;
}

typedef GanachaudMokraniPolaronTrapSM::GanachaudMokraniPolaronTrapSM GanachaudMokraniPolaronTrapSMT;

namespace TabulatedInelasticSM
{
   class TabulatedInelasticSM;
}

typedef TabulatedInelasticSM::TabulatedInelasticSM TabulatedInelasticSMT;

namespace GanachaudMokraniPhononInelasticSM
{
   class GanachaudMokraniPhononInelasticSM;
}

typedef GanachaudMokraniPhononInelasticSM::GanachaudMokraniPhononInelasticSM GanachaudMokraniPhononInelasticSMT;

namespace ZeroCSD
{
   class ZeroCSD;
}

typedef ZeroCSD::ZeroCSD ZeroCSDT;


namespace NormalShape
{
   class NormalShape;
}

typedef NormalShape::NormalShape NormalShapeT;

namespace NormalMultiPlaneShape
{
   class NormalMultiPlaneShape;
}

typedef NormalMultiPlaneShape::NormalMultiPlaneShape NormalMultiPlaneShapeT;

namespace NormalIntersectionShape
{
   class NormalIntersectionShape;
}

typedef NormalIntersectionShape::NormalIntersectionShape NormalIntersectionShapeT;

#endif