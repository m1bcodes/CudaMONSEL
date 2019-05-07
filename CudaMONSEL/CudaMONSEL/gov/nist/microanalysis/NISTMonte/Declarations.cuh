#ifndef _DECLARATIONS_CUH_
#define _DECLARATIONS_CUH_

#include <string>
#include <unordered_set>
#include <map>
#include <stack>

typedef std::string StringT;
typedef std::vector<double> PositionVecT;
typedef std::vector<float> VectorXf;
typedef std::vector<VectorXf> MatrixXf;
typedef std::vector<double> VectorXd;
typedef std::vector<VectorXd> MatrixXd;
typedef std::vector<int> VectorXi;
typedef std::vector<VectorXi> MatrixXi;
typedef std::unordered_set<double> Setd;

namespace Element
{
   class Element;
}

typedef Element::Element ElementT;

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

namespace SEmaterial
{
   class SEmaterial;
}

typedef SEmaterial::SEmaterial SEmaterialT;

namespace NormalShape
{
   class NormalShape;
}

typedef NormalShape::NormalShape NormalShapeT;

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

#endif