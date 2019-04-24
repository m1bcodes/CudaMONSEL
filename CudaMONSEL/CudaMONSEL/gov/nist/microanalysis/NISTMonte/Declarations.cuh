#ifndef _DECLARATIONS_CUH_
#define _DECLARATIONS_CUH_

#include <string>
#include <unordered_set>

typedef std::vector<double> PositionVecT;

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

namespace Element
{
   class Element;
}

typedef Element::Element ElementT;

//typedef std::unordered_set<ElementT, HashFcn> UnorderedSetT;
//typedef std::unordered_set<ElementT, HashFcn> OrderedSetT;

typedef std::string StringT;

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

#endif