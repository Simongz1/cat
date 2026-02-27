#include "PolycrystalDensityIC.h"
#include "RandomIC.h"

registerMooseObject("mistApp", PolycrystalDensityIC);

InputParameters
PolycrystalDensityIC::validParams()
{
  InputParameters params = RandomIC::validParams();
  params.addClassDescription("use a voronoi tessellation to assign density_i");
  params.addRequiredParam<Real>("num_grains", "num_grains");
  params.addRequiredParam<Real>("min_in", "min_in");
  params.addRequiredParam<Real>("max_in", "max_in");
  params.addRequiredParam<Real>("min_out", "min_out");
  params.addRequiredParam<Real>("max_out", "max_out");
  params.addRequiredParam<Real>("target_grain", "target_grain");

  return params;
}

PolycrystalDensityIC::PolycrystalDensityIC(const InputParameters & parameters)
  : RandomIC(parameters),
  _num_grains(getParam<Real>("num_grains")),
  _min_in(getParam<Real>("min_in")),
  _max_in(getParam<Real>("max_in")),
  _min_out(getParam<Real>("min_out")),
  _max_out(getParam<Real>("max_out")),
  _target_grain(getParam<Real>("target_grain"))
{
}

void
PolycrystalDensityIC::initialSetup(){
  //get domain size
  const BoundingBox bbox = MeshTools::create_bounding_box(_fe_problem.mesh().getMesh());
  const Point min_corner = bbox.min();
  const Point max_corner = bbox.max();

  //generate clear array for centers
  _centers.clear();
  _centers.reserve(_num_grains);

  //assign random centers
  for (unsigned int i = 0; i < _num_grains; ++i){
    const Real x = MooseRandom::randl() * (max_corner(0) - min_corner(0)) + min_corner(0);
    const Real y = MooseRandom::randl() * (max_corner(1) - min_corner(1)) + min_corner(1);
    //in case of 3D
    const Real z = (_fe_problem.mesh().dimension() == 3) ? MooseRandom::randl() * (max_corner(2) - min_corner(2)) + min_corner(2) : 0.;
    _centers.emplace_back(x, y, z);
  }

  //output that the centers have been generated
  if (_tid == 0.){
    mooseInfo("generated", _num_grains, " centers for voronoi tessellation");
  }
}

Real
PolycrystalDensityIC::value(const Point & p)
{
  //get closest point for voronoi 
  Real min_dist = std::numeric_limits<Real>::max();

  unsigned int nearest = 0;

  for (unsigned int i = 0; i < _centers.size(); ++i){
    const Real d = (p - _centers[i]).norm();
    if (d < min_dist){
      min_dist = d;
      nearest = i;
    }
  }

  //up to this point, we have the grains assigned:
  const unsigned int  grain_id = nearest + 1;

  const Real rand_value = MooseRandom::randl();
  const bool in_grain = MooseUtils::absoluteFuzzyEqual(grain_id, _target_grain);

  if (in_grain){
    return _min_in + rand_value * (_max_in - _min_in);
  }
  else{
    return _min_out + rand_value * (_max_out - _min_out);
  }
}