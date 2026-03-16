#include "CSVFractionIC.h"
#include "FEProblem.h"           // ✅ Full definition of FEProblemBase
#include "MooseMesh.h"
#include "AuxiliarySystem.h"
#include "libmesh/point_locator_base.h"
#include "libmesh/system.h"       // ✅ Full libMesh System type

registerMooseObject("mistApp", CSVFractionIC);

InputParameters
CSVFractionIC::validParams()
{
  InputParameters params = InitialCondition::validParams();
  params.addRequiredParam<VariableName>(
      "fraction_csv",
      "Name of the AuxVariable (elemental) that stores the fraction values from CSV.");
  return params;
}

CSVFractionIC::CSVFractionIC(const InputParameters & parameters)
  : InitialCondition(parameters),
    _fraction_var(getParam<VariableName>("fraction_csv"))
{
}

Real
CSVFractionIC::value(const Point & p)
{
  const libMesh::PointLocatorBase & locator = *_fe_problem.mesh().getPointLocator();
  const Elem * elem = locator(p);

  if (!elem)
  {
    // Out of domain (node may lie on partition boundary)
    return 0.0;
  }

  const auto & aux_sys = _fe_problem.getAuxiliarySystem().system();
  const unsigned int var_num = aux_sys.variable_number(_fraction_var);
  const auto & dof_map = aux_sys.get_dof_map();

  std::vector<dof_id_type> dof_indices;
  dof_map.dof_indices(elem, dof_indices, var_num);

  if (dof_indices.empty())
  {
    // No local DOFs on this processor for that element
    return 0.0;
  }

  const auto & vec = *aux_sys.current_local_solution;

  // Make sure index is valid for this processor
  const dof_id_type dof = dof_indices[0];
  if (dof >= vec.size())
    return 0.0;

  const Real val = vec.el(dof);

  // Optional safety clamp
  if (!std::isfinite(val))
    return 0.0;

  return val;
}
