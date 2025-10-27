#include "PolycrystalDensityAux.h"

registerMooseObject("mistApp", PolycrystalDensityAux);

InputParameters
PolycrystalDensityAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("read density from CSV using tessellation of density_i");
  params.addRequiredParam<std::string>("csv_density", "csv_density");
  params.addRequiredCoupledVar("density_i", "density_i");
  params.addRequiredParam<Real>("density_scaling", "density_scaling");
  return params;
}

PolycrystalDensityAux::PolycrystalDensityAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _csv_density(getParam<std::string>("csv_density")),
    _density_i(coupledValue("density_i")),
    _density_scaling(getParam<Real>("density_scaling"))
{
    _csv_total_density = readCSV(_csv_density);
}

Real
PolycrystalDensityAux::computeValue()
{
    std::vector<Real> first_col = _csv_total_density[0];
    Real pred_density = first_col.at(_density_i[_qp]);
    return _density_scaling * pred_density;
}

std::vector<std::vector<Real>>
PolycrystalDensityAux::readCSV(const std::string csv_file_name){
  std::ifstream file(csv_file_name);
  if (!file.is_open()){
    mooseError("can't open CSV file");
  }

  std::vector<std::vector<Real>> data;
  std::string line;
  while (std::getline(file, line)){
    std::vector<Real> row;
    std::stringstream ss(line);
    std::string value;

    while (std::getline(ss, value, ',')){
      try{
        row.push_back(static_cast<Real>(std::stod(value)));
      }
      catch (const std::invalid_argument &e){
        mooseError("Invalid value in CSV file: " + value);
      }
    }
    data.push_back(row);
  }
  file.close();
  return data;
}