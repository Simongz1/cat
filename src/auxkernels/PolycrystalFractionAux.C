#include "PolycrystalFractionAux.h"

registerMooseObject("mistApp", PolycrystalFractionAux);

InputParameters
PolycrystalFractionAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("read initial RDX fraction from CSV using tessellation of density_i");
  params.addRequiredParam<std::string>("csv_fraction", "csv_fraction");
  params.addRequiredCoupledVar("density_i", "density_i");
  return params;
}

PolycrystalFractionAux::PolycrystalFractionAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _csv_fraction(getParam<std::string>("csv_fraction")),
    _density_i(coupledValue("density_i"))
{
    _csv_total_fraction = readCSV(_csv_fraction);
}

Real
PolycrystalFractionAux::computeValue()
{
    const std::vector<Real> first_col = _csv_total_fraction[0];
    const Real pred_fraction = first_col.at(_density_i[_qp]);
    return  pred_fraction;
}

std::vector<std::vector<Real>>
PolycrystalFractionAux::readCSV(const std::string csv_file_name){
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