#include "PolycrystalFractionAux.h"

registerMooseObject("mlApp", PolycrystalFractionAux);

InputParameters
PolycrystalFractionAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("read initial RDX fraction from CSV using tessellation of density_i");
  params.addRequiredParam<std::string>("csv_fraction", "csv_fraction");
  params.addRequiredParam<std::string>("csv_fraction_pore", "csv_fraction_pore");
  params.addRequiredCoupledVar("density_i", "density_i");

  //we need to add an input that is consistent with the placeholder MicroID assigned in the UO
  params.addRequiredParam<unsigned int>("bulk_MicroID", "bulk_MicroID");
  params.addRequiredParam<Real>("bulk_RDX_fraction", "bulk_RDX_fraction");
  params.addRequiredParam<std::vector<unsigned int>>("range_pore", "range_pore");

  //for using loaded microstructure
  params.addParam<bool>("use_loaded_microstructure", false, "use_loaded_microstructure");
  params.addCoupledVar("loaded_microstructure", "loaded_microstructure");

  //for segmentation
  params.addParam<Real>("pore_limit", 0.25, "pore_limit");
  params.addParam<Real>("bulk_limit", 1.77, "bulk_limit");
  return params;
}

PolycrystalFractionAux::PolycrystalFractionAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _csv_fraction(getParam<std::string>("csv_fraction")),
    _csv_fraction_pore(getParam<std::string>("csv_fraction_pore")),
    _density_i(coupledValue("density_i")),

    //for bulk assignment
    _bulk_MicroID(getParam<unsigned int>("bulk_MicroID")),
    _bulk_RDX_fraction(getParam<Real>("bulk_RDX_fraction")),
    _range_pore(getParam<std::vector<unsigned int>>("range_pore")),
    _use_loaded_microstructure(getParam<bool>("use_loaded_microstructure")),
    _loaded_microstructure(isCoupled("loaded_microstructure") ? &coupledValue("loaded_microstructure") : nullptr),
    _pore_limit(getParam<Real>("pore_limit")),
    _bulk_limit(getParam<Real>("bulk_limit"))
{
    _csv_total_fraction = readCSV(_csv_fraction);
    _csv_total_fraction_pore = readCSV(_csv_fraction_pore);
}

Real
PolycrystalFractionAux::computeValue()
{
    const std::vector<Real> first_col = _csv_total_fraction[0];
    const std::vector<Real> first_col_pore = _csv_total_fraction_pore[0];

    bool is_pore = false;
    bool is_bulk = false;

    is_pore = (_density_i[_qp] >= _range_pore[0] && _density_i[_qp] <= _range_pore[1]);
    is_bulk = (_density_i[_qp] == _bulk_MicroID);

    //retrieve loaded microstructure if needed
    //if (_use_loaded_microstructure){
    //  //use the loaded density to define cases
    //  is_pore = ((*_loaded_microstructure)[_qp] <= _pore_limit);
    //  is_bulk = ((*_loaded_microstructure)[_qp] >= _bulk_limit);
    //}

    //before assigning fraction values, we need to evaluate the MicroID that was assigned by the UO
    Real pred_fraction;

    if (is_pore){
      const int call_density = static_cast<int>(std::round(_density_i[_qp])) - 100;
      pred_fraction = first_col_pore.at(call_density);
    }
    else if (is_bulk){
      pred_fraction = _bulk_RDX_fraction;
    }
    else { //outside, microPBXs
      pred_fraction = first_col.at(_density_i[_qp]);
    }
    return pred_fraction;
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