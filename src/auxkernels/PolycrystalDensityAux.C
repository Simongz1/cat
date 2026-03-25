#include "PolycrystalDensityAux.h"

registerMooseObject("mlApp", PolycrystalDensityAux);

InputParameters
PolycrystalDensityAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("read density from CSV using tessellation of density_i");
  params.addRequiredParam<std::string>("csv_density", "csv_density");
  params.addRequiredCoupledVar("density_i", "density_i");
  params.addRequiredParam<Real>("density_scaling", "density_scaling");

  //parameters for bulk grains
  params.addRequiredParam<unsigned int>("bulk_MicroID", "bulk_MicroID");
  params.addRequiredParam<Real>("bulk_RDX_density", "bulk_RDX_density");
  params.addRequiredParam<std::string>("csv_density_pore", "csv_density_pore");
  params.addRequiredParam<std::vector<unsigned int>>("range_pore", "range_pore");

  //for using loaded microstructure
  params.addParam<bool>("use_loaded_microstructure", false, "use_loaded_microstructure");
  params.addCoupledVar("loaded_microstructure", "loaded_microstructure");

  //for segmentation
  params.addParam<Real>("pore_limit", 0.25, "pore_limit");
  params.addParam<Real>("bulk_limit", 1.77, "bulk_limit");
  return params;
}

PolycrystalDensityAux::PolycrystalDensityAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _csv_density(getParam<std::string>("csv_density")),
    _csv_density_pore(getParam<std::string>("csv_density_pore")),
    _density_i(coupledValue("density_i")),
    _density_scaling(getParam<Real>("density_scaling")),

    //parameters for bulk grains
    _bulk_MicroID(getParam<unsigned int>("bulk_MicroID")),
    _bulk_RDX_density(getParam<Real>("bulk_RDX_density")),

    //get params for pores, these are set as global params
    _range_pore(getParam<std::vector<unsigned int>>("range_pore")),
    _use_loaded_microstructure(getParam<bool>("use_loaded_microstructure")),
    _loaded_microstructure(isCoupled("loaded_microstructure") ? &coupledValue("loaded_microstructure") : nullptr),
    _pore_limit(getParam<Real>("pore_limit")),
    _bulk_limit(getParam<Real>("bulk_limit"))
{
    _csv_total_density = readCSV(_csv_density);
    _csv_total_density_pore = readCSV(_csv_density_pore);
}

Real
PolycrystalDensityAux::computeValue()
{
    std::vector<Real> first_col = _csv_total_density[0];
    std::vector<Real> first_col_pore = _csv_total_density_pore[0];

    //evaluate denisty_i before reading data
    //we need to create the cases

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

    Real pred_density;
    if (is_pore){
      const int call_density = _density_i[_qp] - 100;
      pred_density = first_col_pore.at(call_density);
    }
    else if (is_bulk){
      pred_density = _bulk_RDX_density;
    }
    else {
      pred_density = first_col.at(_density_i[_qp]);
    }
    
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