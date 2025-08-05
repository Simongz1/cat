#include "ComputeYRates.h"
#include "RankTwoTensor.h"

registerMooseObject("catApp", ComputeYRates);

InputParameters
ComputeYRates::validParams()
{
    InputParameters params = Material::validParams();
    params.addClassDescription("compute the rate of reaction of each species");
    params.addRequiredCoupledVar("temperature", "temperature");
    params.addRequiredCoupledVar("Y1", "mass fraction 1");
    params.addRequiredCoupledVar("Y2", "mass fraction 2");
    params.addRequiredCoupledVar("Y3", "mass fraction 3");
    params.addRequiredCoupledVar("Y4", "mass fraction 4");
    params.addRequiredParam<std::vector<Real>>("Z", "list of pre-exponential factors");
    params.addRequiredParam<std::vector<Real>>("E", "list of energies of activation");
    params.addRequiredParam<Real>("Rg", "thermodynamic constan");
    params.addRequiredParam<std::vector<std::string>>("Ydot_names", "name of rates to be declared");
    return params;
}

ComputeYRates::ComputeYRates(const InputParameters & parameters)
  : Material(parameters),
    _T(coupledValue("temperature")),
    _Y1(coupledValue("Y1")),
    _Y2(coupledValue("Y2")),
    _Y3(coupledValue("Y3")),
    _Y4(coupledValue("Y4")),
    _Z(getParam<std::vector<Real>>("Z")),
    _E(getParam<std::vector<Real>>("E")),
    _Rg(getParam<Real>("Rg")),
    _Ydot_names(getParam<std::vector<std::string>>("Ydot_names")), //get a vector that has the names of the material properties to create
    _Ydot(declareProperty<std::vector<Real>>("Ydot")),
    _arrate(declareProperty<std::vector<Real>>("arrate"))

{   
    //create a field vector with the coupled variables pointers
    int ProblemSize = _Ydot_names.size();
    _Y = {&_Y1, &_Y2, &_Y3, &_Y4};
    _Ydot.resize(ProblemSize);
    _arrate.resize(ProblemSize);

    //allocate memory to avoid bottlenecking when appending new elements to the vector
    /////_Ydot.reserve(_Ydot_names.size()); //saves enough memory to at most allocate n entries
    /////_arrate.reserve(_Ydot_names.size());

    //iterate to obtain each material property and declare it
    //for (unsigned int i = 0; i < _Ydot_names.size(); i++)
    //{
        //Ydot is created as a vector of material property pointer values
        /////_Ydot.push_back(&declareProperty<Real>(_Ydot_names[i]));
        //compute the pre-exponential rates for each reaction   
        /////_arrate.push_back(&declareProperty<Real>("arrate" + std::to_string(i + 1)));
    //}
}

void
ComputeYRates::computeQpProperties()
{
    //define the reaction rate terms and reaction evolution
    //for (unsigned int i = 0; i < _arrate.size(); i++)
    //{
    //    (*_arrate[i])[_qp] = _Z[i] * std::exp(- _E[i] / (_Rg * _T[_qp]));
    //}

    //(*_Ydot[0])[_qp] = - (*_arrate[0])[_qp] * (*_Y[0])[_qp];
    //(*_Ydot[1])[_qp] = (*_arrate[0])[_qp] * (*_Y[0])[_qp] - (*_arrate[1])[_qp] * (*_Y[1])[_qp];
    //(*_Ydot[2])[_qp] = (*_arrate[1])[_qp] * (*_Y[1])[_qp] - (*_arrate[2])[_qp] * std::pow((*_Y[2])[_qp], 2.0);
    //(*_Ydot[3])[_qp] = (*_arrate[2])[_qp] * std::pow((*_Y[2])[_qp], 2.0);
    for (unsigned int i; i < _Ydot_names.size(); i++)
    {
        _arrate[_qp][i] = _Z[i] * std::exp(- _E[i] / (_Rg * _T[_qp]));
    }

    _Ydot[_qp][0] = - _arrate[_qp][0] * (*_Y[0])[_qp];
    _Ydot[_qp][1] = _arrate[_qp][0] * (*_Y[0])[_qp] - _arrate[_qp][1] * (*_Y[1])[_qp];
    _Ydot[_qp][2] = _arrate[_qp][1] * (*_Y[1])[_qp] - _arrate[_qp][2] * std::pow((*_Y[2])[_qp], 2.0);
    _Ydot[_qp][3] = _arrate[_qp][2] * std::pow((*_Y[2])[_qp], 2.0);

}