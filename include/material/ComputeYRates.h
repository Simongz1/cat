#include "Material.h"
#include "RankTwoTensor.h"
#include <vector>
#include <cmath>

class ComputeYRates;

class ComputeYRates : public Material
{
public:
    ComputeYRates(const InputParameters & parameters);
    static InputParameters validParams();

protected:
    virtual void computeQpProperties() override;

private:
    const VariableValue & _T;
    const VariableValue & _Y1;
    const VariableValue & _Y2;
    const VariableValue & _Y3;
    const VariableValue & _Y4;
    const std::vector<Real> _Z;
    const std::vector<Real> _E;
    const Real _Rg;
    const std::vector<std::string> _Ydot_names;
    //create vectors to store variable and material values
    std::vector<const VariableValue *> _Y; //pointer value
    MaterialProperty<std::vector<Real>> & _Ydot;
    //std::vector<MaterialProperty<Real> *> _arrate;
    MaterialProperty<std::vector<Real>> & _arrate;
};