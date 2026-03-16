//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "mlTestApp.h"
#include "mlApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
mlTestApp::validParams()
{
  InputParameters params = mlApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

mlTestApp::mlTestApp(const InputParameters & parameters) : MooseApp(parameters)
{
  mlTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

mlTestApp::~mlTestApp() {}

void
mlTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  mlApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"mlTestApp"});
    Registry::registerActionsTo(af, {"mlTestApp"});
  }
}

void
mlTestApp::registerApps()
{
  registerApp(mlApp);
  registerApp(mlTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
mlTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  mlTestApp::registerAll(f, af, s);
}
extern "C" void
mlTestApp__registerApps()
{
  mlTestApp::registerApps();
}
