//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "mistTestApp.h"
#include "mistApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
mistTestApp::validParams()
{
  InputParameters params = mistApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

mistTestApp::mistTestApp(const InputParameters & parameters) : MooseApp(parameters)
{
  mistTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

mistTestApp::~mistTestApp() {}

void
mistTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  mistApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"mistTestApp"});
    Registry::registerActionsTo(af, {"mistTestApp"});
  }
}

void
mistTestApp::registerApps()
{
  registerApp(mistApp);
  registerApp(mistTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
mistTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  mistTestApp::registerAll(f, af, s);
}
extern "C" void
mistTestApp__registerApps()
{
  mistTestApp::registerApps();
}
