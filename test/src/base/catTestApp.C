//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "catTestApp.h"
#include "catApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
catTestApp::validParams()
{
  InputParameters params = catApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

catTestApp::catTestApp(InputParameters parameters) : MooseApp(parameters)
{
  catTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

catTestApp::~catTestApp() {}

void
catTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  catApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"catTestApp"});
    Registry::registerActionsTo(af, {"catTestApp"});
  }
}

void
catTestApp::registerApps()
{
  registerApp(catApp);
  registerApp(catTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
catTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  catTestApp::registerAll(f, af, s);
}
extern "C" void
catTestApp__registerApps()
{
  catTestApp::registerApps();
}
