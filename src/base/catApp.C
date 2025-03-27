#include "catApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
catApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

catApp::catApp(InputParameters parameters) : MooseApp(parameters)
{
  catApp::registerAll(_factory, _action_factory, _syntax);
}

catApp::~catApp() {}

void 
catApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<catApp>(f, af, s);
  Registry::registerObjectsTo(f, {"catApp"});
  Registry::registerActionsTo(af, {"catApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
catApp::registerApps()
{
  registerApp(catApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
catApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  catApp::registerAll(f, af, s);
}
extern "C" void
catApp__registerApps()
{
  catApp::registerApps();
}
