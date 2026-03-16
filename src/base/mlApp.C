#include "mlApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
mlApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

mlApp::mlApp(const InputParameters & parameters) : MooseApp(parameters)
{
  mlApp::registerAll(_factory, _action_factory, _syntax);
}

mlApp::~mlApp() {}

void
mlApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<mlApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"mlApp"});
  Registry::registerActionsTo(af, {"mlApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
mlApp::registerApps()
{
  registerApp(mlApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
mlApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  mlApp::registerAll(f, af, s);
}
extern "C" void
mlApp__registerApps()
{
  mlApp::registerApps();
}
