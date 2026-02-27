#include "mistApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
mistApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

mistApp::mistApp(const InputParameters & parameters) : MooseApp(parameters)
{
  mistApp::registerAll(_factory, _action_factory, _syntax);
}

mistApp::~mistApp() {}

void
mistApp::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<mistApp>(f, af, syntax);
  Registry::registerObjectsTo(f, {"mistApp"});
  Registry::registerActionsTo(af, {"mistApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
mistApp::registerApps()
{
  registerApp(mistApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
mistApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  mistApp::registerAll(f, af, s);
}
extern "C" void
mistApp__registerApps()
{
  mistApp::registerApps();
}
