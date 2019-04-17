#include "RunParameters.h"
#include "TimeStepSizerConstant.h"
#include "TimeStepSizerCFL.h"
#include "utils.h"

#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>

RunParameters::RunParameters(int argc, char* argv[])
  : _t_set(false)
{
  parseInputFile(argc, argv);

  // time step integrator
  _n_stages = getIntParameter("temporal_order");

  // time step sizer
  const std::string dt_option = getStringParameter("dt_option");
  if (dt_option == "constant")
  {
    const double dt = getFloatParameter("dt");
    _time_step_sizer = std::make_shared<TimeStepSizerConstant>(dt);
  }
  else if (dt_option == "cfl")
  {
    const double cfl = getFloatParameter("cfl");
    _time_step_sizer = std::make_shared<TimeStepSizerCFL>(cfl);
  }
  else
    throwInvalidStringParameterValueError("dt_option", dt_option);

  // number of elements
  _n_elem = getIntParameter("n_elem");

  // end time
  if (parameterExists("t"))
  {
    _t = getFloatParameter("t");
    _t_set = true;
  }
}

void RunParameters::parseInputFile(int argc, char* argv[])
{
  if (argc == 2)
    _input_file_name = argv[1];
  else
    throwError("USAGE: " + std::string(argv[0]) + " <input_file_name>");

  std::ifstream input_file(_input_file_name);
  if (!input_file.good())
    throwError("The input file '" + _input_file_name + "' does not exist.");

  std::string line;
  unsigned int line_no = 1;
  while (std::getline(input_file, line))
  {
    std::vector<std::string> elems;
    std::stringstream ss(line);
    std::string item;
    while (std::getline(ss, item, ' '))
      *(std::back_inserter(elems)++) = item;
    if (elems.size() == 2)
    {
      _param_map[elems[0]] = elems[1];
      _param_used_map[elems[0]] = false;
    }
    else
      throwError(_input_file_name + ", line " + std::to_string(line_no) + ": Invalid input line.");

    line_no++;
  }
  input_file.close();
}

void RunParameters::checkAllParametersUsed() const
{
  for (auto const & it: _param_used_map)
    if (!it.second)
      throwError("The parameter '" + it.first + "' was not used.");
}

std::string RunParameters::getStringParameter(const std::string & name)
{
  checkParameterExists(name);
  _param_used_map[name] = true;
  return _param_map.at(name);
}

unsigned int RunParameters::getIntParameter(const std::string & name)
{
  checkParameterExists(name);
  _param_used_map[name] = true;
  return std::stoi(_param_map.at(name));
}

double RunParameters::getFloatParameter(const std::string & name)
{
  checkParameterExists(name);
  _param_used_map[name] = true;
  return std::stod(_param_map.at(name));
}

bool RunParameters::parameterExists(const std::string & name) const
{
  auto it = _param_map.find(name);
  return it != _param_map.end();
}

void RunParameters::checkParameterExists(const std::string & name) const
{
  if (!parameterExists(name))
    throwError(_input_file_name + ": The parameter '" + name + "' was not found.");
}

const unsigned int & RunParameters::getNumberOfStages() const
{
  return _n_stages;
}

const unsigned int & RunParameters::getNumberOfElements() const
{
  return _n_elem;
}

const TimeStepSizer & RunParameters::getTimeStepSizer() const
{
  return *_time_step_sizer;
}

const double & RunParameters::getEndTime() const
{
  if (_t_set)
    return _t;
  else
    throwError("The end time was not set.", __PRETTY_FUNCTION__);
}
