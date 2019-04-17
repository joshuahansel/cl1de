#ifndef RunParameters_H
#define RunParameters_H

#include <map>
#include <memory>
#include <string>

class TimeStepSizer;

class RunParameters
{
public:
  RunParameters(int argc, char* argv[]);

  std::string getStringParameter(const std::string & name);
  unsigned int getIntParameter(const std::string & name);
  double getFloatParameter(const std::string & name);

  bool parameterExists(const std::string & name) const;
  void checkParameterExists(const std::string & name) const;

  const unsigned int & getNumberOfStages() const;
  const unsigned int & getNumberOfElements() const;
  const TimeStepSizer & getTimeStepSizer() const;
  const double & getEndTime() const;

  bool specifiedEndTime() const {return _t_set;}

protected:
  void parseInputFile(int argc, char* argv[]);
  void checkAllParametersUsed() const;

  std::string _input_file_name;
  std::map<std::string, std::string> _param_map;
  std::map<std::string, bool> _param_used_map;

  unsigned int _n_stages;
  unsigned int _n_elem;
  std::shared_ptr<TimeStepSizer> _time_step_sizer;
  double _t;

  bool _t_set;
};

#endif /* RunParameters_H */
