#include "CSVOutput.h"
#include "utils.h"

CSVOutput::CSVOutput(const std::string & name) : _name(name)
{
  _outfile.open(_name);
}

void CSVOutput::addOutput(const std::vector<double> & values, const std::string & name)
{
  _data_names.push_back(name);
  _data_vectors.push_back(&values);
}

void CSVOutput::save()
{
  const unsigned int n_outputs = _data_names.size();

  if (n_outputs == 0)
    throwError("There must be at least one vector added.", __PRETTY_FUNCTION__);

  unsigned int n = _data_vectors[0]->size();

  for (unsigned int i = 0; i < n_outputs; i++)
    if (_data_vectors[i]->size() != n)
      throwError("All vectors must have the same length.", __PRETTY_FUNCTION__);

  // output header row
  _outfile << _data_names[0];
  for (unsigned int k = 1; k < n_outputs; k++)
    _outfile << "," << _data_names[k];
  _outfile << "\n";

  // loop over number of points
  for (unsigned int i = 0; i < n; i++)
  {
    _outfile << (*_data_vectors[0])[i];
    // loop over outputs
    for (unsigned int k = 1; k < n_outputs; k++)
      _outfile << "," << (*_data_vectors[k])[i];
    _outfile << "\n";
  }

  _outfile.close();
}
