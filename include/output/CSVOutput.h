#ifndef CSVOutput_H
#define CSVOutput_H

#include <string>
#include <vector>
#include <fstream>

class CSVOutput
{
public:
  CSVOutput(const std::string & name);

  void addOutput(const std::vector<double> & values, const std::string & name);
  void save();

protected:
  const std::string _name;

  std::vector<std::string> _data_names;
  std::vector<std::vector<double> const *> _data_vectors;

  std::ofstream _outfile;
};

#endif /* CSVOutput_H */
