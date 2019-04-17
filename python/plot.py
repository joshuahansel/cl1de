# Plots variable from the output file 'output.csv' in the current directory.
#
# USAGE: python plot.py <var>
#   <var>  Variable to plot
#

import argparse

from file_utilities import readCSVFile
from error_utilities import raiseException
from PlotterLine import PlotterLine

parser = argparse.ArgumentParser(description='Plots solutions')
parser.add_argument('variable', help='Variable to plot')
parsed = parser.parse_args()
variable = parsed.variable

data = readCSVFile("output.csv")

data_names = {
  "A": ("Area", "A"),
  "r": ("Density", "\\rho"),
  "u": ("Velocity", "u"),
  "p": ("Pressure", "p"),
  "T": ("Temperature", "T"),
  }

if variable in data_names:
  desc, symbol = data_names[variable]
else:
  raiseException("Variable '" + variable + "' is invalid.")

plotter = PlotterLine("$x$", desc + ", $" + symbol + "$")
plotter.addSet(data["x"], data[variable], "$" + symbol + "$")
plotter.show()
