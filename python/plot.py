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
  "a_liq": ("Liquid Volume Fraction", "\\alpha_\\ell"),
  "a_vap": ("Vapor Volume Fraction", "\\alpha_v"),
  "r_liq": ("Liquid Density", "\\rho_\\ell"),
  "r_vap": ("Vapor Density", "\\rho_v"),
  "u_liq": ("Liquid Velocity", "u_\\ell"),
  "u_vap": ("Vapor Velocity", "u_v"),
  "p_liq": ("Liquid Pressure", "p_\\ell"),
  "p_vap": ("Vapor Pressure", "p_v"),
  "T_liq": ("Liquid Temperature", "T_\\ell"),
  "T_vap": ("Vapor Temperature", "T_v")
  }

if variable in data_names:
  desc, symbol = data_names[variable]
else:
  raiseException("Variable '" + variable + "' is invalid.")

plotter = PlotterLine("$x$", desc + ", $" + symbol + "$")
plotter.addSet(data["x"], data[variable], "$" + symbol + "$")
plotter.show()
