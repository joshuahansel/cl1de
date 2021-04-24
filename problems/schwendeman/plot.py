#!/usr/bin/python

from file_utilities import readCSVFile
from PlotterLine import PlotterLine

solution = readCSVFile("output.csv")

r_liq_plotter = PlotterLine("Position", "Density [kg/m$^3$]")
r_liq_plotter.addSet(solution["x"], solution["r_liq"], "$\\rho_\ell$", color=4, linetype=2)
r_liq_plotter.save("density_liquid.png")

r_vap_plotter = PlotterLine("Position", "Density [kg/m$^3$]")
r_vap_plotter.addSet(solution["x"], solution["r_vap"], "$\\rho_v$", color=1, linetype=2)
r_vap_plotter.save("density_vapor.png")

p_plotter = PlotterLine("$\\alpha_\ell$", "Pressure [Pa]")
p_plotter.addSet(solution["a_liq"], solution["p_liq"], "$p_\\ell$", color=4, linetype=0, marker=1)
p_plotter.save("pressure.png")

a_plotter = PlotterLine("Position", "Volume Fraction")
a_plotter.addSet(solution["x"], solution["a_liq"], "$\\alpha_\ell$", color=4, linetype=2)
a_plotter.save("volume_fraction.png")
