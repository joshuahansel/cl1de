#!/usr/bin/python

from file_utilities import readCSVFile
from PlotterLine import PlotterLine

solution = readCSVFile("output.csv")

r_plotter = PlotterLine("Position", "Density [kg/m$^3$]")
r_plotter.addSet(solution["x"], solution["r_mix"], "$\\rho_{mix}$", color=1, linetype=1)
r_plotter.save("density.png")

u_plotter = PlotterLine("Position", "Velocity [m/s]")
u_plotter.addSet(solution["x"], solution["u_liq"], "$u_\\ell$", color=4, linetype=2)
u_plotter.addSet(solution["x"], solution["u_vap"], "$u_v$", color=1, linetype=2)
u_plotter.save("velocity.png")

p_plotter = PlotterLine("Position", "Pressure [MPa]")
p_plotter.addSet(solution["x"], solution["p_liq"], "$p_\\ell$", color=4, linetype=2, scale=1e-6)
p_plotter.addSet(solution["x"], solution["p_vap"], "$p_v$", color=1, linetype=2, scale=1e-6)
p_plotter.save("pressure.png")

a_plotter = PlotterLine("Position", "Volume Fraction")
a_plotter.addSet(solution["x"], solution["a_vap"], "$\\alpha_v$", color=1, linetype=2)
a_plotter.save("volume_fraction.png")
