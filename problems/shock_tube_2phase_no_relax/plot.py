#!/usr/bin/python

from file_utilities import readCSVFile
from PlotterLine import PlotterLine

solution = readCSVFile("output.csv")

u_plotter = PlotterLine("Position", "Velocity [m/s]")
u_plotter.addSet(solution["x"], solution["u_liq"], "$u_\\ell$, Computed", color=1, linetype=2)
u_plotter.addSet(solution["x"], solution["u_vap"], "$u_v$, Computed", color=4, linetype=2)
u_plotter.save("velocity.pdf")

p_plotter = PlotterLine("Position", "Pressure [MPa]")
p_plotter.addSet(solution["x"], solution["p_liq"], "$p_\\ell$, Computed", color=1, linetype=2, scale=1e-6)
p_plotter.addSet(solution["x"], solution["p_vap"], "$p_v$, Computed", color=4, linetype=2, scale=1e-6)
p_plotter.save("pressure.pdf")
