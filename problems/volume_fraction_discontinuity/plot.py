#!/usr/bin/python

from file_utilities import readCSVFile
from PlotterLine import PlotterLine

solution = readCSVFile("output.csv")

u_plotter = PlotterLine("Position", "Velocity [m/s]")
u_plotter.addSet(solution["x"], solution["u_mix"], "$u_{mix}$", color=0, linetype=1)
u_plotter.setYRange(0.0, 160.0)
u_plotter.setLegendLocation("upper left")
u_plotter.save("velocity.pdf")

p_plotter = PlotterLine("Position", "Pressure [MPa]")
p_plotter.addSet(solution["x"], solution["p_mix"], "$p_{mix}$", color=0, linetype=1, scale=1e-6)
p_plotter.setYRange(0.0, 200.0)
p_plotter.save("pressure.pdf")

a_plotter = PlotterLine("Position", "Volume Fraction")
a_plotter.addSet(solution["x"], solution["a_liq"], "$\\alpha_\\ell$, Computed", color=4, linetype=2)
a_plotter.addSet(solution["x"], solution["a_vap"], "$\\alpha_v$, Computed", color=1, linetype=2)
a_plotter.setLegendLocation("center left")
a_plotter.save("volume_fraction.pdf")
