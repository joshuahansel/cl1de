#!/usr/bin/python

import numpy as np

from file_utilities import readCSVFile
from PlotterLine import PlotterLine

exact = readCSVFile("exact.csv")
solution = readCSVFile("output.csv")

r_plotter = PlotterLine("Position", "Density")
r_plotter.addSet(exact["x"], exact["r"], "$\\rho_\\ell$, Exact", color=4, linetype=1)
r_plotter.addSet(exact["x"], np.flip(exact["r"]), "$\\rho_v$, Exact", color=1, linetype=1)
r_plotter.addSet(solution["x"], solution["r_liq"], "$\\rho_\\ell$, Computed", color=4, linetype=2)
r_plotter.addSet(solution["x"], solution["r_vap"], "$\\rho_v$, Computed", color=1, linetype=2)
r_plotter.setLegendLocation("right")
r_plotter.save("density.pdf")

u_plotter = PlotterLine("Position", "Velocity")
u_plotter.addSet(exact["x"], exact["u"], "$u_\\ell$, Exact", color=4, linetype=1)
u_plotter.addSet(exact["x"], -np.flip(exact["u"]), "$u_v$, Exact", color=1, linetype=1)
u_plotter.addSet(solution["x"], solution["u_liq"], "$u_\\ell$, Computed", color=4, linetype=2)
u_plotter.addSet(solution["x"], solution["u_vap"], "$u_v$, Computed", color=1, linetype=2)
u_plotter.setLegendLocation("lower right")
u_plotter.save("velocity.pdf")

p_plotter = PlotterLine("Position", "Pressure")
p_plotter.addSet(exact["x"], exact["p"], "$p_\\ell$, Exact", color=4, linetype=1)
p_plotter.addSet(exact["x"], np.flip(exact["p"]), "$p_v$, Exact", color=1, linetype=1)
p_plotter.addSet(solution["x"], solution["p_liq"], "$p_\\ell$, Computed", color=4, linetype=2)
p_plotter.addSet(solution["x"], solution["p_vap"], "$p_v$, Computed", color=1, linetype=2)
p_plotter.setLegendLocation("right")
p_plotter.save("pressure.pdf")
