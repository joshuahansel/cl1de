#!/usr/bin/python

from file_utilities import readCSVFile
from PlotterLine import PlotterLine

solution = readCSVFile("output.csv")

u_plotter = PlotterLine("Position", "Velocity [m/s]")
u_plotter.addSet(solution["x"], solution["u_liq"], "$u_\\ell$, Computed", color=4, linetype=2)
u_plotter.addSet(solution["x"], solution["u_vap"], "$u_v$, Computed", color=1, linetype=2)
u_plotter.save("velocity.png")

p_plotter = PlotterLine("Position", "Pressure [MPa]")
p_plotter.addSet(solution["x"], solution["p_liq"], "$p_\\ell$, Computed", color=4, linetype=2, scale=1e-6)
p_plotter.addSet(solution["x"], solution["p_vap"], "$p_v$, Computed", color=1, linetype=2, scale=1e-6)
p_plotter.save("pressure.png")

r_liq_plotter = PlotterLine("Position", "Density [kg/m$^3$]")
r_liq_plotter.addSet(solution["x"], solution["r_liq"], "$\\rho_\\ell$, Computed", color=4, linetype=2)
r_liq_plotter.save("density_liquid.png")

r_vap_plotter = PlotterLine("Position", "Density [kg/m$^3$]")
r_vap_plotter.addSet(solution["x"], solution["r_vap"], "$\\rho_v$, Computed", color=1, linetype=2)
r_vap_plotter.save("density_vapor.png")

a_plotter = PlotterLine("Position", "Volume Fraction")
a_plotter.addSet(solution["x"], solution["a_liq"], "$\\alpha_\\ell$, Computed", color=4, linetype=2)
a_plotter.addSet(solution["x"], solution["a_vap"], "$\\alpha_v$, Computed", color=1, linetype=2)
a_plotter.setYRange(0.0, 1.0)
a_plotter.save("volume_fraction.png")
