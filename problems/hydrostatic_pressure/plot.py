#!/usr/bin/python

from file_utilities import readCSVFile
from PlotterLine import PlotterLine

solution = readCSVFile("output.csv")

u_plotter = PlotterLine("Position", "Velocity [m/s]")
u_plotter.addSet(solution["x"], solution["u"], "$u$")
u_plotter.save("velocity.png")

p_plotter = PlotterLine("Position", "Pressure [MPa]")
p_plotter.addSet(solution["x"], solution["p"], "$p$", scale=1e-6)
p_plotter.save("pressure.png")

r_plotter = PlotterLine("Position", "Density [kg/m$^3$]")
r_plotter.addSet(solution["x"], solution["r"], "$\\rho$")
r_plotter.save("density.png")
