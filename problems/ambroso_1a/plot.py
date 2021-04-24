from file_utilities import readCSVFile
from PlotterLine import PlotterLine

solution = readCSVFile("output.csv")

a_liq_plotter = PlotterLine("Position", "Volume Fraction")
a_liq_plotter.addSet(solution["x"], solution["a_liq"], "$\\alpha_\ell$", color=4, linetype=2)
a_liq_plotter.save("volume_fraction.png")

r_liq_plotter = PlotterLine("Position", "Density [kg/m$^3$]")
r_liq_plotter.addSet(solution["x"], solution["r_liq"], "$\\rho_\ell$", color=4, linetype=2)
r_liq_plotter.save("density_liquid.png")

r_vap_plotter = PlotterLine("Position", "Density [kg/m$^3$]")
r_vap_plotter.addSet(solution["x"], solution["r_vap"], "$\\rho_v$", color=1, linetype=2)
r_vap_plotter.save("density_vapor.png")

u_liq_plotter = PlotterLine("Position", "Velocity [m/s]")
u_liq_plotter.addSet(solution["x"], solution["u_liq"], "$u_\ell$", color=4, linetype=2)
u_liq_plotter.save("velocity_liquid.png")

u_vap_plotter = PlotterLine("Position", "Velocity [m/s]")
u_vap_plotter.addSet(solution["x"], solution["u_vap"], "$u_v$", color=1, linetype=2)
u_vap_plotter.save("velocity_vapor.png")

p_liq_plotter = PlotterLine("Position", "Pressure [m/s]")
p_liq_plotter.addSet(solution["x"], solution["p_liq"], "$p_\ell$", color=4, linetype=2)
p_liq_plotter.save("pressure_liquid.png")

p_vap_plotter = PlotterLine("Position", "Pressure [m/s]")
p_vap_plotter.addSet(solution["x"], solution["p_vap"], "$p_v$", color=1, linetype=2)
p_vap_plotter.save("pressure_vapor.png")
