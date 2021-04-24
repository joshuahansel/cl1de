from file_utilities import readCSVFile
from PlotterLine import PlotterLine

def plotDataSets(data_sets, data_names):
  for var in data_names:
    desc, symbol = data_names[var]
    plotter = PlotterLine("$x$", desc + ", $" + symbol + "$")
    for i, data_set in enumerate(data_sets):
      set_name, data = data_set
      plotter.addSet(data["x"], data[var], set_name, color=i)
    plotter.save(var + ".png")

euler1phase_data_names = {
  "r": ("Density", "\\rho"),
  "u": ("Velocity", "u"),
  "p": ("Pressure", "p")
  }

def plotEuler1PhaseDataSets(data_sets):
  plotDataSets(data_sets, euler1phase_data_names)
