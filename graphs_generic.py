import sys
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from PyQt5.QtCore import Qt

class genericGraph(QDialog):
    def __init__(self, parent=None):
        super(genericGraph, self).__init__(parent)
        self.data = None
        self.parameter = None
        # a figure instance to plot on
        self.figure = plt.figure()
  
        # this is the Canvas Widget that 
        # displays the 'figure'it takes the
        # 'figure' instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)
  
        # this is the Navigation widget
        self.toolbar = NavigationToolbar(self.canvas, self)
  
        # creating a Vertical Box layout
        layout = QVBoxLayout()
          
        # adding tool bar to the layout
        layout.addWidget(self.toolbar)
          
        # adding canvas to the layout
        layout.addWidget(self.canvas)
          
        # setting layout to the main window
        self.plot_initial()
    def plot_initial(self):
        ax = self.figure.add_subplot(111)
        ax.clear()
        ax.plot([1, 2, 3], [4, 5, 6], label="Initial Plot")
        ax.set_title("Initial Plot")
        ax.set_xlabel("X-axis")
        ax.set_ylabel("Y-axis")
        ax.legend()
        self.canvas.draw()
    def update(self, data, parameter):
        ax = self.figure.add_subplot(111)
        ax.clear()
        ax.plot([1, 2, 3, 4,5,6,7,8,9,10], data, label=new_label)
        ax.set_title(title)
        ax.set_xlabel(x_units)
        ax.set_ylabel(y_unit)
        ax.legend()
        ax.grid(True)
        self.canvas.draw()
        if parameter == 0:
            x_units = 'Time (ps)'
            y_units = 'Radius of Gyration (Å)'
            title = 'Radius of Gyration Over Time'
            new_label = 'Radius of Gyration'
        if parameter == 1:
            new_label='Center of Geometry'
            x_units = 'Time (ps)'
            y_units = 'Center of Geometry (Å)'
            title = 'Center of Geometry of Protein Over Time'
        if parameter == 2:
            new_label='Total Mass'
            x_units = 'Time (ps)'
            y_units = 'Total Mass (Å)'
            new_title= 'Total Mass of Protein Over Time'
        if parameter == 3:
            new_label='Center of Mass'
            x_units='Time (ps)'
            y_units='Center of Mass (Å)'
            title ='Canter of Mass of Protein Over Time'
        if parameter == 4:
            new_label='Total Charge'
            x_units = 'Time (ps)'
            y_units ='Total Charge (Å)'
            title='Total Charge of Protein Over Time'
        if parameter == 4:
            new_label='Total Charge'
            x_units = 'Time (ps)'
            y_units ='Total Charge (Å)'
            title='Total Charge of Protein Over Time'
        if parameter == 4:
            new_label='Total Charge'
            x_units = 'Time (ps)'
            y_units ='Total Charge (Å)'
            title='Total Charge of Protein Over Time'
        if parameter == 4:
            new_label='Total Charge'
            x_units = 'Time (ps)'
            y_units ='Total Charge (Å)'
            title='Total Charge of Protein Over Time'
        
  
