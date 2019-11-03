# Matplotlib :
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.ion()
# Numpy :
import numpy as np
class Cameleon:
    """
    This is a class to create a matplotlib graphic using a non-pyplot format.
    It consist in a simplification of the normally hard implementation of the
    figure into a tkinter frame.

    Attributes:
        parent : tkinter frame in which the graphic is placed in.
        Fig : Figure object of the matplotlib class analogue to pyplot.figure
        axes : Matplotlib Axis object created to plot data
        Line : Matplotlib Line object created to update data in the given axis
        canvas : Matplotlib object that generates the 'drawing' in a tkinter
        frame it is mainly used to update the graphic in real time.
        toolbar : Matplotlib toolbar normally under any pyplot graph.
    """
    def __init__(self, parent, axis_name=['', ''], figsize=[1, 1]):
        """
        The constructor for the GraphicFrame Class.

        Parameters:
            parent : tkinter Frame object where the object is placed in.
            axis_name : This is a list of two strings that will be respectivly
            the x and y axis name. (Should be Latex friendly)
            figsize : This is the initial figure size (The figure size is
            automaticly updated when the window is changed in size)
        """
        self.parent = parent
        self.Fig = Figure(dpi=100, figsize=figsize)
        # Adjusting the axis in the figure and setting up necessary parameters
        # for naming the axis
        self.axes = self.Fig.add_axes([0.11, 0.05, 0.85, 0.88])
        self.axes.set_aspect('auto', adjustable='box')
        self.axes.set_adjustable('box')
        self.Line, = self.axes.plot([], [])
        self.axes.tick_params(axis='both', which='major', labelsize=8)
        self.axes.grid()
        #Creating toolbar for convinience reason and canvas to host the figure
        self.canvas = FigureCanvasTkAgg(self.Fig, parent)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(expand=True, fill='both')
#        self.canvas._tkcanvas.pack()

    def change_dimensions(self, event):
        """
        This function is a way to update the size of the figure when you change
        the size of your window automaticly it takes the width of your parent
        and the dpi of your figure to update the height and width.

        How to set it up : Your_Frame.bind('<Configure>', Your_Graph.change_dimensions)

        Parameters:
            event: An event is an object in tkinter that is created when you
            click on the screen. See the given documentation for the specifics.
            This parameter automaticly sent through when you click on the
            line.
        """
        width = event.width/self.Fig.get_dpi()
        height = event.height/self.Fig.get_dpi()
        self.Fig.set_size_inches(w=width, h=height)

    def update_graph(self):
        """
        This function is a compilation of two line to update the figure canvas
        so it update the values displayed whitout recreating the figure in the
        tkinter frame.

        """
        self.Fig.canvas.draw()
        self.Fig.canvas.flush_events()

    def destroy_graph(self):
        """
        This function is a compilation of two line to destroy a graph ie if
        you want to replace it or just get rid of it. It does not destroy the
        class itself so creating the canvas and the toolbar would make it
        appear the same way it was before.

        """
        self.canvas.get_tk_widget().destroy()
        self.toolbar.destroy()

    def log_scale(self):
        """
        This function is changing the y axis to make it a logarithmic scale.
        """
        self.axes.set_yscale('log')
        self.update_graph()

    def lin_scale(self):
        """
        This function is changing/reverting the y axis back to a linear scale.
        """
        self.axes.set_yscale('linear')
        self.update_graph()

class Cameleon2:
    """
    This is a class to create an vertical line that can be drag in
    a vertical fashion. Documentation of this part is not torough as it does
    not come from me. I took some part of another code and adapted it to this
    one.

    Attributes:
        parent : tkinter frame in which the matplotlib pyplot is gridded into.
        press : This attribute is a state of the draggable line it describe
        it's current state and determines if it's available to be moved around.
        background : ??
        x : This is the position on the y axis of the draggable line it is
        updated everytimes it is clicked on.
        line : This is an object from matplotlib. It consist in a straight line
        going from minus inf to inf. It inherite from the Axis class of
        matplotlib.
    """
    Lock = None

    def __init__(self, parent=None, x=0.5, y=0.5, axes=None,
                velx=0.01, vely=0.01):
        """
        The constructor for the HorizontalDraggableLine Class.

        Parameters:
            parent : tkinter Frame object where the object is placed in.
            y : Initial position of your horizontal line on the y axis
            Axis : axis of your pyplot graphic that will contain this line.
        """

        self.parent = parent
        self.press = None
        self.background = None
        self.x = x
        self.y = y
        #self.line, = axes.plot([x.get(), x.get()+velx.get()],
        #                       [y.get(), y.get()+vely.get()],
        #                       'k')
        self.line, = axes.plot([x.get(), x.get()],
                               [y.get(), y.get()],
                               '*k')
        #self.line.axes.annotate('',
        #xytext=(x.get(), y.get()),
        #xy=(x.get()+velx.get(), y.get()+vely.get()),
        #arrowprops=dict(arrowstyle="->", color='k'),
        #size=12)
        self.connect()

    def connect(self):
        """
        This function connect all the events we need to control the line
        vertically

        """
        self.cidpress1 = self.line.figure.canvas.mpl_connect('button_press_event', self.on_press)
        self.cidrelease1 = self.line.figure.canvas.mpl_connect('button_release_event', self.on_release)
        self.cidmotion1 = self.line.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def on_press(self, event):
        """
        This function activate itself when there is a click on the line

        Parameters:
            event: An event is an object in tkinter that is created when you
            click on the screen. See the given documentation for the specifics.
            This parameter automaticly sent through when you click on the
            line.
        """
        if event.inaxes != self.line.axes: return
        if self.Lock is not None: return
        contains, attrd = self.line.contains(event)
        if contains != True: return
        self.press = (self.line.get_xydata()), event.xdata, event.ydata
        self.Lock = self

        canvas = self.line.figure.canvas
        axes = self.line.axes
        self.line.set_animated(True)

        canvas.draw()
        self.background = canvas.copy_from_bbox(self.line.axes.bbox)
        axes.draw_artist(self.line)

        canvas.blit(axes.bbox)

    def on_release(self, event):
        """
        This function activate itself when you release the click on the line

        Parameters:
            event: An event is an object in tkinter that is created when you
            click on the screen. See the given documentation for the specifics.
            This parameter automaticly sent through when you click on the
            line.
        """

        if self.Lock is not self: return

        self.press = None
        self.Lock = None

        self.line.set_animated(False)
        self.BOX.set_animated(False)

        self.background = None
        self.line.figure.canvas.draw()

        self.x.set(self.line.get_xdata())
        self.y.set(self.line.get_ydata())

    def on_motion(self, event):
        """
        This function activate itself once you have clicked a line and you move
        it around.

        Parameters:
            event: An event is an object in tkinter that is created when you
            click on the screen. See the given documentation for the specifics.
            This parameter automaticly sent through when you click on the
            line.
        """
        if self.Lock is not self: return
        if event.inaxes != self.line.axes: return
        array, xpress, ypress = self.press
        x0 = array[0][0]
        dx = event.xdata - xpress
        self.line.set_xdata(x0 + dx)
        y0 = array[0][0]
        dy = event.ydata - xpress
        self.line.set_ydata(y0 + yx)

        canvas = self.line.figure.canvas
        axes = self.line.axes
        canvas.restore_region(self.background)

        axes.draw_artist(self.line)
        axes.draw_artist(self.parent.line_list[1].BOX)
        self.x.set(self.line.get_xdata())
        self.y.set(self.line.get_ydata())

        canvas.blit(axes.bbox)

    def disconnect(self):
        """
        This function if my understanding is right is to disconnect all of the
        function above from the line. Once you call this it wont be possible to
        control the line anymore.
        """
        self.cidpress1 = self.line.figure.canvas.mpl_disconnect('button_press_event', self.on_press)
        self.cidrelease1 = self.line.figure.canvas.mpl_disconnect('button_release_event', self.on_release)
        self.cidmotion1 = self.line.figure.canvas.mpl_disconnect('motion_notify_event', self.on_motion)

    def update_position(self, graph):
        self.line.set_ydata(self.y.get())
        self.line.set_xdata(self.x.get())
        #self.line.axes.annotate('',
        #xytext=(x.get(), y.get()),
        #xy=(x.get()+velx.get(), y.get()+vely.get()),
        #arrowprops=dict(arrowstyle="->", color='k'),
        #size=12)
        graph.update_graph()
