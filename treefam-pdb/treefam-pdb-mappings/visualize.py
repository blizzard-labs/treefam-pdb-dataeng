import random as rnd
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from matplotlib.text import Annotation
import seaborn as sns

#! WORK IN PROGRESS
#* Right now just some code from a sample tutorial, creating a simple visualizer for the final data

sns.set()
labels_color_map = {0: '#20b2aa', 1: '#ff7373'}
no_examples = 50
generated_data = [(x, rnd.randint(0, no_examples)) for x in range(0, no_examples)]
generated_labels = ["Label for instance #{0}".format(i) for i in range(0, no_examples)]

instances_colors = []
axis_values_x = []
axis_values_y = []
for index, instance in enumerate(generated_data):
    coordinate_x, coordinate_y = instance
    color = labels_color_map[index % 2]

    instances_colors.append(color)
    axis_values_x.append(coordinate_x)
    axis_values_y.append(coordinate_y)

fig = plt.figure(figsize=(20, 16))
ax = plt.subplot()

def draw_scatterplot():
    ax.scatter(
        axis_values_x,
        axis_values_y,
        c=instances_colors,
        picker=True
    )

draw_scatterplot()

def annotate(axis, text, x, y):
    text_annotation = Annotation(text, xy=(x, y), xycoords='data')
    axis.add_artist(text_annotation)


def onpick(event):
    ind = event.ind
    label_pos_x = event.mouseevent.xdata
    label_pos_y = event.mouseevent.ydata
    offset = 0

    for i in ind:
        label = generated_labels[i]
        annotate(
            ax,
            label,
            label_pos_x + offset,
            label_pos_y + offset
        )
        ax.figure.canvas.draw_idle()
        offset += 0.01

fig.canvas.mpl_connect('pick_event', onpick)
ax_clear_all = plt.axes([0.0, 0.0, 0.1, 0.05])
button_clear_all = Button(ax_clear_all, 'Clear all')


def onclick(event):
    ax.cla()
    draw_scatterplot()
    ax.figure.canvas.draw_idle()

button_clear_all.on_clicked(onclick)
plt.plot()
plt.show()