#Krishna Bhatt @ Holmes Lab (UC Berkeley) 2024

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from matplotlib.text import Annotation
import seaborn as sns

from utils import csv2numpy, numpy2csv, numpy2json, json2numpy

#Import and Load Alignment Contacts

alignment_contacts = json2numpy("mapping.json")
contact_threshold = 10.0

scatter_data = []
scatter_info = []

for i in range(len(alignment_contacts)):
    for j in range(len(alignment_contacts[i])):
        if alignment_contacts[i, j] < contact_threshold and alignment_contacts[i, j] >= 0 and i != j:
            scatter_data += [(i+1, j+1)]
            scatter_info += [alignment_contacts[i][j]]

color_map = {0: '#fa4428', 1: '#fa9057', 2: '#fcd590'}

colors = []
axis_values_x = []
axis_values_y = []

interval = contact_threshold/3
for index, contact in enumerate(scatter_data):
    res_i, res_j = contact
    if alignment_contacts[res_i - 1, res_j - 1] >= 0 and alignment_contacts[res_i - 1, res_j - 1] < interval:
        color = color_map[0]
    elif alignment_contacts[res_i - 1, res_j - 1] >= interval and alignment_contacts[res_i - 1, res_j - 1] < interval*2:
        color = color_map[1]
    else:
        color = color_map[2]
        
    
    colors.append(color)
    axis_values_x.append(res_i)
    axis_values_y.append(res_j)

fig = plt.figure(figsize=(20, 16))
ax = plt.subplot()

def drawplot():
    ax.scatter(
        axis_values_x,
        axis_values_y,
        c=colors,
        picker=True
    )    
    
drawplot()

def annotate(axis, text, x, y):
    annotation = Annotation(text, xy=(x, y), xycoords='data')
    axis.add_artist(annotation)

def onpick(event):
    ind = event.ind
    label_pos_x = event.mouseevent.xdata
    label_pos_y = event.mouseevent.ydata
    offset = 0
    
    for i in ind:
        label = scatter_info[i]
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
    drawplot()
    ax.figure.canvas.draw_idle()

button_clear_all.on_clicked(onclick)

plt.plot()
plt.show()
