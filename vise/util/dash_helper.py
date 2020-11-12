# -*- coding: utf-8 -*-
#  Copyright (c) 2020 Kumagai group.
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from plotly.graph_objs import Figure


def show_png(fig: Figure):
    fig.write_image("fig.png")
    img = mpimg.imread('fig.png')
    imageplot = plt.imshow(img)
    plt.show()
