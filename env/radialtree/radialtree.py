#!/usr/bin/env python3

import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
from matplotlib.axes import Axes
import matplotlib
from matplotlib.lines import Line2D


colormap_list = [
    "nipy_spectral",
    "terrain",
    "gist_rainbow",
    "CMRmap",
    "coolwarm",
    "gnuplot",
    "gist_stern",
    "brg",
    "rainbow",
]


def radialTreee(
    Z2,
    fontsize=8,
    ax: Axes = None,
    pallete="gist_rainbow",
    addlabels=True,
    sample_classes=None,
    colorlabels=None,
    colorlabels_legend=None,
):
    """
    Drawing a radial dendrogram from a scipy dendrogram output.
    Parameters
    ----------
    Z2 : dictionary
        A dictionary returned by scipy.cluster.hierarchy.dendrogram
    addlabels: bool
        A bool to choose if labels are shown.
    fontsize : float
        A float to specify the font size

    ax : Axes or None:
        Axes in which to draw the plot, otherwise use the currently-active Axes.
    pallete : string or None
        Matlab colormap name.
        If `None` is provided then `color_list` from Z2 is used as is.
    sample_classes : dict
        A dictionary that contains lists of sample subtypes or classes. These classes appear
        as color labels of each leaf. Colormaps are automatically assigned. Not compatible
        with options "colorlabels" and "colorlabels_legend".
        e.g., {"color1":["Class1","Class2","Class1","Class3", ....]}
    colorlabels : dict
        A dictionary to set color labels to leaves. The key is the name of the color label.
        The value is the list of RGB color codes, each corresponds to the color of a leaf.
        e.g., {"color1":[[1,0,0,1], ....]}
    colorlabels_legend : dict
        A nested dictionary to generate the legends of color labels. The key is the name of
        the color label. The value is a dictionary that has two keys "colors" and "labels".
        The value of "colors" is the list of RGB color codes, each corresponds to the class of a leaf.
        e.g., {"color1":{"colors":[[1,0,0,1], ....], "labels":["label1","label2",...]}}
    Returns
    -------
    Raises
    ------
    Notes
    -----
    References
    ----------
    See Also
    --------
    Examples
    --------
    """
    if ax is None:
        ax: Axes = plt.gca()

    linewidth = 0.5
    R = 1
    width = R * 0.1
    space = R * 0.05

    if colorlabels != None:
        offset = (
            width * len(colorlabels) / R + space * (len(colorlabels) - 1) / R + 0.05
        )
    elif sample_classes != None:
        offset = (
            width * len(sample_classes) / R
            + space * (len(sample_classes) - 1) / R
            + 0.05
        )
    else:
        offset = 0

    xmax = np.amax(Z2["icoord"])
    xmin = np.amin(Z2["icoord"])
    ymax = np.amax(Z2["dcoord"])

    # number of leaves
    n_leaves = len(Z2["ivl"])

    ucolors = sorted(set(Z2["color_list"]))
    if pallete:
        cmp = cm.get_cmap(pallete, len(ucolors))
        if type(cmp) == matplotlib.colors.LinearSegmentedColormap:
            cmap = cmp(np.linspace(0, 1, len(ucolors)))
        else:
            cmap = cmp.colors
        def get_color(c):
            return cmap[ucolors.index(c)]
    else:
        def get_color(c):
            return c
        
    # create a mapping from original x-coordinates to evenly spaced angles
    x_to_angle = {}
    leaf_positions = []

    # calculate positions for leaves
    for i, label in enumerate(Z2["ivl"]):
        # original x position in scipy dendrogram
        orig_x = 5.0 + i * 10.0
        # convert to radians, evenly spaced around the circle
        angle = 2 * np.pi * i / n_leaves
        x_to_angle[orig_x] = angle
        leaf_positions.append((orig_x, angle))

    # for internal nodes, interpolate angles based on their children's positions
    for coords in Z2["icoord"]:
        x1, x2 = coords[0], coords[2]  # the x positions of the two children
        
        if x1 not in x_to_angle:
            
            # find closest leaf positions and interpolate
            closest_positions = sorted(leaf_positions, key=lambda pos: abs(pos[0] - x1))
            if closest_positions:
                closest_x, closest_angle = closest_positions[0]
                x_to_angle[x1] = closest_angle
        
        if x2 not in x_to_angle:
            
            closest_positions = sorted(leaf_positions, key=lambda pos: abs(pos[0] - x2))
            if closest_positions:
                closest_x, closest_angle = closest_positions[0]
                x_to_angle[x2] = closest_angle
    
    # draw the dendrogram with the angular mapping
    for icoord, dcoord, c in sorted(zip(Z2["icoord"], Z2["dcoord"], Z2["color_list"])):
        _color = get_color(c)        
        if c == "C0":
            _color = "black"
        
        # map original x coordinates to angles and then to x,y positions
        r = R * (1 - np.array(dcoord) / ymax)
        
        # get the angles for the two endpoints using our mapping
        angle1 = x_to_angle.get(icoord[0], 2 * np.pi * icoord[0] / xmax)
        angle2 = x_to_angle.get(icoord[2], 2 * np.pi * icoord[2] / xmax)
        
        # calculate x,y coordinates for the four points
        _xr0 = np.cos(angle1) * r[0]
        _yr0 = np.sin(angle1) * r[0]
        _xr1 = np.cos(angle1) * r[1]
        _yr1 = np.sin(angle1) * r[1]
        _xr2 = np.cos(angle2) * r[2]
        _yr2 = np.sin(angle2) * r[2]
        _xr3 = np.cos(angle2) * r[3]
        _yr3 = np.sin(angle2) * r[3]
        
        # plot radial lines
        ax.plot([_xr0, _xr1], [_yr0, _yr1], c=_color, linewidth=linewidth)
        ax.plot([_xr2, _xr3], [_yr2, _yr3], c=_color, linewidth=linewidth)
        
        # determine how to draw the arc connecting branches
        arc_radius = r[1]  # use the inner radius for the arc
        
        # calculate angular difference, ensuring we take the shorter path
        angular_diff = angle2 - angle1
        if angular_diff > np.pi:
            angular_diff -= 2 * np.pi
        elif angular_diff < -np.pi:
            angular_diff += 2 * np.pi
        
        # draw the arc with multiple segments for smoothness
        num_segments = 50
        angles = np.linspace(angle1, angle2, num_segments)
        x_arc = arc_radius * np.cos(angles)
        y_arc = arc_radius * np.sin(angles)
        ax.plot(x_arc, y_arc, c=_color, linewidth=linewidth)
    
    # calculate evenly spaced label positions
    label_coords = []
    for i, label in enumerate(Z2["ivl"]):
        # evenly distribute angles around the circle
        angle = 2 * np.pi * i / n_leaves
        # calculate base position with offset for labels
        x_pos = np.cos(angle) * (1.05 + offset)
        y_pos = np.sin(angle) * (1.05 + offset)
        # determine rotation based on position
        # for left side of circle, adjust text rotation for readability
        if x_pos < 0:
            rotation = angle * 180 / np.pi + 180  # Convert to degrees and flip
        else:
            rotation = angle * 180 / np.pi  # Convert to degrees
        label_coords.append([x_pos, y_pos, rotation])
    
    if addlabels == True:
        assert len(Z2["ivl"]) == len(label_coords), (
            f'Internal error, label numbers for Z2 ({len(Z2["ivl"])})'
            f" and for calculated labels ({len(label_coords)}) must be equal!"
        )
        for (_x, _y, _rot), label in zip(label_coords, Z2["ivl"]):
            # determine text alignment based on position
            if _x < 0:
                ha = "right"
            else:
                ha = "left"
            ax.text(
                _x,
                _y,
                label,
                {"va": "center", "ha": ha},
                rotation_mode="anchor",
                rotation=_rot,
                fontsize=fontsize,
            )
    # handle color labels (wedges)
    if colorlabels != None or sample_classes != None:
        
        # compute even pie intervals for wedges
        intervals = np.ones(n_leaves) * (2 * np.pi / n_leaves)
        
        if colorlabels != None:
            j = 0
            outerrad = R * 1.05 + width * len(colorlabels) + space * (len(colorlabels) - 1)
            labelnames = []
        
            for labelname, colorlist in colorlabels.items():
                colorlist = np.array(colorlist)[Z2["leaves"]]
                outerrad = outerrad - width * j - space * j
                
                # draw the pie chart with even wedges
                patches, texts = ax.pie(
                    intervals,
                    colors=colorlist,
                    radius=outerrad,
                    counterclock=True,
                    startangle=360 / (2 * n_leaves),
                    wedgeprops=dict(width=width),
                )
                
                labelnames.append(labelname)
                j += 1
            
            # handle legend
            if colorlabels_legend != None:
                for i, labelname in enumerate(labelnames):
                    colorlines = []
                    for c in colorlabels_legend[labelname]["colors"]:
                        colorlines.append(Line2D([0], [0], color=c, lw=4))
                    leg = ax.legend(
                        colorlines,
                        colorlabels_legend[labelname]["labels"],
                        bbox_to_anchor=(1.5 + 0.3 * i, 1.0),
                        title=labelname,
                    )
                    ax.add_artist(leg)
        
        elif sample_classes != None:
            j = 0
            outerrad = R * 1.05 + width * len(sample_classes) + space * (len(sample_classes) - 1)
            labelnames = []
            colorlabels_legend = {}
            
            for labelname, colorlist in sample_classes.items():
                ucolors = sorted(list(np.unique(colorlist)))
                type_num = len(ucolors)
                _cmp = cm.get_cmap(colormap_list[j], type_num)
                _colorlist = [_cmp(ucolors.index(c)) for c in colorlist]
                _colorlist = np.array(_colorlist)[Z2["leaves"]]
                outerrad = outerrad - width * j - space * j
                
                # draw the pie chart with even wedges
                patches, texts = ax.pie(
                    intervals,
                    colors=_colorlist,
                    radius=outerrad,
                    counterclock=True,
                    startangle=360 / (2 * n_leaves),
                    wedgeprops=dict(width=width),
                )
                labelnames.append(labelname)
                colorlabels_legend[labelname] = {}
                colorlabels_legend[labelname]["colors"] = _cmp(np.linspace(0, 1, type_num))
                colorlabels_legend[labelname]["labels"] = ucolors
                j += 1
            
            # handle legend
            for i, labelname in enumerate(labelnames):
                colorlines = []
                for c in colorlabels_legend[labelname]["colors"]:
                    colorlines.append(Line2D([0], [0], color=c, lw=4))
                leg = ax.legend(
                    colorlines,
                    colorlabels_legend[labelname]["labels"],
                    bbox_to_anchor=(1.5 + 0.3 * i, 1.0),
                    title=labelname,
                )
                ax.add_artist(leg)
    
    # remove spines and ticks
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.spines.left.set_visible(False)
    ax.spines.bottom.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    
    # set axis limits
    if colorlabels != None:
        maxr = R * 1.05 + width * len(colorlabels) + space * (len(colorlabels) - 1)
    elif sample_classes != None:
        maxr = R * 1.05 + width * len(sample_classes) + space * (len(sample_classes) - 1)
    else:
        maxr = R * 1.05
    ax.set_xlim(-maxr, maxr)
    ax.set_ylim(-maxr, maxr)
    
    return ax


def plot(
    Z2,
    fontsize=8,
    figsize=None,
    pallete="gist_rainbow",
    addlabels=True,
    show=True,
    sample_classes=None,
    colorlabels=None,
    colorlabels_legend=None,
):
    """
    Drawing a radial dendrogram from a scipy dendrogram output.
    Parameters
    ----------
    Z2 : dictionary
        A dictionary returned by scipy.cluster.hierarchy.dendrogram
    addlabels: bool
        A bool to choose if labels are shown.
    fontsize : float
        A float to specify the font size
    figsize : [x, y] array-like
        1D array-like of floats to specify the figure size
    pallete : string
        Matlab colormap name.
    sample_classes : dict
        A dictionary that contains lists of sample subtypes or classes. These classes appear
        as color labels of each leaf. Colormaps are automatically assigned. Not compatible
        with options "colorlabels" and "colorlabels_legend".
        e.g., {"color1":["Class1","Class2","Class1","Class3", ....]}
    colorlabels : dict
        A dictionary to set color labels to leaves. The key is the name of the color label.
        The value is the list of RGB color codes, each corresponds to the color of a leaf.
        e.g., {"color1":[[1,0,0,1], ....]}
    colorlabels_legend : dict
        A nested dictionary to generate the legends of color labels. The key is the name of
        the color label. The value is a dictionary that has two keys "colors" and "labels".
        The value of "colors" is the list of RGB color codes, each corresponds to the class of a leaf.
        e.g., {"color1":{"colors":[[1,0,0,1], ....], "labels":["label1","label2",...]}}

    Returns
    -------
    Raises
    ------
    Notes
    -----
    References
    ----------
    See Also
    --------
    Examples
    --------
    """
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial"]
    plt.rcParams["svg.fonttype"] = "none"

    if figsize == None and colorlabels != None:
        figsize = [10, 5]
    elif figsize == None and sample_classes != None:
        figsize = [10, 5]
    elif figsize == None:
        figsize = [5, 5]
    
    fig, ax = plt.subplots(figsize=figsize)
    
    ax = radialTreee(
        Z2,
        fontsize=fontsize,
        ax=ax,
        pallete=pallete,
        addlabels=addlabels,
        sample_classes=sample_classes,
        colorlabels=colorlabels,
        colorlabels_legend=colorlabels_legend,
    )
    
    if show == True:
        fig.show()
    else:
        return ax


def mat_plot(mat):
    # Take a matrix data instead of a dendrogram data, calculate dendrogram and draw a circular dendrogram
    pass


def pandas_plot(df):

    pass


def _test_1(Z2):
    # optionally leaves can be labeled by colors
    type_num = 12
    _cmp = cm.get_cmap("bwr", type_num)
    _cmp2 = cm.get_cmap("hot", type_num)
    colors_dict = {
        "example_color": _cmp(np.random.rand(numleaf)),
        "example_color2": _cmp2(np.random.rand(numleaf)),
    }
    colors_legends = {
        "example_color": {
            "colors": _cmp(np.linspace(0, 1, type_num)),
            "labels": ["ex1_" + str(i + 1) for i in range(type_num)],
        },
        "example_color2": {
            "colors": _cmp2(np.linspace(0, 1, type_num)),
            "labels": ["ex2_" + str(i + 1) for i in range(type_num)],
        },
    }
    # fig = pylab.figure(figsize=(8,8))

    # Compute and plot the dendrogram.
    # ax2 = fig.add_axes([0.3,0.71,0.6,0.2])

    fig, ax = plt.subplots(figsize=(10, 5))
    # plot(Z2, colorlabels=colors_dict,colorlabels_legend=colors_legends,show=True)
    radialTreee(Z2, ax=ax, colorlabels=colors_dict, colorlabels_legend=colors_legends)
    fig.show()


def _test_2(Z2):
    type_num = 6
    type_list = ["ex" + str(i) for i in range(type_num)]
    sample_classes = {
        "example_color": [np.random.choice(type_list) for i in range(numleaf)]
    }
    fig, ax = plt.subplots(figsize=(10, 5))
    radialTreee(Z2, ax=ax, sample_classes=sample_classes)
    fig.show()
    # plot(Z2, sample_classes=sample_classes,show=True)


def _test_3(Z2):
    fig, ax = plt.subplots(2, 2, figsize=(10, 10))
    ax = ax.flatten()
    # no arguments
    radialTreee(Z2, ax=ax[0])
    ax[0].set_aspect(1)

    type_num = 12
    _cmp = cm.get_cmap("bwr", type_num)
    _cmp2 = cm.get_cmap("hot", type_num)
    colors_dict = {
        "example_color": _cmp(np.random.rand(numleaf)),
        "example_color2": _cmp2(np.random.rand(numleaf)),
    }
    colors_legends = {
        "example_color": {
            "colors": _cmp(np.linspace(0, 1, type_num)),
            "labels": ["ex1_" + str(i + 1) for i in range(type_num)],
        },
        "example_color2": {
            "colors": _cmp2(np.linspace(0, 1, type_num)),
            "labels": ["ex2_" + str(i + 1) for i in range(type_num)],
        },
    }
    # fig = pylab.figure(figsize=(8,8))

    # Compute and plot the dendrogram.
    # ax2 = fig.add_axes([0.3,0.71,0.6,0.2])

    # like in test_1
    radialTreee(
        Z2, ax=ax[1], colorlabels=colors_dict, colorlabels_legend=colors_legends
    )

    type_num = 6
    type_list = ["ex" + str(i) for i in range(type_num)]
    sample_classes = {
        "example_color": [np.random.choice(type_list) for i in range(numleaf)]
    }
    radialTreee(Z2, ax=ax[2], sample_classes=sample_classes)
    ax[3].axis("off")
    fig.show()


if __name__ == "__main__":
    # Generate random features and distance matrix.

    test = [0, 1, 2, 3]
    np.random.seed(1)
    numleaf = 42
    _alphabets = [chr(i) for i in range(97, 97 + 24)]
    labels = sorted(
        ["".join(list(np.random.choice(_alphabets, 10))) for i in range(numleaf)]
    )

    x = np.random.rand(numleaf)
    D = np.zeros([numleaf, numleaf])
    for i in range(numleaf):
        for j in range(numleaf):

            D[i, j] = abs(x[i] - x[j])
    Y = sch.linkage(D, method="single")
    Z2 = sch.dendrogram(Y, labels=labels, no_plot=True)

    if 3 in test:
        _test_3(Z2)

    if 0 in test:
        plot(Z2, show=True)

    if 1 in test:
        _test_1(Z2)

    if 2 in test:
        _test_2(Z2)

    plt.show()