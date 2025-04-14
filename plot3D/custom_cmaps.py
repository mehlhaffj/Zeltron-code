import numpy as np
import os
import matplotlib as mpl

# Function to import custom XML colormaps, of the kind you can download from
# https://sciviscolor.org/colormaps/divergent/
def cmap_from_xml(xml_cmap_file, cmap_name = "mycmap", reverse=False):
    try:
        from bs4 import BeautifulSoup
    except:
        msg = "func cmap_from_xml requires missing module: "
        msg += "bs4.BeautifulSoup"
        raise ValueError(msg)

    try:
        from matplotlib.colors import LinearSegmentedColormap
    except:
        raise ValueError(msg)
        msg = "func cmap_from_xml requires missing module: "
        msg += "matplotlib.colors.LinearSegmentedColormap"

    try:
        from pathlib import Path
    except:
        msg = "func cmap_from_xml requires missing module: "
        msg = "pathlib.Path; are you running on Python 3?"
        raise ValueError(msg)

    path = Path(xml_cmap_file)
    txt = path.read_text()
    soup = BeautifulSoup(txt, "html.parser")
    points = soup.find_all("point")

    levels = np.array( [ float(points[i].get("x")) for i in range(len(points))] )
    colors = np.array( [[float(points[i].get("r")), float(points[i].get("g")),
                        float(points[i].get("b"))] for i in range(len(points))] )

    if reverse:
        colors = colors[::-1]

    return LinearSegmentedColormap.from_list(cmap_name, list(zip(levels, colors)))


# This file defines custom colormaps that can be used with matplotlib
# To define your own colormap, add the corresponding xml file to the
# directory xml_cmaps (which is probably a subdirectory of the directory
# containing this file). When you do this, you will NOT have to modify
# the following line, which simply tells this file where to look for
# your xml files.
xml_maps_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "xml_cmaps")

# After adding your xml file to the xml_cmaps directory, ADD A LINE below
# specifying the name of the color map, which you can then import from this
# file for use in other scripts.
#
# So, for example, in some other plotting script, you would do something like:
#     from custom_cmaps import kestrel 
#     import matplotlib.pyplot as plt
#
#     ...
#     <EPIC analysis code>
#     ...
#    
#     plt.pcolormesh(..., cmap=kestrel, ...)
kestrel = cmap_from_xml(os.path.join(xml_maps_path, "div1-blue-orange-div.xml"), cmap_name = "kestrel")

def register_custom_cmaps(plt):
    # Add lines below to register new custom colormaps with an instance of
    # matplotlib
    # plt.cm.register_cmap(cmap=kestrel)
    mpl.colormaps.register(cmap=kestrel)
