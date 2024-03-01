# Altair themes
# (initally) By
# Gregory W. Schwartz

'''
Use:

sys.path.append("/path/to/parent/directory/of/altairThemes.py")

if True:  # In order to bypass isort when saving
    import altairThemes

# register the custom theme under a chosen name
alt.themes.register("publishTheme", altairThemes.publishTheme)

# enable the newly registered theme
alt.themes.enable("publishTheme")
'''

from matplotlib.colors import LinearSegmentedColormap, to_hex, rgb2hex
from typing import List
import altair as alt
import matplotlib.pyplot as plt
import numpy as np


def publishTheme():
    """The publish theme, standardizing some sizes."""

    # Typography
    font = "Arial"
    labelFont = "Arial"
    sourceFont = "Arial"
    fontSize = 10.66670036  # Size 8pt for publishing

    # Axes
    axisColor = "#000000"

    return {
        "config": {
            "view": {
                "strokeOpacity": 0,  # Despine
                "strokeWidth": 0.756,
            },
            "style": {
                "guide-title": {
                    "fontSize": fontSize
                },
                "guide-label": {
                    "fontSize": fontSize
                }
            },
            "title": {
                "fontSize": fontSize,
                "font": font,
                "fontColor": "#000000",
                "fontWeight": "normal",
            },
            "legend": {
                "titleFontWeight": "normal"
            },
            "axis": {
                "domainColor": "#000000",
                "domainWidth": 0.756,
                "grid": False,
                "gridWidth": 0.756,
                "labelAngle": 0,
                "labelFont": labelFont,
                "labelFontSize": fontSize,
                "tickColor": axisColor,
                "tickWidth": 0.756,
                "titleFont": font,
                "titleFontSize": fontSize,
                "titleFontWeight": "normal",
            },
            "header": {
                "titleFontWeight": "normal",
            }
            # For individual axis
            # "axisX": {
            #     "grid": False,
            #     "domainColor": "#000000",
            #     "labelFont": labelFont,
            #     "labelFontSize": fontSize,
            #     "labelAngle": 0,
            #     "tickColor": axisColor,
            #     "titleFont": font,
            #     "titleFontSize": fontSize,
            #     "titleFontWeight": "normal",
            # },
        }
    }

def publishThemeStandardSize():
    """The publish theme but with a standard height for single plots."""
    theme = publishTheme()
    theme["config"]["view"]["width"] = 100
    theme["config"]["view"]["height"] = 80

    return theme


def get_colour_scheme(palette_name: str, num_colours: int) -> List[str]:
    """Extend a colour scheme using colour interpolation.

    By Viet Hoang

    Parameters
    ----------
    palette_name: The matplotlib colour scheme name that will be extended.
    num_colours: The number of colours in the output colour scheme.

    Returns
    -------
    New colour scheme containing 'num_colours' of colours. Each colour is a hex
    colour code.

    """
    scheme = [rgb2hex(c) for c in plt.get_cmap(palette_name).colors]
    if len(scheme) >= num_colours:
        return scheme[:num_colours]
    else:
        cmap = LinearSegmentedColormap.from_list("cmap", scheme)
        extended_scheme = cmap(np.linspace(0, 1, num_colours))
        return [to_hex(c, keep_alpha=False) for c in extended_scheme]
