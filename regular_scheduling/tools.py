import re

from bs4 import BeautifulSoup
from IPython.display import SVG, display
from PySimpleAutomata import automata_IO

def draw(dfa,name:str = "dfa",path = "./regular_scheduling/graphics"):
    """Compute a SVG file representing the automaton (see PySimpleAutomata
    documentation for more information).

    Parameters
    ----------
    dfa : dict
        PySimpleAutomata DFA 
    name : str
        name of the future file
    path : str
        path in which the file will appear
    """
    automata_IO.dfa_to_dot(dfa,name,path)

def show(dfa, scale=1.0):
    """Display a representation of the automaton at a chosen scale.

    A SVG file is also computed and save at "./regular_scheduling/graphics/dfa.dot.svg".
    Parameters
    ----------
    dfa : dict
        PySimpleAutomata DFA 
    scale : float
        scale to size the display
    """
    draw(dfa)
    svg = __scale_svg(SVG("./regular_scheduling/graphics/dfa.dot.svg"), scale)
    display(svg)

def __scale_svg(svg_object, scale=1.0):

    soup = BeautifulSoup(svg_object.data, "lxml")
    svg_elt = soup.find("svg")
    w = svg_elt.attrs["width"].rstrip("pt")
    h = svg_elt.attrs["height"].rstrip("pt")

    ws = float(w) * scale
    hs = float(h) * scale

    svg_elt.attrs["width"] = f"{ws}pt"
    svg_elt.attrs["height"] = f"{hs}pt"
    svg_elt.attrs["viewbox"] = f"0.00 0.00 {ws} {hs}"

    g_elt = svg_elt.find("g")
    tf = g_elt.attrs["transform"]
    # non-greedy regex-search-and-replace
    tf2 = re.sub("scale\(.*?\)", f"scale({scale} {scale})", tf)
    g_elt.attrs["transform"] = tf2

    svg_object.data = str(svg_elt)

    return svg_object