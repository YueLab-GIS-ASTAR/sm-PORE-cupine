import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import numpy as np

def get_list_of_arcs(dotbracket):
    list_of_coords = []
    curr_list = []
    for n2, i in enumerate(dotbracket):
        if i == "(":
            curr_list.append(n2)
        elif i == ")":
            n1 = curr_list.pop(-1)
            list_of_coords.append((n1, n2))
    return list_of_coords

def rgb_to_hex(rgb):
    r, g, b = map(int, rgb.split(","))
    def clamp(x): 
        return max(0, min(x, 255))
    return "#{0:02x}{1:02x}{2:02x}".format(clamp(r), clamp(g), clamp(b))

def calc_height(left, right):
    x1, y1 = left, 0.0
    x2, y2 = right, 0.0
    mxmy = mx, my = [(x1 + x2) / 2, (y1 + y2) / 2]
    r = np.sqrt((x1 - mx)**2 + (y1 - my)**2)
    width = 2 * r
    height = 2 * r
    start_angle = np.arctan2(y1 - my, x1 - mx) * 180 / np.pi
    end_angle = np.arctan2(my - y2, mx - x2) * 180 / np.pi
    return height


def plot_arcs_from_list(list_of_coords,
                        list_of_colors=None, p_color="black",
                        list_of_counts=None,
                        p_linewidth=0.5, p_adjust_height=1.5,
                        p_figax=None, p_figsize=(12, 6), p_show_count=False, p_offset=0):
    
    if list_of_colors is None:
        list_of_colors = [ p_color for _ in range(len(list_of_coords)) ]
    
    if list_of_counts is None:
        list_of_counts = [ 1 for _ in range(len(list_of_coords)) ]

    if p_figax is None:
        fig, ax = plt.subplots(figsize=p_figsize)
    else:
        fig, ax = p_figax
        
    heights, starts, ends = [], [], []
    for (start, end), color, count in zip(list_of_coords, list_of_colors, list_of_counts):
        width = end - start
        half_height = calc_height(start-p_offset, end-p_offset) / 2
        
        center = start - p_offset + width / 2
        ax.add_patch(Arc((center, 0), width,
                    2 * half_height, 0, 0, 180, color=color,
                    linewidth=p_linewidth,
                    ls="solid"))

        if p_show_count:
            ax.annotate(count, (center, half_height+5))

        heights.append(half_height)
        starts.append(start-p_offset)
        ends.append(end-p_offset)

    max_height = max(heights) * p_adjust_height
    xlim_begin = min(starts)
    xlim_end = max(ends)
    
    #ax.set_xlim(xlim_begin, xlim_end)
    ax.set_xlim(50, 29800)
    ax.set_ylim(0, max_height)

    return fig, ax



def plot_arcs_from_list_HJ(list_of_coords,
                               list_of_colors=None, p_color="steelblue",
                               list_of_counts=None,
                               p_linewidth=0.5, p_adjust_height=1.5,
                               p_figax=None, p_figsize=(12, 6), p_show_count=False, p_offset=0):
    
    if list_of_colors is None:
        list_of_colors = [ p_color for _ in range(len(list_of_coords)) ]
    
    if list_of_counts is None:
        list_of_counts = [ 1 for _ in range(len(list_of_coords)) ]

    if p_figax is None:
        fig, ax = plt.subplots(figsize=p_figsize)
    else:
        fig, ax = p_figax
        
    heights, starts, ends = [], [], []
    for (start, end), color, count in zip(list_of_coords, list_of_colors, list_of_counts):
        width = end - start
        half_height = calc_height(start-p_offset, end-p_offset) / 2


        custom_widths = count/50
        center = start - p_offset + width / 2
        ax.add_patch(Arc((center, 0), width,
                    2 * half_height, 0, 0, 180, color=p_color, alpha=0.5,
                    linewidth=custom_widths,
                    ls="solid"))

        if p_show_count:
            ax.annotate(count, (center, half_height+5))

        heights.append(half_height)
        starts.append(start-p_offset)
        ends.append(end-p_offset)

    max_height = max(heights) * p_adjust_height
    xlim_begin = min(starts)
    xlim_end = max(ends)
    
    #ax.set_xlim(xlim_begin, xlim_end)
    ax.set_xlim(50, 29800)
    ax.set_ylim(0, max_height)

    return fig, ax



def plot_arcs_from_list_subgenomic(list_of_coords, 
                        list_of_colors=None, p_color="black",
                        list_of_counts=None,
                        p_subgenomic_start=28200,
                        p_linewidth=0.5, p_adjust_height=1.5,
                        p_figax=None, p_figsize=(12, 6), p_show_count=False, p_offset=0):
    
    if list_of_colors is None:
        list_of_colors = [ p_color for _ in range(len(list_of_coords)) ]
    
    if list_of_counts is None:
        list_of_counts = [ 1 for _ in range(len(list_of_coords)) ]

    if p_figax is None:
        fig, ax = plt.subplots(figsize=p_figsize)
    else:
        fig, ax = p_figax
        
    heights, starts, ends = [], [], []
    for (start, end), color, count in zip(list_of_coords, list_of_colors, list_of_counts):
        width = end - start
        half_height = calc_height(start-p_offset, end-p_offset) / 2
        
        center = start - p_offset + width / 2
        ax.add_patch(Arc((center, 0), width,
                    7 * half_height, 0, 0, 180, color=color,
                    linewidth=p_linewidth,
                    ls="solid"))

        if p_show_count:
            ax.annotate(count, (center, half_height+5))

        heights.append(half_height)
        starts.append(start-p_offset)
        ends.append(end-p_offset)

    max_height = max(heights) * p_adjust_height
    xlim_begin = min(starts)
    xlim_end = max(ends)
    
    #ax.set_xlim(xlim_begin, xlim_end)
    ax.set_xlim(p_subgenomic_start, 29800)
    ax.set_ylim(0, max_height)

    return fig, ax




def plot_arcs_from_list_truncate_coordinate(list_of_coords,
                                            list_of_colors=None, p_color="black",
                                            list_of_counts=None,
                                            p_linewidth=0.5, p_adjust_height=1.5,
                                            p_figax=None, p_figsize=(12, 6), p_show_count=False, p_offset=0, p_truncate=50):
           
    if list_of_colors is None:
        list_of_colors = [ p_color for _ in range(len(list_of_coords)) ]
    
    if list_of_counts is None:
        list_of_counts = [ 1 for _ in range(len(list_of_coords)) ]

    if p_figax is None:
        fig, ax = plt.subplots(figsize=p_figsize)
    else:
        fig, ax = p_figax
        
    heights, starts, ends = [], [], []
    for (start, end), color, count in zip(list_of_coords, list_of_colors, list_of_counts):
        width = end - start
        half_height = calc_height(start-p_offset, end-p_offset) / 2
        
        center = start - p_offset + width / 2
        ax.add_patch(Arc((center, 0), width,
                    2 * half_height, 0, 0, 180, color=color,
                    linewidth=p_linewidth,
                    ls="solid"))

        if p_show_count:
            ax.annotate(count, (center, half_height+5))

        heights.append(half_height)
        starts.append(start-p_offset)
        ends.append(end-p_offset)

    max_height = max(heights) * p_adjust_height
    xlim_begin = min(starts)
    xlim_end = max(ends)
    
    #ax.set_xlim(xlim_begin, xlim_end)
    ax.set_xlim(25280, 29800)
    ax.set_ylim(0, max_height)

    return fig, ax



def plot_arcs_from_bed(f_bed, p_figax=None, p_figsize=(12, 6)):

    if p_figax is None:
        fig, ax = plt.subplots(figsize=p_figsize)
    else:
        fig, ax = p_figax
    
    list_of_coords = []
    list_of_colors = []
    with open(f_bed, "r") as f:
        for line in f:
            if line.startswith("track"):
                continue
            row = line.strip("\r\n").split()
            start, end, color = int(row[1]), int(row[2]), row[8]

            list_of_coords.append((start, end))
            list_of_colors.append(rgb_to_hex(color))
    
    fig, ax = plot_arcs_from_list(list_of_coords,
                                  list_of_colors=list_of_colors,
                                  p_figax=(fig, ax))
    return fig, ax
    
