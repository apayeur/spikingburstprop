import numpy as np

# user-defined colors for colorblindness
custom_colors = {'blue': (0, 0.45, 0.7),
                 'red': (0.57, 0, 0),
                 'orange': (0.9, 0.6, 0),
                 'vermilion': (0.8, 0.4, 0),
                 'bluish_green': (0, 0.6, 0.5),
                 'yellow': (0.95, 0.9, 0.25),
                 'sky_blue': (0.35, 0.7, 0.9),
                 'pink': (0.8, 0.6, 0.7)}


def remove_xticklabel(axis):
    """
    Removes xtick labels while keeping the ticks themselves.
    """
    labels = [item.get_text() for item in axis.get_xticklabels()]
    empty_string_labels = [''] * len(labels)
    axis.set_xticklabels(empty_string_labels)


def adjust_range(x, y):
    """
    For plotting, scale array x so that min(x) = min(y) and max(x) = max(y)
    :param x: array whose range we want to adjust
    :param y: array
    :return: adjusted array
    """
    min_x = np.min(x)
    max_x = np.max(x)
    min_y = np.min(y)
    max_y = np.max(y)

    return min_y + (x - min_x)*(max_y - min_y)/(max_x - min_x)

