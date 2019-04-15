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

