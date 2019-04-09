def remove_xticklabel(axis):
    labels = [item.get_text() for item in axis.get_xticklabels()]
    empty_string_labels = [''] * len(labels)
    axis.set_xticklabels(empty_string_labels)

