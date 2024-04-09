import json
import os

import numpy as np

import parasbolv as psv
import matplotlib.pyplot as plt
from collections import namedtuple

"""
Visualization of Genetic Logic Circuits consisting of two plasmids with SBOL Visual on the basis of paraSBOLv.

Reference paraSBOLv
Clark C.J., Scott-Brown J. & Gorochowski T.E. "paraSBOLv: a foundation for standard-compliant genetic design visualisation tools", Synthetic Biology, 2021 doi:10.1093/synbio/ysab022

paraSBOLv has been slightly adapted to allow for the needs of multi plasmid visualizations. In particular, it returns the bound_list of the glyphs included in a construct besides the bounds of the construct itself.
"""

Part = namedtuple('part', ['glyph_type', 'orientation', 'user_parameters', 'style_parameters'])
Interaction = namedtuple('interaction',
                         ['starting_glyph', 'ending_glyph', 'interaction_type', 'interaction_parameters'])


def visualize_plasmids(plasmids, name, output_dir):
    plt.figure()
    pass


def multi_plasmid():
    part_list = []
    part_list.append(Part("Promoter", "forward", None, None))
    part_list.append(Part('RibosomeEntrySite', 'forward', None, None))
    part_list.append(Part('CDS',
                          'forward',
                          None,
                          {'cds': {'facecolor': (1, 0.5, 0.5), 'edgecolor': (1, 0, 0), 'linewidth': 2}}
                          ))

    part_list.append(Part('Terminator', 'forward', None, None))
    part_list.append(Part("Promoter", "forward", None, None))

    part_list.append(Part('RibosomeEntrySite', 'forward', None, None))
    part_list.append(Part('CDS',
                          'forward',
                          None,
                          {'cds': {'facecolor': (0.5, 0.5, 1), 'edgecolor': (0, 0, 1), 'linewidth': 2}}
                          ))

    part_list.append(Part('Terminator', 'forward', None, None))

    # Create renderer
    renderer = psv.GlyphRenderer()

    # Create list of interactions to pass to render_part_list
    interaction_list = []
    interaction_list.append(Interaction(part_list[2], part_list[4], 'inhibition', {'color': (0.75, 0, 0)}))
    # Multi Plasmid Example

    fig, axes = plt.subplots(ncols=1, nrows=1)
    ax = axes

    start_position = (0, 0)
    construct = psv.Construct(part_list, renderer, interaction_list=interaction_list, fig=fig, ax=ax,
                              start_position=start_position)
    fig, ax, baseline_start, baseline_end, bounds1, bounds_list1 = construct.draw(draw_for_bounds=False)
    ax.plot([baseline_start[0], baseline_end[0]], [baseline_start[1], baseline_end[1]], color=(0, 0, 0), linewidth=1.5,
            zorder=0)

    x_lim1 = ax.get_xlim()
    y_lim1 = ax.get_ylim()

    start_position = (0, -80)
    construct = psv.Construct(part_list, renderer, interaction_list=interaction_list, fig=fig, ax=ax,
                              start_position=start_position)
    fig, ax, baseline_start, baseline_end, bounds2, bounds_list2 = construct.draw(draw_for_bounds=False)
    ax.plot([baseline_start[0], baseline_end[0]], [baseline_start[1], baseline_end[1]], color=(0, 0, 0), linewidth=1.5,
            zorder=0)

    # psv.draw_interaction()

    x_lim2 = ax.get_xlim()
    y_lim2 = ax.get_ylim()

    x_lims = [x_lim1, x_lim2]
    y_lims = [y_lim1, y_lim2]

    x_start = np.min(x_lims)
    x_end = np.max(x_lims)
    y_start = np.min(y_lims)
    y_end = np.max(y_lims)

    ax.set_xlim([x_start, x_end])
    ax.set_ylim([y_start, y_end])

    # plt.savefig("multi_plasmid_example.pdf", dpi=300)

    # plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    plasmid_directory = "data/11111101/plasmids/"
    output_dir = plasmid_directory

    multi_plasmid()
    exit(0)

    for file in os.listdir(plasmid_directory):
        file_path = plasmid_directory + file
        name, extension = os.path.splitext(file)

        if extension != ".json":
            continue

        print(file, extension)

        with open(file_path, "r") as file:
            plasmids = json.load(file)
        visualize_plasmids(plasmids, name, output_dir)
    pass
