#!/usr/bin/env python
import collections

import matplotlib.animation as animation
import matplotlib.pyplot as plt

from shocktube1dcalc import cese, helper, solver_analytic
from shocktube1dcalc.helper import convert_format1_to_format2, get_deviation_values

TIME_STEP_SIZE = 0.004
TIME_TOTAL_ELAPSE = 0.4


def get_analytic_solutions(moment, mesh):
    shocktube = solver_analytic.ShockTube()

    return shocktube.get_analytic_solution(mesh, t=moment)


def get_cese_solutions(moment):
    cese_grid_size_t = 0.004
    # multiply 2 for half grids, so total iteration number should be double
    iteration_number = round(moment / cese_grid_size_t * 2)
    shocktube = cese.ShockTube(iteration=iteration_number, grid_size_t=cese_grid_size_t)
    shocktube.run_cese_iteration()

    return shocktube.data.solution


def plot_solution_single_frame_overlapping_init(
    frame,
    titles_ordered,
    time_moment,
    type_name,
    ax_type,
    values_base,
    values_target,
    color_base,
    color_target,
    marker_base,
    marker_target
):
    key_idx = 0
    containers = []
    for key in titles_ordered:
        ax = ax_type[key_idx]

        container_base = ax.scatter(
            values_base["x"], values_base[key], s=14, c=color_base, marker=marker_base, label="Analytic"
        )

        container_target = ax.scatter(
            values_target["x"], values_target[key], s=2, c=color_target, marker=marker_target, label="CESE"
        )

        if type_name == "deviation":
            ax.set(xlim=[-1, 1], ylim=[-0.01, 0.01])
        else:
            ax.set(xlim=[-1, 1], ylim=[0, 1.1])

        ax.set_title(titles_ordered[key] + f" ({type_name})")

        text_time = plt.text(
            0.1, 0.08, f"Time: {time_moment:.2f}", transform=ax.transAxes
        )

        frame.append(container_base)
        frame.append(container_target)
        frame.append(text_time)
        key_idx = key_idx + 1

        containers.append(container_base)
        containers.append(container_target)

    return containers


def plot_solution_single_frame_overlapping(
    frame,
    titles_ordered,
    time_moment,
    type_name,
    ax_type,
    values_base,
    values_target,
    color_base,
    color_target,
    marker_base,
    marker_target
):
    plot_solution_single_frame_overlapping_init(
        frame,
        titles_ordered,
        time_moment,
        type_name,
        ax_type,
        values_base,
        values_target,
        color_base,
        color_target,
        marker_base,
        marker_target
    )


def plot_solution_single_frame(
    frame, titles_ordered, time_moment, type_name, ax_type, values_type, color, marker
):
    key_idx = 0
    for key in titles_ordered:
        ax = ax_type[key_idx]
        subplot = ax.scatter(
            values_type["x"], values_type[key], s=10, c=color, marker=marker
        )

        if type_name == "deviation":
            ax.set(xlim=[-1, 1], ylim=[-0.4, 0.4])
        else:
            ax.set(xlim=[-1, 1], ylim=[0, 1.1])

        ax.set_title(titles_ordered[key] + f" ({type_name})")
        text_time = plt.text(
            0.1, 0.08, f"Time: {time_moment:.2f}", transform=ax.transAxes
        )

        frame.append(subplot)
        frame.append(text_time)
        key_idx = key_idx + 1


def get_solution_single_artist_frame(
    values_base, values_target, ax_overlapping, ax_deviation, moment
):
    frame = []

    title_items = [
        ("p", "Pressure"),
        ("rho", "Density"),
        ("u", "Velocity"),
    ]
    titles_ordered = collections.OrderedDict(title_items)

    helper.DEVIATION_PRECISION = 0.4
    values_deviation = get_deviation_values(values_base, values_target, titles_ordered)

    plot_solution_single_frame_overlapping(
        frame, titles_ordered, moment, "analytic vs. CESE", ax_overlapping, values_base, values_target,
        "g", "b", "X", "x"
    )
    plot_solution_single_frame(
        frame,
        titles_ordered,
        moment,
        "deviation",
        ax_deviation,
        values_deviation,
        "r",
        "o",
    )

    return frame


def get_solution_video_frames(
    time_step, time_total_elapse, ax_overlapping, ax_deviation
):
    artist_frames = []
    time_total_steps = int(time_total_elapse / time_step)
    for idx_step in range(0, time_total_steps):
        moment = TIME_STEP_SIZE * idx_step

        # generate the CESE solutions
        sodtube_cese = get_cese_solutions(moment)
        mesh_cese = helper.get_cese_mesh(sodtube_cese)
        # generate the analytic result from shocktube1d package
        sodtube_analytic = get_analytic_solutions(moment, mesh_cese)

        solution_cese = convert_format1_to_format2(sodtube_cese)[2]
        solution_analytic = convert_format1_to_format2(sodtube_analytic)[2]

        artist_frame = get_solution_single_artist_frame(
            solution_analytic,
            solution_cese,
            ax_overlapping,
            ax_deviation,
            moment,
        )
        artist_frames.append(artist_frame)

    return artist_frames


# Output video
#
# Set up formatting for the movie files
subplot_row = 2
subplot_column = 3
Writer = animation.writers["ffmpeg"]
writer = Writer()
my_dpi = 96
fig4video, (axis_overlapping, axis_deviation) = plt.subplots(
    subplot_row, subplot_column, figsize=(1600 / my_dpi, 1000 / my_dpi), dpi=my_dpi
)
fig4video.subplots_adjust(hspace=0.4, wspace=0.4)
frames = get_solution_video_frames(
    TIME_STEP_SIZE, TIME_TOTAL_ELAPSE, axis_overlapping, axis_deviation
)

ani = animation.ArtistAnimation(
    fig4video, frames, interval=250, repeat_delay=1000
)
ani.save("/tmp/1d-sod-tube-cese-analytic.mp4", writer=writer)
