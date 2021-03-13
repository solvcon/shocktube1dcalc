#!/usr/bin/env python
# Usage:
#     ./run-animation.py
#
# Description:
#     An example to show the iteration animation
import shocktube1dcalc.cese as cese
import shocktube1dcalc.helper_plot as helper_plot
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig = plt.figure()
frame_seq = []
# iterate
for it_nb in range(0, 100):
    solver_cese = cese.ShockTube(iteration=it_nb)
    solver_cese.run_cese_iteration()
    solution = list(solver_cese.data.solution)
    frame_seq.append(helper_plot.get_gas_status_plot(solution))

ani = animation.ArtistAnimation(fig,
                                frame_seq,
                                interval=25,
                                repeat_delay=300,
                                blit=True)

plt.show()
