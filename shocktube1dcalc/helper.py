#!/usr/bin/env python
import numpy as np
from shocktubecalc import sod

from shocktube1dcalc import solver_analytic

TUBE_DIAPHRAGM = 0.5
TUBE_LENGTH = 1.0
TUBE_LEFT_X = 0.0
TUBE_RIGHT_X = TUBE_LEFT_X + TUBE_LENGTH
TUBE_X_COORDINATE_SHIFT = 0.5

STATE_LEFT_RHO = 1.0
STATE_LEFT_U = 0.0
STATE_LEFT_P = 1.0
STATE_RIGHT_RHO = 0.125
STATE_RIGHT_U = 0.0
STATE_RIGHT_P = 0.1

TIME_STEP_SIZE = 0.01
TIME_TOTAL_ELAPSE = 2.0
MESH_POINTS_NUMBER = 50

DEVIATION_PRECISION = 0.000000001


def convert_format1_to_format2(solutions):
    # convert to shocktubecalc compatible format
    mesh = []
    rho_list = []
    u_list = []
    p_list = []

    for solution_point in solutions:
        mesh.append(solution_point[0])
        rho_list.append(solution_point[1])
        u_list.append(solution_point[2])
        p_list.append(solution_point[3])

    shocktube1d_values = [
        {"Positions": "NA"},
        {"States": "NA"},
        {"x": mesh, "rho": rho_list, "u": u_list, "p": p_list},
    ]

    return shocktube1d_values


def get_cese_mesh(solutions):
    mesh = []
    for solution in solutions:
        mesh.append(solution[0])

    return mesh


def get_shocktube1d_values(time_moment):
    mesh_x_array = np.linspace(
        TUBE_LEFT_X - TUBE_X_COORDINATE_SHIFT,
        TUBE_RIGHT_X - TUBE_X_COORDINATE_SHIFT,
        MESH_POINTS_NUMBER,
    )
    shocktube = solver_analytic.ShockTube()
    analytic_solution = shocktube.get_analytic_solution(mesh_x_array, t=time_moment)
    # convert to shocktubecalc compatible format
    ao_rho_list = []
    ao_u_list = []
    ao_p_list = []

    for solution_point in analytic_solution:
        ao_rho_list.append(solution_point[1])
        ao_u_list.append(solution_point[2])
        ao_p_list.append(solution_point[3])

    mesh_x_array_on_shocktubecalc = np.linspace(
        TUBE_LEFT_X, TUBE_RIGHT_X, MESH_POINTS_NUMBER,
    )
    shocktube1d_values = [
        {"Positions": "NA"},
        {"States": "NA"},
        {
            "x": mesh_x_array_on_shocktubecalc,
            "rho": ao_rho_list,
            "u": ao_u_list,
            "p": ao_p_list,
        },
    ]

    return shocktube1d_values


def get_shocktubecalc_values(time_moment):
    shocktubecalc_values = sod.solve(
        left_state=(STATE_LEFT_P, STATE_LEFT_RHO, STATE_LEFT_U),
        right_state=(STATE_RIGHT_P, STATE_RIGHT_RHO, STATE_RIGHT_U),
        geometry=(TUBE_LEFT_X, TUBE_RIGHT_X, TUBE_DIAPHRAGM),
        t=time_moment,
        gamma=1.4,
        npts=MESH_POINTS_NUMBER,
    )

    return shocktubecalc_values


def get_deviation_values(values_base, values_target, ordered_title):
    values_deviation = {"x": values_target["x"]}

    for key in ordered_title:
        values_deviation_derived = []
        for idx_value in range(len(values_target[key])):
            value_deviation = (
                values_target[key][idx_value] - values_base[key][idx_value]
            )
            values_deviation_derived.append(value_deviation)

            if value_deviation > DEVIATION_PRECISION:
                print(
                    f"{key} deviation {value_deviation} at {values_deviation['x'][idx_value]} "
                    f"is larger than {DEVIATION_PRECISION}"
                )

        values_deviation[key] = values_deviation_derived

    return values_deviation


def check_all_derived_values(base, target, tolerance=1e-15):
    """
    Check the deviation of each element of two data arrays in different derived type (rho, u, or p)

    :param base: dict
    :param target: dict
    :param tolerance: max allowed deviation
    :return: boolean, True for pass
    """
    rtn = True
    for derived_type in base.keys():
        rtn = check_all_values(
            base[derived_type], target[derived_type], tolerance=tolerance
        )

    return rtn


def check_all_values(base, target, tolerance=1e-15):
    """
    Check the deviation of each element of two data arrays.

    :param base: usually data arrays
    :param target: usually data arrays
    :param tolerance: max allowed deviation
    :return: boolean, True for pass
    """
    rtn = True
    # sanity check
    if len(base) != len(target):
        return False

    for idx in range(len(base)):
        deviation = abs(base[idx] - target[idx])
        if deviation > tolerance:
            print("Deviation is larger than tolerance!!")
            rtn = False

    return rtn


def compare(time_step, time_total_elapse):
    time_total_steps = int(time_total_elapse / time_step)
    for idx_step in range(0, time_total_steps):
        moment = TIME_STEP_SIZE * idx_step
        # generate the result from shocktube1d package
        shocktube1d_sodtube = get_shocktube1d_values(moment)

        # generate the result from shocktubecalc package
        shocktubecalc_sodtube = get_shocktubecalc_values(moment)

        shocktube1d_sodtube_values = shocktube1d_sodtube[2]
        shocktubecalc_sodtube_values = shocktubecalc_sodtube[2]

        if not check_all_derived_values(
            shocktube1d_sodtube_values, shocktubecalc_sodtube_values
        ):
            return False

    return True


if __name__ == "__main__":
    compare(TIME_STEP_SIZE, TIME_TOTAL_ELAPSE)
