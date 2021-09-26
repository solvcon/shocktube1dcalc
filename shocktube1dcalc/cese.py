# Description:
#   1D Sod shock tube solver based on CESE method.
#
#   The derivation of the equations for the analytic solution
#   is based on the book,
#   Principles of Computational Fluid Dynamics,
#   written by Pieter Wesseling.
#

import numpy as np

from shocktube1dcalc import generator_mesh

gmesh = generator_mesh.Mesher()

GAMMA = 1.4

# Sod tube initial conditions of the left and right
# of the diaphragm
RHO_L = 1.0
U_L = 0.0
P_L = 1.0
RHO_R = 0.125
U_R = 0.0
P_R = 0.1


class ShockTube(object):
    """
    CESE method to generate the 1D Sod tube solution
    @iteration, which is integer
    TODO:
        1. grid_size_x, mesh_x_start, and mesh_x_stop should be integer
        for gmesh.gen_mesh method, which requires integer inputs to be easier
        to generate mesh points along x. This should be smarter and allow user
        to input their own numbers. Users' inputs are 10000 times as
        the real inputs for CESE interation.

    """

    def __init__(
        self,
        iteration=100,
        grid_size_t=0.004,
        grid_size_x=0.01,
        mesh_x_start=-1.0,
        mesh_x_stop=1.0,
        rho_l=RHO_L,
        u_l=U_L,
        p_l=P_L,
        rho_r=RHO_R,
        u_r=U_R,
        p_r=P_R,
    ):
        grid_point_number = round(((mesh_x_stop - mesh_x_start) / grid_size_x) + 1)
        mesh_x = np.linspace(mesh_x_start, mesh_x_stop, grid_point_number)

        # mesh point number along x
        # it we expect there are total N iteration number,
        # we will need at least N+2 grid points
        mesh_pt_number_x_at_half_t = iteration + 1
        mesh_pt_number_x = iteration + 2

        # u_m of 1D Eular equation.
        # this also means the status of the gas
        # on grids.
        # prefix mtx_ means 'matrix of sth.'
        mtx_q = np.asmatrix(np.zeros(shape=(3, mesh_pt_number_x)))
        # u_m, but means 'next' u_m
        # u_m at the next time step
        mtx_qn = np.asmatrix(np.zeros(shape=(3, mesh_pt_number_x)))
        # partial matrix of q partial x
        mtx_qx = np.asmatrix(np.zeros(shape=(3, mesh_pt_number_x)))
        # partial mastrix of q partial t
        mtx_qt = np.asmatrix(np.zeros(shape=(3, mesh_pt_number_x)))
        # matrix s, it is a part of the marching mtx_qn
        mtx_s = np.asmatrix(np.zeros(shape=(3, mesh_pt_number_x)))
        vxl = np.zeros(shape=(3, 1))
        vxr = np.zeros(shape=(3, 1))

        mtx_f = np.asmatrix(np.zeros(shape=(3, 3)))

        # introduce data object to contain data used during the CESE iteration
        self._data = Data(
            iteration=iteration,
            grid_size_t=grid_size_t,
            grid_size_x=grid_size_x,
            mesh_pt_number_x=mesh_pt_number_x,
            mesh_pt_number_x_at_half_t=mesh_pt_number_x_at_half_t,
            mesh_x=mesh_x,
            mesh_x_start=mesh_x_start,
            mesh_x_stop=mesh_x_stop,
            mtx_f=mtx_f,
            mtx_q=mtx_q,
            mtx_qn=mtx_qn,
            mtx_qx=mtx_qx,
            mtx_qt=mtx_qt,
            mtx_s=mtx_s,
            vxl=vxl,
            vxr=vxr,
        )

        # initialize the gas status before the diaphragm was removed.
        # The number used to calculate the status before
        # the half delta t stepping is applied is 2.
        # This means we begin the iteration process by two grid points.
        self._data.it_pt_nb = 2

        # access necessary data member
        mesh_pt_number_x_at_half_t = self._data.mesh_pt_number_x_at_half_t
        mtx_q = self._data.mtx_q
        # set up the status in the left-hand side
        mtx_q[0][0] = rho_l
        mtx_q[1][0] = rho_l * u_l
        mtx_q[2][0] = p_l / (GAMMA - 1.0) + 0.5 * rho_l * (u_l ** 2.0)
        # set up the status in the right-hand side
        for i in range(mesh_pt_number_x_at_half_t):
            mtx_q[0, i + 1] = rho_r
            mtx_q[1, i + 1] = rho_r * u_r
            mtx_q[2, i + 1] = p_r / (GAMMA - 1.0) + 0.5 * rho_r * (u_r ** 2.0)

        self._data.it_nb = -1  # means CESE iteration ready but not begins

    @property
    def data(self):
        return self._data

    def get_cese_solution(self):
        self._data.refresh_solution()
        return list(self._data.solution)

    def run_cese_iteration(self):
        """
        the whole CESE iteration process
        """
        for i in range(self._data.iteration):
            self._data.it_nb = i
            self.calc_cese_status_before_half_dt()
            self.calc_cese_status_after_half_dt()
            self.push_status_along_t()
        # this is not necessary
        # but it has no risk to refresh and make sure
        # our solution is up-to-date.
        self._data.refresh_solution()

    def calc_cese_status_before_half_dt(self):
        """
        the gas current status

        mtx: means matrix or vector
        """
        data = self._data
        m = data.it_pt_nb
        mtx_q = data.mtx_q
        mtx_f = data.mtx_f
        mtx_qt = data.mtx_qt
        mtx_qx = data.mtx_qx
        mtx_s = data.mtx_s

        dx = self._data.grid_size_x
        dt = self._data.grid_size_t

        for j in range(m):
            w2 = mtx_q[1, j] / mtx_q[0, j]
            w3 = mtx_q[2, j] / mtx_q[0, j]
            # chang95 eq 4.7, the f matrix
            mtx_f[0, 0] = 0.0
            mtx_f[0, 1] = 1.0
            mtx_f[0, 2] = 0.0
            mtx_f[1, 0] = -((3.0 - GAMMA) / 2.0) * (w2 ** 2)
            mtx_f[1, 1] = (3.0 - GAMMA) * w2
            mtx_f[1, 2] = GAMMA - 1.0
            mtx_f[2, 0] = (GAMMA - 1.0) * (w2 ** 3) - GAMMA * w2 * w3
            mtx_f[2, 1] = GAMMA * w3 - 1.5 * (GAMMA - 1.0) * (w2 ** 2)
            mtx_f[2, 2] = GAMMA * w2

            # (4.17) in chang95
            mtx_qt[:, j] = -1.0 * mtx_f * mtx_qx[:, j]
            # (4.25) in chang95, which also uses (4.12) and (4.13)
            mtx_s[:, j] = (
                (dx / 4.0) * mtx_qx[:, j]
                + (dt / dx) * mtx_f * mtx_q[:, j]
                + (dt / dx) * (dt / 4.0) * mtx_f * mtx_qt[:, j]
            )

    def calc_cese_status_after_half_dt(self):
        """
        the gas status after half of dt
        """
        # stepping into the next halt delta t
        # m mesh points along t could introduce m - 1 mesh points along t + 0.5*dt
        data = self._data
        m = data.it_pt_nb
        mm = m - 1
        mtx_q = data.mtx_q
        mtx_qn = data.mtx_qn
        mtx_qt = data.mtx_qt
        mtx_qx = data.mtx_qx
        mtx_s = data.mtx_s
        half_dt = self._data.grid_size_t * 0.5
        half_dx = self._data.grid_size_x * 0.5
        for j in range(mm):
            # (4.24) in chang95
            # Please note the j+1 on the left hand side addresses
            # different t (n) from the j/j+1 on the right hand side,
            # which address the other t (n-1/2) on the space-time
            # surface.
            mtx_qn[:, j + 1] = 0.5 * (
                mtx_q[:, j] + mtx_q[:, j + 1] + mtx_s[:, j] - mtx_s[:, j + 1]
            )
            # (4.27) and (4.36) in chang95
            vxl = np.asarray(
                (mtx_qn[:, j + 1] - mtx_q[:, j] - half_dt * mtx_qt[:, j]) / half_dx
            )
            vxr = np.asarray(
                (mtx_q[:, j + 1] + half_dt * mtx_qt[:, j + 1] - mtx_qn[:, j + 1])
                / half_dx
            )
            # (4.39) in chang95
            #
            # use np.nextafter(0, 1) to provide a very small number
            mtx_qx[:, j + 1] = np.asmatrix(
                (vxl * abs(vxr) + vxr * abs(vxl))
                / (abs(vxl) + abs(vxr) + np.nextafter(0, 1))
            )

    def push_status_along_t(self):
        """
        step into the next iteration status
        """
        # ask the status at t + 0.5*dt to be the next status before the half delta t is applied
        # hdt means 0.5*grid_size_t
        data = self._data
        number_mesh_points_before_hdt = data.it_pt_nb
        mtx_q = data.mtx_q
        mtx_qn = data.mtx_qn

        for j in range(1, number_mesh_points_before_hdt):
            mtx_q[:, j] = mtx_qn[:, j]

        data.it_pt_nb = number_mesh_points_before_hdt + 1


class Data(object):
    """
    a container of the data during the CESE iteration process
    """

    _excludes = [
        "__class__",
        "__delattr__",
        "__dict__",
        "__doc__",
        "__format__",
        "__getattribute__",
        "__hash__",
        "__init__",
        "__module__",
        "__new__",
        "__reduce__",
        "__reduce_ex__",
        "__repr__",
        "__setattr__",
        "__sizeof__",
        "__str__",
        "__subclasshook__",
        "__weakref__",
    ]

    _includes = [
        "iteration",  # total iteration number
        "it_nb",  # the current iteration number, 0 <= it_nb < iteration, -1 means the iteration does not start (yet)
        "it_pt_nb",  # how many x point number to iterate
        "grid_size_t",
        "grid_size_x",
        "mesh_pt_number_x",
        "mesh_pt_number_x_at_half_t",
        "mesh_x",
        "mesh_x_start",
        "mesh_x_stop",
        "mtx_f",
        "mtx_q",
        "mtx_qn",
        "mtx_qx",
        "mtx_qt",
        "mtx_s",
        "vxl",
        "vxr",
        "solution",
    ]

    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            if k in self._excludes or k not in self._includes:
                raise TypeError("{0} is not a valid keyword argument".format(k))
            self.__dict__[k] = v

        self.solution = None

    def refresh_solution(self):
        """
        Use the latest mtx_q to calculate the associated rho, v and p.
        Also, apply physic meaning, e.g. location, according to the given parameters.
        """
        # TODO: warning or implementation for refreshing on odd iteration numbers (status at half t)
        solution = []
        # TODO: only for even iteration number status
        #
        # When self.it_nb is -1, it means iteration does not begin (yet).
        # There are two grid points ready to iterate at the beginning of whole marching,
        # There will be it_nb + 3 iterated grid points during marching. That is to say,
        # it_pt_nb = it_nb + 3
        #   Status              N grid points will be used  Reason
        #   Not begin           2                           2 = -1 +3
        #   After 2 iteration   4                           4 = 1 + 3
        #   After 4 iteration   6                           6 = 3 + 3
        #
        # Every two iteration step make one more x step new derived values in left-hand and right-hand side separately,
        # and we begin the iteration with 2 grid points (it_pt_nb = 2 in the initialized state)
        #
        # In general, if there is len(cese.data.mesh_x)/2 grid points created in the left hand side,
        # we could associate x mesh grid location index with the mtx_q by solution_x_start_index as
        # solution_x_start_index = len(self.mesh_x)/2 - total_iterated_pt_nb/2
        # For example, at the very beginning before marching
        # Assume there are 102 grid points from the beginning,
        # solution_x_start_index = 102/2 - (-1 + 3)/2 = 51 - 1 = 50
        # Because, index counting from 0, there are 51 grid points in the left hand side.
        # So, for 100 iterations, totally 51 grid points are 'used' in the left hand side,
        # and a mesh of at least 102 grid points are needed.

        total_iterated_pt_nb = self.it_nb + 3
        solution_x_start_index = int(
            (len(self.mesh_x) / 2) - (total_iterated_pt_nb / 2)
        )
        for i in range(total_iterated_pt_nb):
            solution_x = self.mesh_x[i + solution_x_start_index]
            solution_rho = self.mtx_q[0, i]
            solution_v = self.mtx_q[1, i] / self.mtx_q[0, i]
            solution_p = (GAMMA - 1.0) * (
                self.mtx_q[2, i] - 0.5 * (solution_v ** 2) * self.mtx_q[0, i]
            )
            # solution format
            solution.append((solution_x, solution_rho, solution_v, solution_p))

        self.solution = self.fill_solution(solution)

    def fill_solution(self, solution_list):
        """
        fill solution array with initial values if no solution represents
        """
        solution_list_rtn = [0] * len(self.mesh_x)
        solution_x_list = []

        for solution in solution_list:
            solution_x_list.append(solution[0])

        for idx in range(len(self.mesh_x)):
            location = self.mesh_x[idx]
            # mapping known solutions
            for solution in solution_list:
                solution_x = solution[0]
                if location == solution_x:
                    solution_list_rtn[idx] = solution
            # fill initial values for the remaining mesh points
            if location not in solution_x_list:
                if location < 0:
                    solution_list_rtn[idx] = (location, RHO_L, U_L, P_L)
                else:
                    solution_list_rtn[idx] = (location, RHO_R, U_R, P_R)

        return solution_list_rtn
