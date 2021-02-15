"""
1D Sod Tube analytic solution solver.

The derivation of the equations for the analytic solution
is based on the book,
Principles of Computational Fluid Dynamics,
written by Pieter Wesseling.
Or, people could refer to the SOLVCON website
https://doc.solvcon.net/en/latest/gas/shocktube.html
"""
import scipy.optimize as so


class ShockTube:
    """
    The core to generate the 1D Sod tube test.
    """

    def __init__(
        self,
        rho_left=1.0,
        u_left=0.0,
        p_left=1.0,
        rho_right=0.125,
        u_right=0.0,
        p_right=0.1,
    ):
        # initial condition
        # [(rhol, ul, pl), (rhor, ur, pr)]
        #
        # Sod's initial condition, L for Left, R for Right
        self.RHOL = rho_left
        self.UL = u_left
        self.PL = p_left
        self.RHOR = rho_right
        self.UR = u_right
        self.PR = p_right

        # constants and conventions
        self.GAMMA = 1.4  # ideal gas constant
        self.GAMMA2 = (self.GAMMA - 1.0) / (self.GAMMA + 1.0)
        self.ALPHA = (self.GAMMA + 1.0) / (self.GAMMA - 1.0)
        self.BETA = (self.GAMMA - 1.0) / (2.0 * self.GAMMA)

    def get_analytic_solution(self, mesh, t=0.2):
        """
        Get analytic solution.

        :param mesh: [point_0, point_1, point_2, point_3, ......, point_n]
        :param t: time, in second
        :return: solution in "solution format", [(x0, rho0, u0, p0), (x1, rho1, u1, p1), ......, (xn, rhon, un, pn)]
        """
        rho4 = self.get_analytic_density_region4()
        u4 = self.get_analytic_velocity_region4()
        p4 = self.get_analytic_pressure_region4()

        rho3 = self.get_analytic_density_region3()
        u3 = self.get_analytic_velocity_region3()
        p3 = self.get_analytic_pressure_region3()

        x_shock = self.get_velocity_shock() * t
        x_disconti = u3 * t
        x_fan_right = self.get_velocity_fan_right() * t
        x_fan_left = self.get_velocity_fan_left() * t

        solution = []
        for x in mesh:
            if x < x_fan_left or x == x_fan_left:
                solution.append(
                    (
                        x,
                        self.get_analytic_density_region1(),
                        self.get_analytic_velocity_region1(),
                        self.get_analytic_pressure_region1(),
                    )
                )
            elif x > x_fan_left and (x < x_fan_right or x == x_fan_right):
                d = self.get_analytic_density_region2(float(x), t)
                v = self.get_analytic_velocity_region2(float(x), t)
                p = self.get_analytic_pressure_region2(float(x), t)
                solution.append((x, d, v, p))
            elif x > x_fan_right and (x < x_disconti or x == x_disconti):
                solution.append((x, rho3, u3, p3))
            elif x > x_disconti and (x < x_shock or x == x_shock):
                solution.append((x, rho4, u4, p4))
            elif x > x_shock:
                solution.append(
                    (
                        x,
                        self.get_analytic_density_region5(),
                        self.get_analytic_velocity_region5(),
                        self.get_analytic_pressure_region5(),
                    )
                )
            else:
                print("Something wrong!!!")

        return solution

    def get_velocity_fan_left(self):
        c1 = self.get_velocity_c1()
        return -c1

    def get_velocity_fan_right(self):
        u3 = self.get_analytic_velocity_region3()
        c3 = self.get_velocity_c3()
        return u3 - c3

    def get_velocity_shock(self):
        # P409, Wesseling P.
        c5 = self.get_velocity_c5()  # 1.0583
        gamma = self.GAMMA
        p4 = self.get_analytic_pressure_region4()  # 0.3031
        p5 = self.get_analytic_pressure_region5()  # 0.1
        return c5 * (
            (1.0 + (((gamma + 1.0) * ((p4 / p5) - 1.0)) / (2.0 * gamma))) ** 0.5
        )

    def get_velocity_c1(self):
        return (self.GAMMA * self.PL / self.RHOL) ** 0.5

    def get_velocity_c3(self):
        p3 = self.get_analytic_pressure_region3()
        rho3 = self.get_analytic_density_region3()
        return (self.GAMMA * p3 / rho3) ** 0.5

    def get_velocity_c5(self):
        return (self.GAMMA * self.PR / self.RHOR) ** 0.5

    def solve_analytic_pressure_region4_by_newton(self, x0=1):
        """
        x0 : the guess initial value to be applied in Newton method
        """
        ######################
        # Analytical formula #
        ######################
        def analytic_pressure_region4(x):
            """
            x: the root value we want to know.

            This method return the formula to get the solution
            of the pressure in the region 4.
            It is a equation that could get the solution
            by numerical approaches, e.g. Newton method.

            For details how to derive the equation, someone
            could refer to, for example, the equation (10.51)
            of Pieter Wesseling,
            Principles of Computational Fluid Dynamics

            The method and the return equation will be
            used by scipy numerial method, e.g.
            scipy.newton
            So, the method and the return value format
            follow the request of scipy.
            """
            p1 = self.PL
            p5 = self.get_analytic_pressure_region5()
            c1 = self.get_velocity_c1()
            c5 = self.get_velocity_c5()
            beta = self.BETA
            gamma = self.GAMMA
            return (x / p1) - (
                (
                    1.0
                    - ((gamma - 1.0) * c5 * ((x / p5) - 1.0))
                    / (
                        c1
                        * (
                            (2.0 * gamma * (gamma - 1.0 + (gamma + 1.0) * (x / p5)))
                            ** 0.5
                        )
                    )
                )
                ** (1.0 / beta)
            )

        return so.newton(analytic_pressure_region4, x0)

    ############
    # Velocity #
    ############
    def get_analytic_velocity_region1(self):
        return self.UL

    def get_analytic_velocity_region2(self, x, t):
        c1 = self.get_velocity_c1()
        gamma = self.GAMMA
        return 2.0 / (gamma + 1.0) * (c1 + x / t)

    def get_analytic_velocity_region3(self):
        return self.get_analytic_velocity_region4()

    def get_analytic_velocity_region4(self):
        """
        The equation could be found in the
        equation next to (10.48), Wesseling P.,
        Principles of Computational Fluid Dynamics
        """
        gamma = self.GAMMA
        p4 = self.get_analytic_pressure_region4()
        p5 = self.get_analytic_pressure_region5()
        p = p4 / p5
        c5 = self.get_velocity_c5()
        return (
            c5 * (p - 1.0) * (2.0 / (gamma * (gamma - 1.0 + (gamma + 1.0) * p))) ** 0.5
        )

    def get_analytic_velocity_region5(self):
        return self.UR

    ############
    # Pressure #
    ############
    def get_analytic_pressure_region1(self):
        return self.PL

    def get_analytic_pressure_region2(self, x, t):
        # (10.44) Wesssling P.
        c1 = self.get_velocity_c1()
        u2 = self.get_analytic_velocity_region2(x, t)
        p1 = self.PL
        gamma = self.GAMMA
        beta = self.BETA
        return p1 * (1.0 - (gamma - 1.0) * u2 / 2 / c1) ** (1.0 / beta)

    def get_analytic_pressure_region3(self):
        return self.get_analytic_pressure_region4()

    def get_analytic_pressure_region4(self):
        return self.solve_analytic_pressure_region4_by_newton()

    def get_analytic_pressure_region5(self):
        return self.PR

    ############
    # Density  #
    ############
    def get_analytic_density_region1(self):
        return self.RHOL

    def get_analytic_density_region2(self, x, t):
        # (10.45), Wesseling P.
        # Principles of Computational Fluid Dynamics
        gamma = self.GAMMA
        rho1 = self.RHOL
        p1 = self.get_analytic_pressure_region1()
        p2 = self.get_analytic_pressure_region2(x, t)
        return rho1 * (p2 / p1) ** (1.0 / gamma)

    def get_analytic_density_region3(self):
        # P410, Wesseling P.
        # Principles of Computational Fluid Dynamics
        rho1 = self.get_analytic_density_region1()
        p1 = self.get_analytic_pressure_region1()
        p3 = self.get_analytic_pressure_region3()
        return rho1 * (p3 / p1) ** (1.0 / self.GAMMA)

    def get_analytic_density_region4(self):
        # P410, Wesseling P.
        # Principles of Computational Fluid Dynamics
        alpha = self.ALPHA
        p4 = self.get_analytic_pressure_region4()
        p5 = self.get_analytic_pressure_region5()
        p = p4 / p5
        rho5 = self.get_analytic_density_region5()
        return rho5 * (1.0 + alpha * p) / (alpha + p)

    def get_analytic_density_region5(self):
        return self.RHOR
