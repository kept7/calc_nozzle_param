import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np

from math import *

from scipy.optimize import *
from scipy.interpolate import PchipInterpolator


from pathlib import Path


# with open(f"D:/1_dev/Part_2_python/projects/test/2023.11_november/calc_of_nozzle_parameters/initial_data.txt", "r") as initial_data:
#     str_of_data = str(initial_data.read())

# list_of_data = str_of_data.split("\n")


# def creating_list_of_pressure(list_of_data):
#     str_of_pressure = str(list_of_data[0])
#     list_of_pressure = str_of_pressure.split(" ")
#     list_of_pressure.pop(0)
#     list_of_pressure.pop(0)
#     list_of_pressure = [float(el) for el in list_of_pressure]
#     return(list_of_pressure)


# user must enter the following variables: given_thrust, specific_impulse, double_apostrophe_critical_section_area,
# double_apostrophe_nozzle_exit_section_area, ratio_of_nozzle_expender_cross_area_and_nozzle_exit_section_area,
# alpha_enter, length_of_stay_of_fuel_in_comb_cham, time_of_stay_of_fuel_in_comb_cham =>
# => needs to end the task with adding the func of transform information from user's notebook


def enter_data():
    pressure_in_combustion_chamber = 15
    temperature_in_combustion_chamber = 3297.05
    density_in_combustion_chamber = 13.331

    specific_gas_constant = 341.28
    frex_specific_gas_constant = 341.28

    given_thrust = 1095000
    specific_impulse = 3198.42

    k = 1.1986
    k_frex = 1.2

    specific_critical_section_area = 0.0001112 * (10**6)
    specific_nozzle_exit_section_area = 0.0034731 * (10**6)
    ratio_of_nozzle_expender_cross_area_and_nozzle_exit_section_area = 0.5

    enter_alpha = 50
    beta = 15
    length_of_stay_of_fuel_in_comb_chamber = 1500

    specific_section_area_of_gas = [0.0004448, 0.0002801, 0.0002027, 0.0001696, 0.0001506, 
                                    0.0001382, 0.0001297, 0.0001235, 0.0001191, 0.0001158, 
                                    0.0001136, 0.0001121, 0.0001114, 0.0001112, 0.0001117, 
                                    0.0001129, 0.0001146, 0.0001170, 0.0001203, 0.0001244, 
                                    0.0001297, 0.0001363, 0.0001448, 0.0001558, 0.0001704, 
                                    0.0001906, 0.0002201, 0.0002672, 0.0003554, 0.0005905, 
                                    0.0034731,
    ]

    mach_of_gas_in_sections = [0.00, 0.244439, 0.349113, 0.431928, 0.503973, 0.569567,
                               0.630909, 0.689375, 0.745834, 0.800965, 0.855284, 0.909213,
                               0.963119, 1, 1.07227, 1.12818, 1.18549, 1.24459, 1.30596, 
                               1.37013, 1.43784, 1.50997, 1.58767, 1.67256, 1.76695,
                               1.87438, 2.00071, 2.15679, 2.36665, 2.70559, 3.86319
    ]

    mach_of_frex_gas_in_sections = [0.00, 0.2446, 0.3495, 0.4327, 0.5052, 0.5714,
                                    0.6334, 0.6926, 0.7499, 0.8060, 0.8615, 0.9167,
                                    0.9720, 1.00, 1.0846, 1.1427, 1.2024, 1.2643,
                                    1.3288, 1.3967, 1.4688, 1.5461, 1.6301, 1.7226,
                                    1.8265, 1.9463, 2.0892, 2.2693, 2.5182, 2.9388,
                                    4.4378
    ]


    return (
        pressure_in_combustion_chamber,
        temperature_in_combustion_chamber,
        density_in_combustion_chamber,
        specific_gas_constant,
        frex_specific_gas_constant,
        given_thrust,
        specific_impulse,
        k,
        k_frex,
        specific_critical_section_area,
        specific_nozzle_exit_section_area,
        ratio_of_nozzle_expender_cross_area_and_nozzle_exit_section_area,
        enter_alpha,
        beta,
        length_of_stay_of_fuel_in_comb_chamber,
        specific_section_area_of_gas,
        mach_of_gas_in_sections,
        mach_of_frex_gas_in_sections,
    )


def mass_flow_rate_calculation(given_thrust, specific_impulse):
    mass_flow_rate = given_thrust / specific_impulse

    return mass_flow_rate


def area_and_diametrs_calculation(
    mass_flow_rate,
    specific_critical_section_area,
    specific_nozzle_exit_section_area,
    ratio_of_nozzle_expender_cross_area_and_nozzle_exit_section_area,
    specific_section_area_of_gas,
):
    # critical section calculation
    critical_section_area = mass_flow_rate * specific_critical_section_area
    diametr_of_critical_section = sqrt(4 * critical_section_area / np.pi)
    radius_of_critical_section = diametr_of_critical_section / 2

    # calculation of the cross section at the nozzle exit
    nozzle_exit_section_area = (
        mass_flow_rate * specific_nozzle_exit_section_area
    )
    diametr_of_nozzle_exit_section = sqrt(4 * nozzle_exit_section_area / np.pi)
    radius_of_exit_section = diametr_of_nozzle_exit_section // 2

    # calculation of the nozzle expender cross section
    if ratio_of_nozzle_expender_cross_area_and_nozzle_exit_section_area > 1:
        nozzle_expender_cross_area = (
            nozzle_exit_section_area
            * ratio_of_nozzle_expender_cross_area_and_nozzle_exit_section_area
        )
        diametr_of_nozzle_expender_cross_section = sqrt(
            4 * nozzle_expender_cross_area / np.pi
        )
        radius_of_nozzle_expender_cross_section = (
            diametr_of_nozzle_expender_cross_section // 2
        )
    else:
        radius_of_nozzle_expender_cross_section = 0

    # calculation of the combustuion chamber
    combustion_chamber_area = 4 * critical_section_area
    diametr_of_combustion_chamber = sqrt(4 * combustion_chamber_area / np.pi)
    radius_of_combustion_chamber = diametr_of_combustion_chamber / 2

    # calcualtion of sections areas in which we will consider
    for el in range(len(specific_section_area_of_gas)):
        specific_section_area_of_gas[el] = (
            mass_flow_rate * specific_section_area_of_gas[el]
        )

    section_area_of_gas = specific_section_area_of_gas

    return (
        critical_section_area,
        combustion_chamber_area,
        radius_of_combustion_chamber,
        radius_of_critical_section,
        radius_of_exit_section,
        radius_of_nozzle_expender_cross_section,
        section_area_of_gas,
    )


def calculation_of_parametrs_up_to_the_critical_section(
    radius_of_combustion_chamber,
    radius_of_critical_section,
    enter_alpha,
    length_of_stay_of_fuel_in_comb_chamber,
    critical_section_area,
    combustion_chamber_area,
):
    first_radius_of_tapering_part = 0.8 * radius_of_combustion_chamber
    second_radius_of_tapering_part = 0.8 * radius_of_critical_section

    ya = radius_of_combustion_chamber
    yc1 = ya - first_radius_of_tapering_part
    yc = yc1 + first_radius_of_tapering_part * sin(radians(90 - enter_alpha))

    ye = radius_of_critical_section + second_radius_of_tapering_part
    yd = ye - second_radius_of_tapering_part * sin(radians(90 - enter_alpha))

    # delta_x = xd - xc
    delta_x = (yc - yd) / tan(radians(enter_alpha))
    length_of_the_tapering_part = (
        first_radius_of_tapering_part * sin(radians(enter_alpha))
        + second_radius_of_tapering_part * sin(radians(enter_alpha))
        + delta_x
    )
    volume_of_the_tapering_part = (
        pi
        / 3
        * length_of_the_tapering_part
        * (
            radius_of_combustion_chamber**2
            + radius_of_combustion_chamber * radius_of_critical_section
            + radius_of_critical_section**2
        )
    )
    length_of_combustion_chamber = (
        length_of_stay_of_fuel_in_comb_chamber * critical_section_area
        + volume_of_the_tapering_part
    ) / combustion_chamber_area

    return (
        first_radius_of_tapering_part,
        second_radius_of_tapering_part,
        length_of_combustion_chamber,
    )


def coordinates_for_nozzle_graph(
    enter_alpha,
    radius_of_critical_section,
    radius_of_combustion_chamber,
    radius_of_nozzle_expender_cross_section,
    first_radius_of_tapering_part,
    second_radius_of_tapering_part,
    length_of_combustion_chamber,
):
    x0 = 0
    y0 = radius_of_combustion_chamber

    xa = length_of_combustion_chamber
    ya = y0

    xc1 = xa
    yc1 = ya - first_radius_of_tapering_part

    xc = xa + first_radius_of_tapering_part * cos(radians(90 - enter_alpha))
    yc = yc1 + first_radius_of_tapering_part * sin(radians(90 - enter_alpha))

    ye = radius_of_critical_section + second_radius_of_tapering_part

    yd = ye - second_radius_of_tapering_part * sin(radians(90 - enter_alpha))
    xd = xc + (yc - yd) / tan(radians(enter_alpha))

    xc2 = xd + second_radius_of_tapering_part * sin(radians(enter_alpha))

    xe = xc2

    xf = xc2

    if radius_of_nozzle_expender_cross_section != 0:
        yg0 = radius_of_exit_section
        xg0 = (yg0 - radius_of_critical_section) / tan(radians(beta)) + xf
        yg = radius_of_nozzle_expender_cross_section
        xg = (yg - radius_of_critical_section) / tan(radians(beta)) + xf
    else:
        yg = radius_of_exit_section
        xg = (yg - radius_of_critical_section) / tan(radians(beta)) + xf
        yg0 = yg
        xg0 = xg

    def creation_of_upper_graph_of_nozzle(x):
        if x0 == x:
            return 0
        if x0 < x <= xa:
            return ya
        elif xa < x <= xc:
            return yc1 + sqrt(
                first_radius_of_tapering_part**2 - (x - xc1) ** 2
            )
        elif xc < x <= xd:
            return -tan(radians(enter_alpha)) * (x - xc) + yc
        elif xd < x <= xf:
            return ye - sqrt(
                second_radius_of_tapering_part**2 - (x - xe) ** 2
            )
        elif xf < x <= xg:
            # calculation of y for conical part of the nozzle with beta -> replace with profiled nozzle
            return radius_of_critical_section + tan(radians(beta)) * (x - xf)
        elif x > xg:
            return 0

    y_for_upper_nozzle_graph = [
        creation_of_upper_graph_of_nozzle(x + xe)
        for x in np.arange(-xe, xg - xe + 0.01, 0.01)
    ]
    x_for_upper_nozzle_graph = np.arange(-xe, xg - xe + 0.01, 0.01)

    def creation_of_lower_graph_of_nozzle(x):
        if x0 == x:
            return 0
        if x0 < x <= xa:
            return -ya
        elif xa < x <= xc:
            return -(
                yc1 + sqrt(first_radius_of_tapering_part**2 - (x - xc1) ** 2)
            )
        elif xc < x <= xd:
            return -(-tan(radians(enter_alpha)) * (x - xc) + yc)
        elif xd < x <= xf:
            return -(
                ye - sqrt(second_radius_of_tapering_part**2 - (x - xe) ** 2)
            )
        elif xf < x <= xg:
            # calculation of y for conical part of the nozzle with beta -> replace with profiled nozzle
            return -(radius_of_critical_section + tan(radians(beta)) * (x - xf))
        elif x > xg:
            return 0

    y_for_lower_nozzle_graph = [
        creation_of_lower_graph_of_nozzle(x + xe)
        for x in np.arange(-xe, xg - xe + 0.01, 0.01)
    ]
    x_for_lower_nozzle_graph = x_for_upper_nozzle_graph

    return (
        x0,
        xa,
        ya,
        xc,
        yc,
        xd,
        yd,
        ye,
        xf,
        xe,
        xc1,
        yc1,
        xg,
        yg,
        x_for_upper_nozzle_graph,
        y_for_upper_nozzle_graph,
        x_for_lower_nozzle_graph,
        y_for_lower_nozzle_graph,
        xg0,
        yg0,
    )


def nozzle_graph_creation(
    x_for_upper_nozzle_graph,
    y_for_upper_nozzle_graph,
    x_for_lower_nozzle_graph,
    y_for_lower_nozzle_graph,
    xg,
    yg,
    xg0,
    yg0,
):
    plt.figure()
    plt.title("Nozzle Profile")
    plt.ylabel("$y(x_i)$ - Nozzle height y, mm")
    plt.xlabel("$x_i$ - Coordinate along the nozzle length, mm")

    plt.axis([-6 * xe / 5, 4 * xg / 5, -2 * yg, 2 * yg])

    if radius_of_nozzle_expender_cross_section == 0:
        plt.plot(
            x_for_upper_nozzle_graph,
            y_for_upper_nozzle_graph,
            x_for_lower_nozzle_graph,
            y_for_lower_nozzle_graph,
            color="black",
        )
    else:
        y_cross = [yg0, 0, -yg0]
        x_cross = [xg0 - xe, xg0 - xe, xg0 - xe]
        plt.plot(
            x_for_upper_nozzle_graph,
            y_for_upper_nozzle_graph,
            x_for_lower_nozzle_graph,
            y_for_lower_nozzle_graph,
            color="black",
        )
        plt.plot(x_cross, y_cross, linestyle="dashed", color="darkgray")

    plt.grid(True)

    plt.show()


def calculation_of_x_coordinate_for_gas_parameters(
    section_area_of_gas,
    first_radius_of_tapering_part,
    second_radius_of_tapering_part,
    enter_alpha,
):
    y_coordinates_for_gas_parameters_in_rocket_engine = section_area_of_gas

    for el in range(len(section_area_of_gas)):
        y_coordinates_for_gas_parameters_in_rocket_engine[el] = 500 * sqrt(
            (4 * y_coordinates_for_gas_parameters_in_rocket_engine[el]) / pi
        )

    list_of_y_coordinates_before_critical_section = (
        y_coordinates_for_gas_parameters_in_rocket_engine[
            : (
                y_coordinates_for_gas_parameters_in_rocket_engine.index(
                    min(y_coordinates_for_gas_parameters_in_rocket_engine)
                )
                + 1
            )
        ]
    )
    list_of_y_coordinates_after_critical_section = (
        y_coordinates_for_gas_parameters_in_rocket_engine[
            (
                y_coordinates_for_gas_parameters_in_rocket_engine.index(
                    min(y_coordinates_for_gas_parameters_in_rocket_engine)
                )
                + 1
            ) :
        ]
    )

    def x_coordinate_before_critical_section(y):
        if ya >= y >= yc:
            return round(
                (
                    sqrt(first_radius_of_tapering_part**2 - (y - yc1) ** 2)
                    + xc1
                    - xe
                ),
                2,
            )
        if yc >= y >= yd:
            x = fsolve(
                lambda x: (-tan(radians(enter_alpha)) * (x + xe - xc) + yc) - y,
                0.5,
            )
            return round((x[0]), 2)
        if yd >= y >= radius_of_critical_section:
            return round(
                (sqrt(second_radius_of_tapering_part**2 - (y - ye) ** 2)), 2
            )
        if y <= radius_of_critical_section:
            return 0

    def x_coordinate_after_critical_section(y):
        if radius_of_critical_section <= y:
            x = fsolve(
                lambda x: (
                    radius_of_critical_section
                    + tan(radians(beta)) * (x + xe - xf)
                )
                - y,
                1.5,
            )
            return round((x[0]), 2)

    list_of_x_coordinates_before_critical_section = [
        x_coordinate_before_critical_section(el)
        for el in list_of_y_coordinates_before_critical_section
    ]
    list_of_x_coordinates_after_critical_section = [
        x_coordinate_after_critical_section(el)
        for el in list_of_y_coordinates_after_critical_section
    ]

    for el in range(len(list_of_x_coordinates_before_critical_section)):
        if list_of_x_coordinates_before_critical_section[el] >= 0:
            list_of_x_coordinates_before_critical_section[el] = (
                list_of_x_coordinates_before_critical_section[el] * -1
            )
        else:
            list_of_x_coordinates_before_critical_section[
                el
            ] = list_of_x_coordinates_before_critical_section[el]

    x_coordinate_for_gas_parameters = (
        list_of_x_coordinates_before_critical_section
        + list_of_x_coordinates_after_critical_section
    )

    return x_coordinate_for_gas_parameters


def calculation_of_pressure_along_of_the_nozzle(
    pressure_in_combustion_chamber, k, mach_of_gas_in_sections
):
    def y_coordinate_for_pressure_along_of_the_nozzle(M):
        p_value_in_section = fsolve(
            lambda p_value_in_section: (
                (p_value_in_section / pressure_in_combustion_chamber)
                - (1 / (1 + ((k - 1) / 2)* M**2)) **(k / (k - 1))
            ),
            1.5,
        )
        return round(p_value_in_section[0], 5)

    y_coordinate_for_pressure_graph = [
        y_coordinate_for_pressure_along_of_the_nozzle(M)
        for M in mach_of_gas_in_sections
    ]
    return y_coordinate_for_pressure_graph


def pressure_graph_creation(
    x_coordinate_for_gas_parameters,
    y_coordinate_for_pressure_graph,
    y_coordinate_for_frex_pressure_graph,
):
    plt.figure()
    plt.title("Pressure distribution along the nozzle length")
    plt.ylabel(r"$p(x_i)$ - Pressure, MPa")
    plt.xlabel(r"$x_i$ - Coordinate along the nozzle length, mm")

    plt.axis(
        [
            x_coordinate_for_gas_parameters[0] - 50,
            x_coordinate_for_gas_parameters[-1] + 50,
            -1,
            16 * y_coordinate_for_pressure_graph[0] / 15,
        ]
    )

    if radius_of_nozzle_expender_cross_section == 0:
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_pressure_graph,
            color="red",
        )
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_frex_pressure_graph,
            color="steelblue",
        )
    else:
        y_cross = [yg0, 0, -yg0]
        x_cross = [xg0 - xe, xg0 - xe, xg0 - xe]
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_pressure_graph,
            color="red",
        )
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_frex_pressure_graph,
            color="steelblue",
        )
        plt.plot(x_cross, y_cross, linestyle="dashed", color="darkgray")

    plt.grid(True)

    return plt.show()


def calculation_of_temperature_along_of_the_nozzle(
    temperature_in_combustion_chamber, k, mach_of_gas_in_sections
):
    def y_coordinate_for_temperature_along_of_the_nozzle(M):
        t_value_in_section = fsolve(
            lambda t_value_in_section: (
                t_value_in_section / temperature_in_combustion_chamber
            )
            - 1 / (1 + (k - 1) / 2 * M**2),
            1.5,
        )
        return round(t_value_in_section[0], 5)

    y_coordinate_for_temperature_graph = [
        y_coordinate_for_temperature_along_of_the_nozzle(M)
        for M in mach_of_gas_in_sections
    ]

    return y_coordinate_for_temperature_graph


def temperature_graph_creation(
    x_coordinate_for_gas_parameters,
    y_coordinate_for_temperature_graph,
    y_coordinate_for_frex_temperature_graph,
):
    plt.figure()
    plt.title("Temperature distribution along the nozzle length")
    plt.ylabel(r"$T(x_i)$ - Temperature, K")
    plt.xlabel(r"$x_i$ - Coordinate along the nozzle length, mm")

    plt.axis(
        [
            x_coordinate_for_gas_parameters[0] - 50,
            x_coordinate_for_gas_parameters[-1] + 50,
            y_coordinate_for_frex_temperature_graph[-1] - 100,
            16 * y_coordinate_for_temperature_graph[0] / 15,
        ]
    )

    ys = PchipInterpolator(
        x_coordinate_for_gas_parameters, y_coordinate_for_temperature_graph
    )
    xs = np.arange(
        x_coordinate_for_gas_parameters[0],
        x_coordinate_for_gas_parameters[-1],
        0.1,
    )

    if radius_of_nozzle_expender_cross_section == 0:
        # plt.plot(xs, ys(xs),  color='red')
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_temperature_graph,
            color="red",
        )
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_frex_temperature_graph,
            color="steelblue",
        )
    else:
        y_cross = [yg0, 0, -yg0]
        x_cross = [xg0 - xe, xg0 - xe, xg0 - xe]
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_temperature_graph,
            color="red",
        )
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_frex_temperature_graph,
            color="steelblue",
        )
        plt.plot(x_cross, y_cross, linestyle="dashed", color="darkgray")

    plt.grid(True)

    return plt.show()


def calculation_of_density_along_of_the_nozzle(
    density_in_combustion_chamber, k, mach_of_gas_in_sections
):
    def y_coordinate_for_density_along_of_the_nozzle(M):
        density_value_in_section = fsolve(
            lambda density_value_in_section: (
                density_value_in_section / density_in_combustion_chamber
            )
            - (1 / (1 + (k - 1) / 2 * M**2))
            ** (1 / (k - 1)),
            1.5,
        )
        return round(density_value_in_section[0], 5)

    y_coordinate_for_density_graph = [
        y_coordinate_for_density_along_of_the_nozzle(M)
        for M in mach_of_gas_in_sections
    ]

    return y_coordinate_for_density_graph


def density_graph_creation(
    x_coordinate_for_gas_parameters,
    y_coordinate_for_density_graph,
    y_coordinate_for_frex_density_graph,
):
    plt.figure()
    plt.title("Density distribution along the nozzle length")
    plt.ylabel(r"$\rho(x_i)$ - Density, kg/m^3")
    plt.xlabel(r"$x_i$ - Coordinate along the nozzle length, mm")

    plt.axis(
        [
            x_coordinate_for_gas_parameters[0] - 50,
            x_coordinate_for_gas_parameters[-1] + 50,
            -1,
            16 * y_coordinate_for_density_graph[0] / 15,
        ]
    )

    if radius_of_nozzle_expender_cross_section == 0:
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_density_graph,
            color="red",
        )
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_frex_density_graph,
            color="steelblue",
        )
    else:
        y_cross = [yg0, 0, -yg0]
        x_cross = [xg0 - xe, xg0 - xe, xg0 - xe]
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_density_graph,
            color="red",
        )
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_frex_density_graph,
            color="steelblue",
        )
        plt.plot(x_cross, y_cross, linestyle="dashed", color="darkgray")

    plt.grid(True)

    return plt.show()


def calculation_of_velocity_along_of_the_nozzle(
    specific_gas_constant,
    k,
    mach_of_gas_in_sections,
    y_coordinate_for_temperature_graph,
):
    def kRT(M):
        return (
            k
            * specific_gas_constant
            * y_coordinate_for_temperature_graph[
                mach_of_gas_in_sections.index(M)
            ]
        )

    def y_coordinate_for_velocity_along_of_the_nozzle(M):
        w_value_in_section = fsolve(
            lambda w_value_in_section: w_value_in_section - M * sqrt(kRT(M)),
            1.5,
        )
        return round(w_value_in_section[0], 5)

    y_coordinate_for_velocity_graph = [
        y_coordinate_for_velocity_along_of_the_nozzle(M)
        for M in mach_of_gas_in_sections
    ]

    return y_coordinate_for_velocity_graph


def velocity_graph_creation(
    x_coordinate_for_gas_parameters,
    y_coordinate_for_velocity_graph,
    y_coordinate_for_frex_velocity_graph,
):
    plt.figure()
    plt.title("Velocity distribution along the nozzle length")
    plt.ylabel("$W(x_i)$ - Velocity, m/s")
    plt.xlabel(r"$x_i$ - Coordinate along the nozzle length, mm")

    plt.axis(
        [
            x_coordinate_for_gas_parameters[0] - 50,
            x_coordinate_for_gas_parameters[-1] + 50,
            y_coordinate_for_velocity_graph[0] - 100,
            16 * y_coordinate_for_velocity_graph[-1] / 15,
        ]
    )

    if radius_of_nozzle_expender_cross_section == 0:
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_velocity_graph,
            color="red",
        )
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_frex_velocity_graph,
            color="steelblue",
        )
    else:
        y_cross = [yg0, 0, -yg0]
        x_cross = [xg0 - xe, xg0 - xe, xg0 - xe]
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_velocity_graph,
            color="black",
        )
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_frex_velocity_graph,
            color="steelblue",
        )
        plt.plot(x_cross, y_cross, linestyle="dashed", color="darkgray")

    plt.grid(True)

    return plt.show()


def mach_graph_creation(
    x_coordinate_for_gas_parameters,
    mach_of_gas_in_sections,
    mach_of_frex_gas_in_sections,
):
    plt.figure()
    plt.title("Distribution of Mach number along the nozzle length")
    plt.ylabel("$M(x_i)$ - Mach number")
    plt.xlabel(r"$x_i$ - Coordinate along the nozzle length, mm")

    plt.axis(
        [
            x_coordinate_for_gas_parameters[0] - 50,
            x_coordinate_for_gas_parameters[-1] + 50,
            mach_of_gas_in_sections[0] - 0.25,
            round(mach_of_frex_gas_in_sections[-1], 1) + 0.25,
        ]
    )

    if radius_of_nozzle_expender_cross_section == 0:
        plt.plot(
            x_coordinate_for_gas_parameters,
            mach_of_gas_in_sections,
            color="red",
        )
        plt.plot(
            x_coordinate_for_gas_parameters,
            mach_of_frex_gas_in_sections,
            color="steelblue",
        )
    else:
        y_cross = [yg0, 0, -yg0]
        x_cross = [xg0 - xe, xg0 - xe, xg0 - xe]
        plt.plot(
            x_coordinate_for_gas_parameters,
            mach_of_gas_in_sections,
            color="black",
        )
        plt.plot(
            x_coordinate_for_gas_parameters,
            mach_of_frex_gas_in_sections,
            color="steelblue",
        )
        plt.plot(x_cross, y_cross, linestyle="dashed", color="darkgray")

    plt.grid(True)

    return plt.show()


def calculation_of_lambda_along_of_the_nozzle(k, mach_of_gas_in_sections):
    def y_coordinate_for_lambda_along_of_the_nozzle(M):
        lambda_value_in_section = fsolve(
            lambda lambda_value_in_section: (
                M * sqrt((k + 1) / 2)
            )
            / sqrt(1 + (k - 1) / 2 * M**2)
            - lambda_value_in_section,
            1.5,
        )
        return round(lambda_value_in_section[0], 5)

    y_coordinate_for_lambda_graph = [
        y_coordinate_for_lambda_along_of_the_nozzle(M)
        for M in mach_of_gas_in_sections
    ]

    return y_coordinate_for_lambda_graph


def lambda_graph_creation(
    x_coordinate_for_gas_parameters,
    y_coordinate_for_lambda_graph,
    y_coordinate_for_frex_lambda_graph,
):
    plt.figure()
    plt.title("Distribution of reduced velocity along the nozzle length")
    plt.ylabel(r"$\lambda(x_i)$ - Reduced speed")
    plt.xlabel(r"$x_i$ - coordinate along the nozzle length, mm")

    plt.axis(
        [
            x_coordinate_for_gas_parameters[0] - 50,
            x_coordinate_for_gas_parameters[-1] + 50,
            y_coordinate_for_lambda_graph[0] - 0.25,
            round(y_coordinate_for_lambda_graph[-1], 0),
        ]
    )

    if radius_of_nozzle_expender_cross_section == 0:
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_lambda_graph,
            color="red",
        )
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_frex_lambda_graph,
            color="steelblue",
        )
    else:
        y_cross = [yg0, 0, -yg0]
        x_cross = [xg0 - xe, xg0 - xe, xg0 - xe]
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_lambda_graph,
            color="black",
        )
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_frex_lambda_graph,
            color="steelblue",
        )
        plt.plot(x_cross, y_cross, linestyle="dashed", color="darkgray")

    plt.grid(True)

    return plt.show()


def calculation_of_multiplication_rho_and_velocity_along_of_the_nozzle(
    y_coordinate_for_velocity_graph, y_coordinate_for_density_graph
):
    y_coordinate_for_multiplication_rho_and_velocity_graph = []

    for el in range(len(y_coordinate_for_velocity_graph)):
        y_coordinate_for_velocity_graph[el] = (
            y_coordinate_for_velocity_graph[el]
            * y_coordinate_for_density_graph[el]
        )

    y_coordinate_for_multiplication_rho_and_velocity_graph = (
        y_coordinate_for_velocity_graph
    )

    return y_coordinate_for_multiplication_rho_and_velocity_graph


def multiplication_rho_and_velocity_graph_creation(
    x_coordinate_for_gas_parameters,
    y_coordinate_for_multiplication_rho_and_velocity_graph,
    y_coordinate_for_frex_multiplication_rho_and_velocity_graph,
):
    plt.figure()
    plt.title("Distribution of reduced velocity along the nozzle length")
    plt.ylabel(
        r"$\rho(x_i) * W(x_i)$ - Density multiplied by velocity, kg/(c*m^2)"
    )
    plt.xlabel(r"$x_i$ - Coordinate along the nozzle length, mm")

    plt.axis(
        [
            x_coordinate_for_gas_parameters[0] - 50,
            x_coordinate_for_gas_parameters[-1] + 50,
            y_coordinate_for_multiplication_rho_and_velocity_graph[0] - 100,
            round(y_coordinate_for_multiplication_rho_and_velocity_graph[14], 0)
            + 250,
        ]
    )

    if radius_of_nozzle_expender_cross_section == 0:
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_multiplication_rho_and_velocity_graph,
            color="red",
        )
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_frex_multiplication_rho_and_velocity_graph,
            color="steelblue",
        )
    else:
        y_cross = [yg0, 0, -yg0]
        x_cross = [xg0 - xe, xg0 - xe, xg0 - xe]
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_multiplication_rho_and_velocity_graph,
            color="black",
        )
        plt.plot(
            x_coordinate_for_gas_parameters,
            y_coordinate_for_frex_multiplication_rho_and_velocity_graph,
            color="steelblue",
        )
        plt.plot(x_cross, y_cross, linestyle="dashed", color="darkgray")

    plt.grid(True)

    return plt.show()


if __name__ == "__main__":
    (
        pressure_in_combustion_chamber,
        temperature_in_combustion_chamber,
        density_in_combustion_chamber,
        specific_gas_constant,
        frex_specific_gas_constant,
        given_thrust,
        specific_impulse,
        k,
        k_frex,
        specific_critical_section_area,
        specific_nozzle_exit_section_area,
        ratio_of_nozzle_expender_cross_area_and_nozzle_exit_section_area,
        enter_alpha,
        beta,
        length_of_stay_of_fuel_in_comb_chamber,
        specific_section_area_of_gas,
        mach_of_gas_in_sections,
        mach_of_frex_gas_in_sections,
    ) = enter_data()
    mass_flow_rate = mass_flow_rate_calculation(given_thrust, specific_impulse)

    (
        critical_section_area,
        combustion_chamber_area,
        radius_of_combustion_chamber,
        radius_of_critical_section,
        radius_of_exit_section,
        radius_of_nozzle_expender_cross_section,
        section_area_of_gas,
    ) = area_and_diametrs_calculation(
        mass_flow_rate,
        specific_critical_section_area,
        specific_nozzle_exit_section_area,
        ratio_of_nozzle_expender_cross_area_and_nozzle_exit_section_area,
        specific_section_area_of_gas,
    )
    (
        first_radius_of_tapering_part,
        second_radius_of_tapering_part,
        length_of_combustion_chamber,
    ) = calculation_of_parametrs_up_to_the_critical_section(
        radius_of_combustion_chamber,
        radius_of_critical_section,
        enter_alpha,
        length_of_stay_of_fuel_in_comb_chamber,
        critical_section_area,
        combustion_chamber_area,
    )

    (
        x0,
        xa,
        ya,
        xc,
        yc,
        xd,
        yd,
        ye,
        xf,
        xe,
        xc1,
        yc1,
        xg,
        yg,
        x_for_upper_nozzle_graph,
        y_for_upper_nozzle_graph,
        x_for_lower_nozzle_graph,
        y_for_lower_nozzle_graph,
        xg0,
        yg0,
    ) = coordinates_for_nozzle_graph(
        enter_alpha,
        radius_of_critical_section,
        radius_of_combustion_chamber,
        radius_of_nozzle_expender_cross_section,
        first_radius_of_tapering_part,
        second_radius_of_tapering_part,
        length_of_combustion_chamber,
    )
    nozzle_graph_creation(
        x_for_upper_nozzle_graph,
        y_for_upper_nozzle_graph,
        x_for_lower_nozzle_graph,
        y_for_lower_nozzle_graph,
        xg,
        yg,
        xg0,
        yg0,
    )

    x_coordinate_for_gas_parameters = (
        calculation_of_x_coordinate_for_gas_parameters(
            section_area_of_gas,
            first_radius_of_tapering_part,
            second_radius_of_tapering_part,
            enter_alpha,
        )
    )

    y_coordinate_for_pressure_graph = (
        calculation_of_pressure_along_of_the_nozzle(
            pressure_in_combustion_chamber, k, mach_of_gas_in_sections
        )
    )
    y_coordinate_for_frex_pressure_graph = (
        calculation_of_pressure_along_of_the_nozzle(
            pressure_in_combustion_chamber, k_frex, mach_of_frex_gas_in_sections
        )
    )
    pressure_graph_creation(
        x_coordinate_for_gas_parameters,
        y_coordinate_for_pressure_graph,
        y_coordinate_for_frex_pressure_graph,
    )

    y_coordinate_for_temperature_graph = (
        calculation_of_temperature_along_of_the_nozzle(
            temperature_in_combustion_chamber, k, mach_of_gas_in_sections
        )
    )
    y_coordinate_for_frex_temperature_graph = (
        calculation_of_temperature_along_of_the_nozzle(
            temperature_in_combustion_chamber,
            k_frex,
            mach_of_frex_gas_in_sections,
        )
    )
    temperature_graph_creation(
        x_coordinate_for_gas_parameters,
        y_coordinate_for_temperature_graph,
        y_coordinate_for_frex_temperature_graph,
    )

    y_coordinate_for_density_graph = calculation_of_density_along_of_the_nozzle(
        density_in_combustion_chamber, k, mach_of_gas_in_sections
    )
    y_coordinate_for_frex_density_graph = (
        calculation_of_density_along_of_the_nozzle(
            density_in_combustion_chamber, k, mach_of_frex_gas_in_sections
        )
    )
    density_graph_creation(
        x_coordinate_for_gas_parameters,
        y_coordinate_for_density_graph,
        y_coordinate_for_frex_density_graph,
    )

    y_coordinate_for_velocity_graph = (
        calculation_of_velocity_along_of_the_nozzle(
            specific_gas_constant,
            k,
            mach_of_gas_in_sections,
            y_coordinate_for_temperature_graph,
        )
    )
    y_coordinate_for_frex_velocity_graph = (
        calculation_of_velocity_along_of_the_nozzle(
            frex_specific_gas_constant,
            k_frex,
            mach_of_frex_gas_in_sections,
            y_coordinate_for_frex_temperature_graph,
        )
    )
    velocity_graph_creation(
        x_coordinate_for_gas_parameters,
        y_coordinate_for_velocity_graph,
        y_coordinate_for_frex_velocity_graph,
    )

    mach_graph_creation(
        x_coordinate_for_gas_parameters,
        mach_of_gas_in_sections,
        mach_of_frex_gas_in_sections,
    )

    y_coordinate_for_lambda_graph = calculation_of_lambda_along_of_the_nozzle(
        k, mach_of_gas_in_sections
    )
    y_coordinate_for_frex_lambda_graph = (
        calculation_of_lambda_along_of_the_nozzle(
            k_frex, mach_of_frex_gas_in_sections
        )
    )
    lambda_graph_creation(
        x_coordinate_for_gas_parameters,
        y_coordinate_for_lambda_graph,
        y_coordinate_for_frex_lambda_graph,
    )

    y_coordinate_for_multiplication_rho_and_velocity_graph = (
        calculation_of_multiplication_rho_and_velocity_along_of_the_nozzle(
            y_coordinate_for_velocity_graph, y_coordinate_for_density_graph
        )
    )
    y_coordinate_for_frex_multiplication_rho_and_velocity_graph = (
        calculation_of_multiplication_rho_and_velocity_along_of_the_nozzle(
            y_coordinate_for_frex_velocity_graph,
            y_coordinate_for_frex_density_graph,
        )
    )
    multiplication_rho_and_velocity_graph_creation(
        x_coordinate_for_gas_parameters,
        y_coordinate_for_multiplication_rho_and_velocity_graph,
        y_coordinate_for_frex_multiplication_rho_and_velocity_graph,
    )

    print(y_coordinate_for_pressure_graph)
    print(y_coordinate_for_temperature_graph)
    print(y_coordinate_for_density_graph)
    print(y_coordinate_for_velocity_graph)
    print(mach_of_gas_in_sections)
    print(y_coordinate_for_lambda_graph)
    print(y_coordinate_for_multiplication_rho_and_velocity_graph)

    # list_of_pressure = creating_list_of_pressure(list_of_data)

    # print(list_of_pressure, list_of_temperature, list_of_ro, list_of_velocity, list_of_lymbda, list_of_m, list_of_row)
