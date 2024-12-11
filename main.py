"Launch module. Connects all modules together"

import res.calculations as calc
import res.constants as consts
from res import draw


def main() -> None:
    """main-function"""

    user_choice: str = input()
    if user_choice == "0":
        pass

    latitudes: list[float]
    longitudes: list[float]

    if user_choice == "1":
        latitudes, longitudes = calc.CoordinatesCalculator.coordinates_for_sphere()
        plot_builder: draw.Plot3D = draw.SphericalPlot(10, 8)
        plot_builder.show_sphere(consts.EARTH_RADIUS)
        plot_builder.show_trajectory(longitudes, latitudes, consts.EARTH_RADIUS, "red")
    elif user_choice == "2":
        latitudes, longitudes = calc.CoordinatesCalculator.coordinates_for_map()
        plot_builder: draw.Plot3D = draw.MapPlot(10, 8)
        plot_builder.show_map()
        plot_builder.show_trajectory(longitudes, latitudes, "red")
    elif user_choice == "3":
        plot_builder: draw.Plot3D = draw.EllipsoidsPlot(10, 8)
        plot_builder.show_ellipsoid(consts.a, consts.b)
        plot_builder.show_trajectory(
            "red", **dict(zip(("x", "y", "z"), calc.CoordinatesCalculator.coordinates_for_ellipsoid()))
        )
    else:
        print("Incorrect input")


if __name__ == "__main__":
    main()
