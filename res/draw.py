"""Module for drawing calculated plots"""

import abc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from res import constants


class Plot(abc.ABC):
    """Base class for all plot builders"""

    def __init__(self, width: int, height: int) -> None:
        self._figure: plt.plot = plt.figure(figsize=(width, height))
        self._axes: plt.axis = self._figure.add_subplot(projection="3d")


class SphericalPlot(Plot):
    """Class for drawing sphere plot and satellites's trajectory on it"""

    def __init__(self, width: int, height: int) -> None:
        super().__init__(width, height)

        axis: str
        for axis in "xyz":
            getattr(self._axes, f"set_{axis}label")(axis)

        self._axes.set_title("Проекция траектории на Землю(шар)")
        self._axes.grid(True)
        self._axes.set_box_aspect([1, 1, 1])

    def show_sphere(self, radius: float) -> None:
        """Create a plot with sphere (Earth surface)

        Args:
            radius (float, optional): radius of Earth sphere.
        """

        from_zero_to_2pi: constants.NDArrayFloat = np.linspace(0, 2 * np.pi, 100)
        from_zero_to_pi: constants.NDArrayFloat = np.linspace(0, np.pi, 100)
        self._axes.plot_surface(
            radius * np.outer(np.cos(from_zero_to_2pi), np.sin(from_zero_to_pi)),
            radius * np.outer(np.sin(from_zero_to_2pi), np.sin(from_zero_to_pi)),
            radius * np.outer(np.ones(np.size(from_zero_to_2pi)), np.cos(from_zero_to_pi)),
            color="cyan",
            alpha=0.3,
            edgecolor="none",
        )

        max_lim: float = radius * 1.1
        axis: str
        for axis in "xyz":
            getattr(self._axes, f"set_{axis}lim")((-max_lim, max_lim))

    def show_trajectory(self, longitudes: list[float], latitudes: list[float], radius: float, point_color: str) -> None:
        """Shows the trajectory on the Earth's surface (sphere)

        Args:
            longitudes (list[float]): longitudes of satellite's orbit
            latitudes (list[float]): latitudes of satellite's orbit
            radius (float): radius of Earth sphere
            point_color (str): color of satellite's trajectory
        """
        longitude_radian: constants.NDArrayFloat = np.radians(longitudes)
        latitude_radian: constants.NDArrayFloat = np.radians(latitudes)
        self._axes.scatter(
            radius * np.cos(latitude_radian) * np.cos(longitude_radian),
            radius * np.cos(latitude_radian) * np.sin(longitude_radian),
            radius * np.sin(latitude_radian),
            color=point_color,
            s=50,
        )

        plt.show()


class MapPlot(Plot):
    """Class for drawing World's map and satellites's trajectory on it"""

    def __init__(self, width: int, height: int) -> None:
        super().__init__(width, height)

        plt.title("Проекция орбиты спутника на Землю")
        self.__map_plot: Basemap = Basemap(projection="cyl", resolution="l")

    def show_map(self) -> None:
        """Crates a map plot"""
        self.__map_plot.drawcoastlines()
        self.__map_plot.drawcountries()
        self.__map_plot.drawparallels(np.arange(-90.0, 91.0, 30.0), labels=[True, False, False, False])
        self.__map_plot.drawmeridians(np.arange(-180.0, 181.0, 60.0), labels=[False, False, False, True])

    def show_trajectory(self, longitudes: list[float], latitudes: list[float], point_color: str) -> None:
        """Show the trajectory on the Earth's surface (flat)

        Args:
            longitudes (list[float]): longitudes of satellite's orbit
            latitudes (list[float]): latitudes of satellite's orbit
        """
        longitudes = list(map(lambda longitude: (longitude + 180) % 360 - 180, longitudes))

        self.__map_plot.plot(longitudes, latitudes, point_color, linewidth=1)

        plt.show()


class EllipsoidsPlot(Plot):
    """Class for drawing ellipsoid and satellite's orbit on it"""

    def __init__(self, width, height):
        super().__init__(width, height)

        axis: str
        for axis in "xyz":
            getattr(self._axes, f"set_{axis}label")(axis)
        self._axes.set_title("Проекция траектории на Землю(эллипс)")

    def show_ellipsoid(self, semimajor_axis: float, semiminor_axis: float) -> None:
        """Creates a plot containing an ellipsoid

        Args:
            semimajor_axis (float): semimajor axis of ellipsoid
            semiminor_axis (float): semiminor axis of ellipsoid
        """
        from_zero_to_2pi: constants.NDArrayFloat
        from_minus_half_pi_to_half_pi: constants.NDArrayFloat
        from_zero_to_2pi, from_minus_half_pi_to_half_pi = np.meshgrid(
            np.linspace(0, 2 * np.pi, 100), np.linspace(-np.pi / 2, np.pi / 2, 100)
        )

        self._axes.plot_surface(
            semimajor_axis * np.cos(from_minus_half_pi_to_half_pi) * np.cos(from_zero_to_2pi),
            semimajor_axis * np.cos(from_minus_half_pi_to_half_pi) * np.sin(from_zero_to_2pi),
            semiminor_axis * np.sin(from_minus_half_pi_to_half_pi),
            rstride=4,
            cstride=4,
            color="lightblue",
            alpha=0.5,
            edgecolor="none",
        )

        axis: str
        for axis in "xyz":
            getattr(self._axes, f"set_{axis}lim")((-semimajor_axis, semimajor_axis))

    def show_trajectory(self, point_color: str, **kwargs: constants.NDArrayFloat) -> None:
        """Show the trajectory on the Earth's surface (ellipsoid)

        Args:
            point_color (str): color of satellite's orbit
        """
        self._axes.scatter(kwargs["x"], kwargs["y"], kwargs["z"], color=point_color, s=50)
        plt.show()
