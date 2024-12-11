"""Module for base astronomic calculations"""

import numpy as np
from res import constants


class CoordinatesCalculator:
    """Class for calculating result coordinates"""

    @staticmethod
    def coordinates_for_map() -> tuple[list[float]]:
        """Calculates coordinates for map model

        Returns:
            tuple[list[float]]: calculated coordinates
        """
        period_of_rotation: float = 2 * np.pi * np.sqrt(constants.SEMIMAJOR_AXIS**3 / constants.MU)

        latitudes: list[float] = []
        longitudes: list[float] = []

        time: float
        for time in np.linspace(0, period_of_rotation, constants.AMOUNT_OF_POINTS):

            true_anomaly: float = MathsEquations.find_true_anomaly(time)
            radius_vector: float = MathsEquations.find_radius_vector(time)

            geocentric_coordinates: constants.NDArrayFloat = SystemConverter.perifocal_to_geocentric(
                np.array([radius_vector * np.cos(true_anomaly), radius_vector * np.sin(true_anomaly), 0]),
                constants.INCLINATION,
                constants.RAAN,
                constants.PERIGEE_ARGUMENT,
            )

            latitude, longitude = SystemConverter.geocentric_to_geographic(
                geocentric_coordinates * (constants.EARTH_RADIUS / np.linalg.norm(geocentric_coordinates)), time
            )
            latitudes.append(latitude)
            longitudes.append(longitude)

        return latitudes, longitudes

    @staticmethod
    def coordinates_for_ellipsoid() -> tuple[list[float]]:
        """Calculates coordinates for ellipsoid model

        Returns:
            tuple[list[float]]: calculated coordinates
        """
        period_of_rotation: float = MathsEquations.find_rotation_period()
        latitudes: list[float] = []
        longitudes: list[float] = []

        time: float
        for time in np.linspace(0, period_of_rotation, constants.AMOUNT_OF_POINTS):

            true_anomaly: float = MathsEquations.find_true_anomaly(time)
            radius_vector: float = MathsEquations.find_radius_vector(time)

            latitude, longitude = MathsEquations.project_trajectory_vectorized(
                SystemConverter.perifocal_to_geocentric(
                    np.array([radius_vector * np.cos(true_anomaly), radius_vector * np.sin(true_anomaly), 0]),
                    constants.INCLINATION,
                    constants.RAAN,
                    constants.PERIGEE_ARGUMENT,
                ),
                time,
            )
            latitudes.append(latitude)
            longitudes.append(longitude)

        return SystemConverter.geodetic_to_ecef(latitudes, longitudes, constants.a, constants.b)

    @staticmethod
    def coordinates_for_sphere() -> tuple[list[float]]:
        """Calculates coordinates for sphere model

        Returns:
            tuple[list[float]]: calculated coordinates
        """
        period_of_rotation: float = MathsEquations.find_rotation_period()
        latitudes: list[float] = []
        longitudes: list[float] = []

        time: float
        for time in np.linspace(0, period_of_rotation, constants.AMOUNT_OF_POINTS):

            true_anomaly: float = MathsEquations.find_true_anomaly(time)
            radius_vector: float = MathsEquations.find_radius_vector(time)
            geocentric_coordinates: constants.NDArrayFloat = SystemConverter.perifocal_to_geocentric(
                np.array([radius_vector * np.cos(true_anomaly), radius_vector * np.sin(true_anomaly), 0]),
                constants.INCLINATION,
                constants.RAAN,
                constants.PERIGEE_ARGUMENT,
            )
            latitude: float
            longitude: float
            latitude, longitude = SystemConverter.geocentric_to_geographic(
                geocentric_coordinates * (constants.EARTH_RADIUS / np.linalg.norm(geocentric_coordinates)), time
            )
            latitudes.append(latitude)
            longitudes.append(longitude)

        return latitudes, longitudes
class SystemConverter:
    """Class for converting from one system of coordinates to another"""

    @staticmethod
    def perifocal_to_geocentric(
        perifocal_coordinates: constants.NDArrayFloat, inclination: float, raan: float, perigee_argument: float
    ) -> constants.NDArrayFloat:
        """Refactors perifocal to geocentric coordinates

        Args:
            perifocal_coordinates (float): perifocal coordinates
            inclination (float): inclination of satellite's orbit (in degrees)
            raan (float): longitude of ascending node (in degrees)
            perigee_argument (float): argument of perigee (in degrees)

        Returns:
            NDArrayFloat: geocentric coordinates
        """
        cos_o: float = np.cos(-raan)
        sin_o: float = np.sin(-raan)
        cos_i: float = np.cos(-inclination)
        sin_i: float = np.sin(-inclination)
        cos_w: float = np.cos(-perigee_argument)
        sin_w: float = np.sin(-perigee_argument)

        rotation_matrix: constants.NDArrayFloat = np.array(
            [
                [cos_o * cos_w - sin_o * sin_w * cos_i, -cos_o * sin_w - sin_o * cos_w * cos_i, sin_o * sin_i],
                [sin_o * cos_w + cos_o * sin_w * cos_i, -sin_o * sin_w + cos_o * cos_w * cos_i, -cos_o * sin_i],
                [sin_w * sin_i, cos_w * sin_i, cos_i],
            ]
        )
        return rotation_matrix @ perifocal_coordinates

    @staticmethod
    def geocentric_to_geographic(geocentric_coordinates: float, time: float) -> constants.NDArrayFloat:
        """Calculates geographic coordinates from geocentric coordinates

        Args:
            geocentric_coordinates (float): geocentric coordinates (in degrees)
            time (float): time-argument

        Returns:
            NDArrayFloat: geographical coordinates
        """
        x: float
        y: float
        z: float
        x, y, z = geocentric_coordinates
        return (
            np.degrees(np.arctan(z / np.sqrt(x**2 + y**2))),
            (np.degrees(np.arctan2(y, x)) - np.degrees(constants.OMEGA_EARTH * time)) % 360,
        )

    @staticmethod
    def geodetic_to_ecef(
        latitude: constants.NDArrayFloat,
        longitude: constants.NDArrayFloat,
        semimajor_axis: float,
        semiminor_axis: float,
    ) -> tuple[float]:
        """Converts geodetic coordinates to ecef coordinates (Earth - Centered, Earth - Fixed)

        Args:
            latitude (NDArrayFloat): latitudes of geodetic coordinates (in degrees)
            longitude (NDArrayFloat): longitudes of geodetic coordinates (in degrees)
            semimajor_axis (float): semimajor axis of an ellipsoid Earth surface
            semiminor_axis (float): semiminor axis of an ellipsoid Earth surface

        Returns:
            tuple[float]: coordinates in ecef
        """
        latitude: constants.NDArrayFloat = np.radians(latitude)
        longitude: constants.NDArrayFloat = np.radians(longitude)

        radius_in_current_latitude: constants.NDArrayFloat = constants.a / np.sqrt(
            1 - (1 - (semiminor_axis**2 / semimajor_axis**2)) * np.sin(latitude) ** 2
        )

        return (
            radius_in_current_latitude * np.cos(latitude) * np.cos(longitude),
            radius_in_current_latitude * np.cos(latitude) * np.sin(longitude),
            (semiminor_axis**2 / semimajor_axis**2 * radius_in_current_latitude) * np.sin(latitude),
        )


class MathsEquations:
    """Class for solving mathematical equations"""

    @staticmethod
    def solve_kepler_equation(anomaly: float) -> float:
        """Solves Kepler equation

        Args:
            anomaly (float): argument of Kepler's equation

        Returns:
            float: equation's solution
        """
        root: float
        last_root: float
        root = last_root = anomaly

        while abs(root - last_root) >= constants.PRECISION_FOR_KEPLER_EQUATION:
            last_root = root
            root = constants.ECCENTRICITY * np.sin(last_root) + anomaly

        return root

    @staticmethod
    def project_to_ellipsoid_vectorized(time: float, **kwargs: float) -> tuple[float]:
        """Projects point to ellipsoid surface

        Args:
            time (float): time-argument

        Returns:
            tuple[float]: projected point
        """
        scale: float = 1.0 / np.sqrt(
            ((kwargs["x"] ** 2 + kwargs["y"] ** 2) / constants.a**2 + (kwargs["z"] ** 2) / constants.b**2)
        )
        ellipsoid_x: float = scale * kwargs["x"]
        ellipsoid_y: float = scale * kwargs["y"]
        ellipsoid_z: float = scale * kwargs["z"]

        latitude: float = np.arctan2(ellipsoid_z, np.sqrt(ellipsoid_x**2 + ellipsoid_y**2) * (1 - constants.e))
        for _ in range(5):
            latitude = np.arctan2(
                ellipsoid_z + constants.e * constants.b * np.sin(latitude), np.sqrt(ellipsoid_x**2 + ellipsoid_y**2)
            )

        latitude: float = np.degrees(latitude)
        return (
            latitude,
            (np.degrees(np.arctan2(ellipsoid_y, ellipsoid_x)) - np.degrees(constants.OMEGA_EARTH * time)) % 360,
        )

    @staticmethod
    def project_trajectory_vectorized(trajectory: constants.NDArrayFloat, time: float) -> constants.NDArrayFloat:
        """Projects trajectory to ellipsoid surface

        Args:
            trajectory (NDArrayFloat): trajectory of satellite
            time (float): time-argument

        Returns:
            NDArrayFloat: projected trajectory
        """
        return MathsEquations.project_to_ellipsoid_vectorized(time, **dict(zip(("x", "y", "z"), trajectory)))

    @staticmethod
    def find_rotation_period() -> float:
        """Calculates period of rotation

        Returns:
            float: period of rotation
        """
        return 2 * np.pi * np.sqrt(constants.SEMIMAJOR_AXIS**3 / constants.MU)

    @staticmethod
    def find_true_anomaly(time: float) -> float:
        """Calculates true anomaly of satellite in current time

        Args:
            time (float): time-argument

        Returns:
            float: true anomaly
        """
        period_of_rotation: float = MathsEquations.find_rotation_period()
        average_anomaly = 2 * np.pi * time / period_of_rotation
        true_anomaly = 2 * np.arctan(
            np.sqrt((1 + constants.ECCENTRICITY) / (1 - constants.ECCENTRICITY))
            * np.tan(MathsEquations.solve_kepler_equation(average_anomaly) / 2)
        )
        if true_anomaly < 0 or average_anomaly > np.pi:
            true_anomaly += 2 * np.pi

        return true_anomaly

    @staticmethod
    def find_radius_vector(time: float) -> float:
        """Calculates radius of satellite in current time

        Args:
            time (float): time-argument

        Returns:
            float: radius-vector of satellite
        """
        return (
            constants.SEMIMAJOR_AXIS
            * (1 - constants.ECCENTRICITY**2)
            / (1 + constants.ECCENTRICITY * np.cos(MathsEquations.find_true_anomaly(time)))
        )
