"""Creates iceberg object."""

from typing import NamedTuple
import numpy as np

WATERLINE_LENGTH_RANGE_BY_SIZE = {'LG': (120, 200)}
SAIL_HEIGHT_RANGE_BY_SIZE = {'LG': (45, 75)}
HEIGHT_TO_DRAFT_RATIO_BY_SHAPE = {'TAB': 0.2}
SHAPE_FACTOR_BY_SHAPE = {'TAB': 0.2}


class Position:

    def __init__(self, latlon):
        self._latlon = latlon
        self._xy = (latlon[1]*2, latlon[0]*2)

    @property
    def latlon(self):
        return self._latlon

    @latlon.setter
    def latlon(self, value):
        self._xy = (value[1]*2, value[0]*2)
        self._latlon = value

    @property
    def xy(self):
        return self._xy

    @xy.setter
    def xy(self, value):
        self._latlon = (value[1]/2, value[0]/2)
        self._xy = value


class IcebergGeometry:

    def __init__(self, size, shape):

        self.size = size
        self.shape = shape

    @property
    def waterline_length(self):
        """Return the mean waterline length for the size declared."""
        return np.mean(WATERLINE_LENGTH_RANGE_BY_SIZE[self.size])

    @property
    def sail_height(self):
        """Return the mean sail height for the size declared."""
        return np.mean(SAIL_HEIGHT_RANGE_BY_SIZE[self.size])

    @property
    def sail_area(self):
        """Return the area of rectangular sail."""
        return self.waterline_length * self.sail_height

    @property
    def height_to_draft_ratio(self):
        """Return the height to draft ratio for the shape declared."""
        return HEIGHT_TO_DRAFT_RATIO_BY_SHAPE[self.shape]

    @property
    def shape_factor(self):
        """Return the shape factor for the shape declared."""
        return SHAPE_FACTOR_BY_SHAPE[self.shape]

    @property
    def keel_depth(self):
        """Return the keel depth for the shape and sail height declared."""
        h2d_ratio = HEIGHT_TO_DRAFT_RATIO_BY_SHAPE[self.shape]
        sail_height = self.sail_height
        return sail_height / h2d_ratio

    @property
    def mass(self):
        """Return the mass using formula from Rudkin, 2005."""
        factor = self.shape_factor
        length = self.waterline_length
        height = self.sail_height
        return 7.12e3 * factor * length ** 2 * height

    @property
    def keel_area(self):
        """Return the rectangular keel area."""
        return self.waterline_length * self.keel_depth


class Iceberg:
    """Creates iceberg object."""

    def __init__(self, time, position, velocity, geometry, **kwargs):

        self._time = time
        self.position = position
        self.velocity = velocity
        self.geometry = geometry
        self.name = kwargs.get('name', None)
        self.history = {'time': []}

    @property
    def time(self):
        return self._time

    @time.setter
    def time(self, value):
        self.history['time'].append(self._time)
        self._time = value


def quickstart(time, latlon, **kwargs):

    velocity = kwargs.get('velocity', (0, 0))
    size = kwargs.get('size', 'LG')
    shape = kwargs.get('shape', 'TAB')

    geometry = IcebergGeometry(size, shape)
    position = Position(latlon)
    iceberg = Iceberg(time, position, velocity, geometry)

    return iceberg


# class Velocity(NamedTuple):
#     """Creates NamedTuple object with iceberg velocity information."""
#
#     x: 'float'
#     y: 'float'

# class Position(NamedTuple):
#     """Creates NamedTuple object with iceberg position information."""
#
#     # init_longitude: 'float'
#     # _longitude: 'float' = init_longitude
#     longitude: 'float'
#
#     @property
#     def longitude(self):
#         return self.x / 2
#
#     @longitude.setter
#     def longitude(self, value):
#         self.longitude = value
#
#     @property
#     def x(self):
#         return self.longitude * 2
