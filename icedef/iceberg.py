"""Creates iceberg object."""

from typing import NamedTuple
import numpy as np


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

class Position:

    def __init__(self, longitude):
        self._longitude = longitude
        self._x = longitude * 2

    @property
    def longitude(self):
        return self._longitude

    @longitude.setter
    def longitude(self, value):
        self._x = value * 2
        self._longitude = value

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, value):
        self._longitude = value / 2
        self._x = value


class Velocity(NamedTuple):
    """Creates NamedTuple object with iceberg velocity information."""

    x: 'float'
    y: 'float'


class Iceberg:
    """Creates iceberg object."""

    WATERLINE_LENGTH_RANGE_BY_SIZE = {'LG': (120, 200)}
    SAIL_HEIGHT_RANGE_BY_SIZE = {'LG': (45, 75)}
    HEIGHT_TO_DRAFT_RATIO_BY_SHAPE = {'TAB': 0.2}
    SHAPE_FACTOR_BY_SHAPE = {'TAB': 0.2}

    def __init__(self, time, position, velocity, **kwargs):

        self._time = time
        self._position = position
        self.velocity = velocity

        self.name = kwargs.get('name', None)
        self.size = kwargs.get('size', 'LG')
        self.shape = kwargs.get('shape', 'TAB')

        self.history = {'time': []}

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, lonlat):
        self._position.longitude = lonlat[0]
        self._position.latitude = lonlat[1]

    @property
    def time(self):
        return self._time

    @time.setter
    def time(self, value):
        self.history['time'].append(self._time)
        self._time = value

    @property
    def waterline_length(self):
        """Return the mean waterline length for the size declared."""
        return np.mean(self.WATERLINE_LENGTH_RANGE_BY_SIZE[self.size])

    @property
    def sail_height(self):
        """Return the mean sail height for the size declared."""
        return np.mean(self.SAIL_HEIGHT_RANGE_BY_SIZE[self.size])

    @property
    def sail_area(self):
        """Return the area of rectangular sail."""
        return self.waterline_length * self.sail_height

    @property
    def height_to_draft_ratio(self):
        """Return the height to draft ratio for the shape declared."""
        return self.HEIGHT_TO_DRAFT_RATIO_BY_SHAPE[self.shape]

    @property
    def shape_factor(self):
        """Return the shape factor for the shape declared."""
        return self.SHAPE_FACTOR_BY_SHAPE[self.shape]

    @property
    def keel_depth(self):
        """Return the keel depth for the shape and sail height declared."""
        h2d_ratio = self.HEIGHT_TO_DRAFT_RATIO_BY_SHAPE[self.shape]
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
