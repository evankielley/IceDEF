"""Creates iceberg object."""

import numpy as np

WATERLINE_LENGTH_RANGE_BY_SIZE = {'LG': (120, 200)}
SAIL_HEIGHT_RANGE_BY_SIZE = {'LG': (45, 75)}
HEIGHT_TO_DRAFT_RATIO_BY_SHAPE = {'TAB': 0.2}
SHAPE_FACTOR_BY_SHAPE = {'TAB': 0.5}

EARTH_RADIUS = 6371e3  # meters


class Position:

    def __init__(self, latitude, longitude):
        self._latitude = latitude
        self._longitude = longitude
        self._x = 0
        self._y = 0

    @property
    def latitude(self):
        return self._latitude

    @latitude.setter
    def latitude(self, value):
        # dlat = value - self._latitude
        # self._y += dlat_to_dy(dlat)
        self._latitude = value

    @property
    def longitude(self):
        return self._longitude

    @longitude.setter
    def longitude(self, value):
        # dlon = value - self._longitude
        # lat = self._latitude
        # self._x += dlon_to_dx(dlon, lat)
        self._longitude = value

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, value):
        # dx = value - self._x
        # lat = self._latitude
        # self._longitude += dx_to_dlon(dx, lat)
        self._x = value

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, value):
        # dy = value - self._y
        # self._latitude += dy_to_dlat(dy)
        self._y = value


class Velocity:

    def __init__(self, vx, vy):
        self.x = vx
        self.y = vy


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
        self.time = time
        self.latitude, self.longitude = position
        self.eastward_velocity, self.northward_velocity = velocity
        self.geometry = geometry
        self.name = kwargs.get('name', None)
        self.history = {'time': [], 'latitude': [], 'longitude': []}

    def update_history(self):
        self.history['time'].append(self.time)
        self.history['latitude'].append(self.latitude)
        self.history['longitude'].append(self.longitude)


def quickstart(time, position, **kwargs):
    # position must be tuple of (latitude, longitude)
    velocity = kwargs.get('velocity', (0, 0))
    size = kwargs.get('size', 'LG')
    shape = kwargs.get('shape', 'TAB')
    geometry = IcebergGeometry(size, shape)
    iceberg = Iceberg(time, position, velocity, geometry)
    return iceberg


def old_quickstart(time, latitude, longitude, **kwargs):
    velocity = kwargs.get('velocity', (0, 0))
    velocity = Velocity(velocity[0], velocity[1])
    size = kwargs.get('size', 'LG')
    shape = kwargs.get('shape', 'TAB')

    geometry = IcebergGeometry(size, shape)
    position = Position(latitude, longitude)
    iceberg = Iceberg(time, position, velocity, geometry)

    return iceberg


def dx_to_dlon(dx, lat):
    return dx / (np.cos((lat * np.pi) / 180)) * (180 / (np.pi * EARTH_RADIUS))


def dy_to_dlat(dy):
    return dy * (180 / (np.pi * EARTH_RADIUS))