import numpy as np

class Image:

    def __init__(self, image, resolution, width, dtype=float):
        self.image = np.zeros(resolution, dtype=dtype)
        self.resolution = resolution
        self.width = width

class Intensity(Image):

    def __init__(self, resolution, width, defocus):
        super.__init__(self, resolution, width)
        self.defocus = defocus

class Phase(Image):

    super.__init__(self, resolution, width)

class Wavefield(Image):

    super.__init__(self, resolution, width, dtype=complex)