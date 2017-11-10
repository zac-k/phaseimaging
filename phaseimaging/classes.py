import numpy as np
from .utils import import_specimen
from .utils import lambda_to_accel_volt
import copy
from .imaging import *
from .plot import *

class Image:
    def __init__(self, resolution, width, dtype=float):
        self.resolution = resolution
        self.width = width

    def plot(self, limits=None):
        plot_image(self.image, limits)


class Intensity(Image):
    def __init__(self, resolution, width, defocus, incident=1):
        assert len(resolution) == len(width) == 2
        super().__init__(resolution, width)
        self.defocus = defocus
        self.incident = incident
        self.image = None
        self.shape = self.resolution

    def add_noise(self, sigma):
        self.image = add_noise(self.image, self.incident, sigma)

    def transfer(self, phase, beam):
        self.image = transfer_image(self.defocus, beam.wavelength, self.width, phase.image)


class Phase(Image):
    def __init__(self, resolution, width, defocus=0):
        assert len(resolution) == len(width) == 2
        super().__init__(resolution, width)
        self.defocus = defocus
        self.k_kernel = None
        self.inverse_k_squared_kernel = None
        self.image = np.zeros(resolution, dtype=complex)

    def construct_kernels(self, reg_param=None):
        self.k_kernel = construct_k_kernel(self.resolution,self.width)
        self.inverse_k_squared_kernel = construct_inverse_k_squared_kernel(self.resolution, self.width, reg_param)

    def project_electrostatic(self, specimen, beam):
        self.image += project_electrostatic_phase(specimen.image, beam.accel_volt, specimen.mean_inner_potential, specimen.width)

    def project_magnetic(self, specimen, beam):
        self.image += project_magnetic_phase(specimen.image, specimen.mhat, specimen.magnetisation, specimen.width, k_kernel=self.k_kernel, inverse_k_squared_kernel=self.inverse_k_squared_kernel)

    def retrieve_phase_tie(self, tfs, beam):
        if tfs.len % 2 == 1:
            image_in = tfs.intensities[int((tfs.len-1) / 2)].image
        else:
            image_in = None

        # Determine defocus of retrieved phase
        self.defocus = (tfs.defoci[0] + tfs.defoci[-1]) / 2
        self.image = retrieve_phase_tie(beam.wavelength,
                                        self.width,
                                        tfs.derivative,
                                        image_in)

class Wavefield(Image):
    def __init__(self, resolution, width, defocus=0):
        assert len(resolution) == len(width) == 2
        super().__init__(resolution, width, dtype=complex)
        self.defocus = defocus


class Specimen(Image):
    def __init__(self, width, mean_inner_potential=0, magnetisation=0, mhat=(0, 0, 1), specimen_file=None):

        if specimen_file is not None:
            self.image = import_specimen(specimen_file)
        self.resolution = self.image.shape
        super().__init__(self.resolution, width)
        self.mean_inner_potential = mean_inner_potential
        self.magnetisation = magnetisation
        self.mhat = mhat

    def plot(self, limits=None):
        warnings.warn("Specimens can't be plotted in full; plotting slice only.")
        plot_image(self.image[:, :, int(self.resolution[2]/2 - 1)], limits)

class Beam():
    def __init__(self, wavelength):
        self.wavelength = wavelength
        self.accel_volt = lambda_to_accel_volt(wavelength)


class ThroughFocalSeries():
    def __init__(self, resolution, width, defoci, incident=1):
        self.intensities = []
        self.derivative = None
        self.defoci = defoci
        self.len = len(defoci)
        for defocus in defoci:
            self.intensities.append(Intensity(resolution, width, defocus, incident))
        spacing = copy.copy(defoci)
        if self.len % 2 == 0:
            spacing.insert((self.len - 1) / 2, (defoci[int(self.len / 2)] + defoci[int(self.len / 2 - 1)]) / 2)
        for i in range(1, len(spacing)-1):
            assert spacing[i] - spacing[i-1] == spacing[i+1] - spacing[i]


    def transfer_images(self, phase, beam):
        for intensity in self.intensities:
            intensity.transfer(phase, beam)

    def add_noise(self, sigma):
        for intensity in self.intensities:
            intensity.add_noise(sigma)

    def compute_derivative(self):
        # todo: this needs to be generalised
        self. derivative = intensity_derivative(self.intensities[0].image,
                                                self.intensities[-1].image,
                                        (self.intensities[-1].defocus - self.intensities[0].defocus) / 2)


