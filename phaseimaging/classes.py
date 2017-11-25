import numpy as np
from .utils import import_specimen
from .utils import lambda_to_accel_volt
import copy
from .imaging import *
from .plot import *


class Image:
    """ Generic parent class for images.

    Attributes:
        resolution (list of int): Two/three element list giving the dimensions of the image in pixels/voxels.
        width (list of float): Two/three element list giving the dimensions of the image in length units.
        image (ndarray): The image data.
        """
    def __init__(self, resolution, width, dtype=float, name=None):
        self.resolution = resolution
        self.width = width
        self.image = np.zeros(resolution, dtype=dtype)
        self.name = name

    def plot(self, limits=None):
        plot_image(self.image, limits, title=self.name)

    def apodise(self, rad_sup=0.5):
        self.image = apodise(self.image, rad_sup)

    def flip(self, axis=0):
        self.image = np.flip(self.image, axis=axis)


class Intensity(Image):
    """ Intensity measurement. """
    def __init__(self, resolution, width, defocus, incident=1, name=None):
        assert len(resolution) == len(width) == 2
        super().__init__(resolution, width, name=name)
        self.defocus = defocus
        self.incident = incident
        self.shape = self.resolution

    def add_noise(self, sigma):
        self.image = add_noise(self.image, self.incident, sigma)

    def transfer(self, phase, beam):
        self.image = transfer_image(self.defocus, beam.wavelength, self.width, phase.image)


class Phase(Image):
    def __init__(self, resolution, width, defocus=0, name=None):
        assert len(resolution) == len(width) == 2
        super().__init__(resolution, width, dtype=complex, name=name)
        self.defocus = defocus
        self.k_kernel = None
        self.inverse_k_squared_kernel = None


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

    def remove_offset(self, rad_inf=0.45, rad_sup=0.5):
        self.image = remove_offset(self.image, rad_inf=rad_inf, rad_sup=rad_sup)

    def normalised_rms_error(self, exact, display=False):
        return normalised_rms_error(exact.image, self.image, display=display, name=self.name)

class Wavefield(Image):
    def __init__(self, resolution, width, defocus=0, name=None):
        assert len(resolution) == len(width) == 2
        super().__init__(resolution, width, dtype=complex, name=name)
        self.defocus = defocus


class Specimen(Image):
    def __init__(self, width, mean_inner_potential=0, magnetisation=0, mhat=(0, 0, 1), specimen_file=None, moment_file=None, name=None):

        if specimen_file is not None:
            self.image = import_specimen(specimen_file)
        if moment_file is not None:
            self.moment = import_moment(moment_file)
        self.resolution = self.image.shape
        self.width = width
        self.mean_inner_potential = mean_inner_potential
        self.magnetisation = magnetisation
        self.mhat = mhat
        self.name=name

    def plot(self, limits=None):
        warnings.warn("Specimens can't be plotted in full; plotting slice only.")
        plot_image(self.image[:, :, int(self.resolution[2]/2 - 1)], limits, title=self.name)

    def rotate(self, angle, axis=0):
        axes = tuple(np.array((0, 1, 2))[np.where(np.array((0, 1, 2)) != axis)])
        self.image = rotate(input=self.image, angle=angle, reshape=False, cval=0, axes=axes)
        self.mhat = rotate_vector(self.mhat, angle=np.deg2rad(angle), axis=axis)

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

    def apodise(self, rad_sup=0.5):
        for intensity in self.intensities:
            intensity.apodise(rad_sup=rad_sup)


