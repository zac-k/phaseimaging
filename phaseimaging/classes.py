import numpy as np
from .utils import import_specimen
from .utils import lambda_to_accel_volt
import copy
from .imaging import *
from .plot import *


class Image:
    """ Generic parent class for images.

    Attributes:
        resolution (sequence of ints): Two/three element list giving the dimensions of the image
                                       in pixels/voxels.
        width (sequence of floats): Two/three element list giving the dimensions of the image in length units.
        image (ndarray): The image data.
        """
    def __init__(self, resolution, width, dtype=float, name=None):
        """
        Args:
            resolution (sequence of ints): Dimensions of the image in pixels.
            width (sequence of floats): Dimensions of the image in units of length.
            dtype (type|optional): The data type of the image pixels/voxels. Default
                                   is float.
            name (string|optional): Name of the image. Will be used for titles on plots.
                                    Default is None.
        """
        self.resolution = resolution
        self.width = width
        self.image = np.zeros(resolution, dtype=dtype)
        self.name = name

    def plot(self, limits=None):
        """
        Plot the image on screen.

        Args:
            limits (sequence of floats|optional): Two-element sequence providing the minimum
                                              and maximum values of the image. Values outside
                                              this range will be clipped at the min/max. If
                                              limits are not provided, the color axis will be
                                              scaled according to the max/min values in the
                                              image array.

        Returns:
            None
        """
        plot_image(self.image, limits, title=self.name)

    def apodise(self, rad_sup=0.5):
        """
        Apodise the image using a circular (or elliptical) window function.

        Args:
            rad_sup (float|optional): The fractional radius (or semi axes) of the window function.
                                      Default is 0.5, which is the largest window function that
                                      remains circular (elliptical) within the image.

        Returns:
            None
        """
        self.image = apodise(self.image, rad_sup)

    def flip(self, axis=0):
        """
        Flip the image about an axis.

        Args:
            axis (int|optional): The axis about which to flip the image. Default is the first axis.

        Returns:
            None
        """
        self.image = np.flip(self.image, axis=axis)


class Intensity(Image):
    """
    An intensity measurement.
    """
    def __init__(self, resolution, width, defocus, incident=1, name=None):
        """
        Args:
            resolution (sequence of ints): Dimensions of the image in pixels.
            width (sequence of floats): Dimensions of the image in units of length.
            defocus (float): The defocus of the image.
            incident (float|optional): Intensity of the incident beam. Used for computing
                                       noise. Default is unity.
            name (string|optional): Name of the image. Will be used for titles on plots.
                                    Default is None.
        """
        assert len(resolution) == len(width) == 2
        super().__init__(resolution, width, name=name)
        self.defocus = defocus
        self.incident = incident
        self.shape = self.resolution

    def add_noise(self, sigma):
        """
        Add Poisson noise to the micrographs (see https://doi.org/10.1103/PhysRevA.90.023859).

        Args:
            sigma (float): The fractional noise level in the incident beam; e.g., for 5% noise,
                           set sigma=0.05.

        Returns:
            None
        """
        self.image = add_noise(self.image, self.incident, sigma)

    def transfer(self, phase, beam):
        """
        Generate the image from an exit phase and incident beam.

        Args:
            phase (Phase): The exit phase at the image plane.
            beam (Beam): The incident beam.

        Returns:
            None
        """
        self.image = transfer_image(self.defocus, beam.wavelength, self.width, phase.image)


class Phase(Image):
    """
    A phase measurement.
    """
    def __init__(self, resolution, width, defocus=0, name=None):
        """
        Args:
            resolution (sequence of ints): Dimensions of the phase map in pixels.
            width (sequence of floats): Dimensions of the phase map in units of length.
            defocus (float|optional): The defocus of the phase map. Default is zero, which
                                      is typical in most circumstances.
            name (string|optional): Name of the phase map. Will be used for titles on plots.
                                    Default is None.
        """
        assert len(resolution) == len(width) == 2
        super().__init__(resolution, width, dtype=complex, name=name)
        self.defocus = defocus
        self.k_kernel = None
        self.inverse_k_squared_kernel = None


    def construct_kernels(self, reg_param=None):
        """
        Generate k_kernel and inverse_k_squared_kernel.

        Args:
            reg_param (float|optional): The regularisation parameter to use for constructing the
                                        inverse |k|^2 kernel. If none is provided, a default
                                        value based on the dimensions of the image will be used.

        Returns:
            None
        """
        self.k_kernel = construct_k_kernel(self.resolution, self.width)
        self.inverse_k_squared_kernel = construct_inverse_k_squared_kernel(self.resolution,
                                                                           self.width,
                                                                           reg_param)

    def project_electrostatic(self, specimen, beam):
        """
        Project a beam through a specimen and add the phase shift to the phase map.

        Args:
            specimen (Specimen): The specimen to project.
            beam (Beam): The beam to project through the specimen.

        Returns:
            None
        """
        self.image += project_electrostatic_phase(specimen.image,
                                                  beam.accel_volt,
                                                  specimen.mean_inner_potential,
                                                  specimen.width)

    def project_magnetic(self, specimen):
        """
        Project a beam through a magnetisation configuration and add the phase shift to the
        phase map.

        Args:
            specimen (Specimen): The specimen to project.

        Returns:
            None
        """
        self.image += project_magnetic_phase(specimen.image,
                                             specimen.magnetisation,
                                             specimen.width,
                                             mhat=specimen.mhat,
                                             moment=specimen.moment,
                                             k_kernel=self.k_kernel,
                                             inverse_k_squared_kernel=self.inverse_k_squared_kernel)

    def retrieve_phase_tie(self, tfs, beam):
        """
        Retrieve the phase from a through-focal-series.
        Args:
            tfs (ThroughFocalSeries): The through-focal series of intensities.
            beam (Beam): The incident beam.

        Returns:
            None
        """

        # If there is an odd number of images, use the central one as the in-focus image.
        if tfs.len % 2 == 1:
            image_in = tfs.intensities[int((tfs.len-1) / 2)].image
        else:
            image_in = None

        # Compute the through-focal derivative if it does not already exist.
        if tfs.derivative is None:
            tfs.compute_derivative()

        # Determine defocus of retrieved phase (usually zero).
        self.defocus = (tfs.defoci[0] + tfs.defoci[-1]) / 2

        # Set the image data to the retrieved phase map.
        self.image = retrieve_phase_tie(beam.wavelength,
                                        self.width,
                                        tfs.derivative,
                                        image_in)

    def remove_offset(self, rad_inf=0.45, rad_sup=0.5):
        """
        Remove any offset in the phase by assuming it to be zero over an annulus outside the specimen.

        Args:
            rad_inf (float|optional): Fractional inner radius of the annulus over which to compute
                                      the mean offset. Default is 0.45.
            rad_sup (float|optional): Fractional outer radius of the annulus over which to compute
                                      the mean offset. Default is 0.5.

        Returns:
            None
        """
        self.image = remove_offset(self.image, rad_inf=rad_inf, rad_sup=rad_sup)

    def normalised_rms_error(self, exact, display=False):
        """
        Compute the normalised rms error in the phase, relative to a reference phase.

        Args:
            exact (ndarray): The exact (reference) phase.
            display (bool|optional): If True, the error is displayed in the console.

        Returns:
            normalised_rms_error (float): The calculated error as a fraction.
        """
        return normalised_rms_error(exact.image, self.image, display=display, name=self.name)

class Wavefield(Image):
    """
    Slice of the wavefield.
    """
    def __init__(self, resolution, width, defocus=0, name=None):
        """
        Args:
            resolution (sequence of ints): Dimensions of the slice in pixels.
            width (sequence of floats): Dimensions of the slice in units of length.
            defocus (float|optional): The defocus of the slice. Default is zero, which
                                      is typical in most circumstances.
            name (string|optional): Name of the slice. Will be used for titles on plots.
                                    Default is None.
        """
        assert len(resolution) == len(width) == 2
        super().__init__(resolution, width, dtype=complex, name=name)
        self.defocus = defocus


class Specimen(Image):
    """
    Specimen
    """
    def __init__(self, width, mean_inner_potential=0, magnetisation=0, mhat=(0, 0, 1),
                 specimen_file=None,
                 moment_file=None,
                 name=None):
        """
        Args:
            width (sequence of floats): Dimensions of the specimen in units of length.
            mean_inner_potential (complex float): The mean inner potential of the specimen. An
                                                  imaginary component can be used to simulate
                                                  attenuation.
            magnetisation (float): Volume magnetisation strength within the specimen.
            mhat (sequence of floats|optional): Unit vector in direction of magnetisation (for
                                                uniform magnetisations. Ignored if moment is not
                                                None. Default is (0, 0, 1). Will be ignored if a
                                                moment distribution is provided.
            specimen_file (str): Path of specimen file.:
            moment_file (str): Path of moment file.
            name (string): Name of the specimen. Used for titling plots, etc.
        """
        self.image = None
        self._moment = None
        if specimen_file is not None:
            self.image = import_specimen(specimen_file)
        if moment_file is not None:
            self.moment = import_moment(moment_file)

        self.resolution = self.image.shape
        self.width = width
        self.mean_inner_potential = mean_inner_potential
        self.magnetisation = magnetisation
        self.mhat = mhat
        self.name = name

    def plot(self, limits=None):
        """
        Plot a slice of the specimen.
        Args:
            limits (sequence of floats|optional): Two-element sequence providing the minimum
                                                  and maximum values of the plot. Values outside
                                                  this range will be clipped at the min/max. If
                                                  limits are not provided, the color axis will be
                                                  scaled according to the max/min values in the
                                                  image array.

        Returns:
            None
        """
        warnings.warn("Specimens can't be plotted in full; plotting slice only.")
        plot_image(self.image[:, :, int(self.resolution[2]/2 - 1)], limits, title=self.name)

    def rotate(self, angle, axis=0):
        """
        Rotate the specimen mask, as well as any magnetisation or moment vector array, about
        one of the Cartesian axes.

        Args:
            angle (float): The angle of rotation in degrees.
            axis (int|optional): The axis about which to rotate. Default is zero.

        Returns:
            None
        """

        # Convert axis integer to orthogonal plane in R3.
        axes = tuple(np.array((0, 1, 2))[np.where(np.array((0, 1, 2)) != axis)])

        # Rotate specimen mask if it exists
        if self.image is not None:
            self.image = rotate(input=self.image, angle=angle, reshape=False, cval=0, axes=axes)

        # Rotate uniform magnetisation direction if it exists
        if self.mhat is not None:
            self.mhat = rotate_vector(self.mhat, angle=np.deg2rad(angle), axis=axis)

        # Rotate moment distribution array if it exists
        if self.moment is not None:
            self.moment = rotate(input=self.moment, angle=angle, reshape=False, cval=0, axes=axes)
            self.moment = rotate_vector(self.moment, angle=np.deg2rad(angle), axis=axis)

    def _mask_moment(self):
        """
        Set the moment distribution to zero outside the specimen. Executed automatically when
        the Specimen object is initialised with a moment file, or when the moment array property
        is set.

        Returns:
            None
        """
        assert self.image is not None
        self._moment = self.moment * self.image[..., np.newaxis]

    @property
    def moment(self):
        return self._moment

    @moment.setter
    def moment(self, moment):
        self._moment = moment
        if self.image is not None:
            self._mask_moment()


class Beam:
    """
    An indident beam.
    """
    def __init__(self, wavelength):
        """
        Args:
            wavelength (float): The wavelength of the incident beam.
        """
        self.wavelength = wavelength
        self.accel_volt = lambda_to_accel_volt(wavelength)


class ThroughFocalSeries:
    """
    A series of intensity measurements and associated parameters and properties.
    """

    def __init__(self, resolution, width, defoci, incident=1):
        """
        Args:
            resolution (sequence of ints): Two element sequence giving the dimensions of the images
                                           in pixels/voxels.
            width (sequence of floats): Two element sequence giving the dimensions of the images in
                                        units of length.

            defoci (sequence of floats): Sequence of equally spaced defoci, typically centred on zero.
            incident (float|optional): Incident intensity of the beam. Default is unity.
        """

        self.intensities = []
        self.derivative = None
        self.defoci = defoci
        self.len = len(defoci)

        # Construct a list of intensities using the supplied defoci
        for defocus in defoci:
            self.intensities.append(Intensity(resolution, width, defocus, incident))
        spacing = copy.copy(defoci)

        # Ensure that the defoci spacing is even, with the exception of a larger
        # distance between the two central defoci if there is no middle (usually
        # in-focus) intensity image.
        if self.len % 2 == 0:
            spacing.insert((self.len - 1) / 2, (defoci[int(self.len / 2)]
                                                + defoci[int(self.len / 2 - 1)]) / 2)
        for i in range(1, len(spacing)-1):
            assert spacing[i] - spacing[i-1] == spacing[i+1] - spacing[i]


    def transfer_images(self, phase, beam):
        """
        Generate the image data in the intensities list, by transferring the phase.

        Args:
            phase (Phase): The image plane phase.
            beam (Beam): The incident beam.

        Returns:
            None
        """
        for intensity in self.intensities:
            intensity.transfer(phase, beam)

    def add_noise(self, sigma):
        """
        Add Poisson noise to the images in the series (see https://doi.org/10.1103/PhysRevA.90.023859).

        Args:
            sigma (float): The fractional noise level in the incident beam; e.g., for 5% noise,
                       set sigma=0.05.

        Returns:
            None
        """
        for intensity in self.intensities:
            intensity.add_noise(sigma)

    def compute_derivative(self):
        """
        Compute the longitudinal derivative of the intensity using a two image central difference
        approximation.

        Returns:
            None
        """
        # todo: this needs to be generalised
        self. derivative = intensity_derivative(self.intensities[0].image,
                                                self.intensities[-1].image,
                                        (self.intensities[-1].defocus - self.intensities[0].defocus) / 2)

    def apodise(self, rad_sup=0.5):
        """
        Apodise the images in the series using a circular (or elliptical) window function.

        Args:
            rad_sup (float|optional): The fractional radius (or semi axes) of the window function.
                                      Default is 0.5, which is the largest window function that
                                      remains circular (elliptical) within the image.

        Returns:
            None
        """
        for intensity in self.intensities:
            intensity.apodise(rad_sup=rad_sup)


