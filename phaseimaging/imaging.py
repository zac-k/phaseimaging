import numpy as np
from .utils import *
from numpy import pi as PI
import warnings
from scipy.ndimage import rotate, zoom, shift
import copy

h = 6.63e-34  # Planck's constant
m0 = 9.11e-31  # Electron rest mass
mu0 = 4 * PI * 10**-7  # Free space permeability
e = -1.6e-19  # Electron charge
c = 3e8  # Speed of light
hbar = h / (2 * PI)  # Reduced Planck constant
phi0 = h / (2 * e)  # Flux quantum


def construct_k_squared_kernel(resolution, image_width):
    """
    Construct a Fourier space kernel of |k|^2 values.

    Args:
        resolution (sequence of ints): The dimensions of the required kernel in pixels.
        image_width (sequence of floats): The dimensions of the real space image in
                                          units of length.

    Returns:
        kernel (ndarray): The Fourier space |k|^2 kernel.
    """

    # todo: vectorise this function
    kernel = np.zeros(resolution, dtype=complex)
    for i in range(resolution[0]):
        for j in range(resolution[1]):
            i0 = i - resolution[0] / 2
            j0 = j - resolution[1] / 2
            da = [1 / x for x in image_width]
            kernel[i, j] = (i0 * i0) * da[0] * da[0] + j0 * j0 * da[1] * da[1]

    return kernel

def _set_transfer_function(defocus, wavelength, resolution=None, image_width=None, k_squared_kernel=None):
    """
    Sets the imaging system's transfer function for a given defocus.

    Args:
        defocus (float): Defocus that the transfer function with produce.
        wavelength (float): Wavelength of the beam.
        resolution (sequence of ints|optional): The dimensions of the transfer function
                                                in pixels. Must be provided if no k_squared_kernel
                                                is given.
        image_width (sequence of floats|optional): The dimensions of the point-spread-function
                                                   (the transfer function in real space) in
                                                   units of length. Must be provided if no
                                                   k_squared_kernel is given.
        k_squared_kernel (ndarray|optional): The k_squared_kernel can be provided here to save
                                             some computation time at high resolution.

    Returns:
        transfer_function (ndarray): The computed transfer function.
    """

    # Build |k|^2 kernel if it is not supplied.
    if k_squared_kernel is None:
        assert resolution is not None and image_width is not None
        k_squared_kernel = construct_k_squared_kernel(resolution, image_width)

    # Calculate the transfer function.
    return np.exp(1j * (PI * defocus * wavelength * k_squared_kernel))


def transfer_image(defocus, wavelength, image_width, phase=None, wavefield=None, is_image=True):
    # todo: change this function so it only transfers wavefields, and put the phase-to-wave and
    # todo: wave to intensity in different functions. Will also need changes to classes.
    """
    Uses the defocus and exact_phase (at the image plane) to produce an out of focus
    image or wavefield.

    Args:
        defocus (float): Defocus to transfer the image to.
        wavelength (float): Wavelength of the incident beam.
        image_width (sequence of floats): Dimensions of the image in units of length.
        phase (ndarray): Exit phase of the beam.
        is_image (bool|optional): Whether to return the intensity image. Returns the
                                  complex wavefield if False. Default is true.

    Returns:
        intensity|wavefunction (ndarray|complex ndarray): The intensity or wavefield
                                                          ---dependent on the value
                                                          of is_image---at the plane
                                                          of defocus.
    """

    # Make sure exactly one of the input types is provided.
    # The != is effectively an xor here.
    assert (phase is not None) != (wavefield is not None)

    # Compute the complex wavefield at the image plane.
    if phase is not None:
        wavefunction = np.exp(1j*phase)
    elif wavefield is not None:
        wavefunction = wavefield

    # Set the transfer function.
    transfer_function = _set_transfer_function(defocus,
                                               wavelength,
                                               resolution=wavefunction.shape,
                                               image_width=image_width)

    # Transfer the wavefield to the plane of defocus.
    wavefunction = convolve(wavefunction, transfer_function)

    # Return the intensity or complex wavefield.
    if is_image:
        return np.absolute(wavefunction) * np.absolute(wavefunction)
    else:
        return wavefunction


def construct_inverse_k_squared_kernel(resolution, image_width, reg_param=None):
    """
    Generate a Fourier space kernel of inverse |k|^2 values using Tikhonov regularisation.

    Args:
        resolution (sequence of ints): The dimensions of the kernel in pixels.
        image_width (sequence of floats): The real-space dimensions in units of length.
        reg_param (float|optional): Regularisation parameter in inverse units of length.
                                    A default value scaled according to the dimensions
                                    of the kernel will be generated if none is provided.

    Returns:
        kernel (ndarray): The regularised inverse |k|^2 kernel.
    """
    # todo: vectorise this function.
    # Generate a default regularisation parameter if none is provided.
    if reg_param is None:
        reg_param = 0.1 / (np.mean(image_width) * np.mean(resolution))

    # Construct a |k|^2 kernel
    kernel = np.zeros(resolution, dtype=complex)
    da = tuple([1 / x for x in image_width])
    for i in range(resolution[0]):
        for j in range(resolution[1]):
            i0 = i - resolution[0] / 2
            j0 = j - resolution[1] / 2

            kernel[i, j] = (i0 * i0 * da[0] * da[0]) + (j0 * j0 * da[1] * da[1])

    # Use the |k|^2 kernel to generate the regularised inverse |k|^2 kernel
    kernel = kernel / (kernel * kernel + reg_param ** 4)
    return kernel


def construct_k_kernel(resolution, image_width):
    """
    Generate a kernel of k-vectors.

    Args:
        resolution (sequence of ints): The dimensions of the kernel in pixels.
        image_width (sequence of floats): The (real-space) dimensions in units
                                          of length.

    Returns:
        kernel (ndarray): The kernel of k-vectors
    """
    # todo: vectorise this function.

    kernel = np.zeros(list((resolution[0], resolution[1], 2)), dtype=complex)
    da = tuple([1 / x for x in image_width])
    for i in range(resolution[0]):
        for j in range(resolution[1]):
            i0 = i - resolution[0] / 2
            j0 = j - resolution[1] / 2
            kernel[i, j, 0] = i0 * da[0]
            kernel[i, j, 1] = j0 * da[1]
    return kernel

def retrieve_phase_tie(
                       wavelength,
                       image_width,
                       intensity_derivative,
                       image_in=None,
                       image_intensity=1,
                       k_kernel=None,
                       inverse_k_squared_kernel=None,
                       reg_param=0.1,
                       reg_param_tie=None
                       ):
    """
    Utilise the transport-of-intensity equation to compute the phase from
    the longitudinal derivative of the intensity.

    Args:
        wavelength (float): Wavelength of the imaging beam.
        image_width (sequence of floats): Width of the images in units of length.
        intensity_derivative (ndarray): Longitudinal derivative of the intensity.
        image_in (ndarray|optional): The intensity at the focal plane. This is used to
                                     improve phase reconstruction by accounting for
                                     absorption contrast in images. If no image is
                                     supplied here, a less sophisticated algorithm
                                     will be used, which may be inaccurate for highly
                                     attenuating specimens.
        image_intensity (float|optional): If no image_in is supplied, this average
                                     intensity will be used to compute the phase.
                                     For weakly attenuating specimens, this will
                                     be closely approximated by the incident
                                     intensity, which will typically be normalised
                                     to unity; this is the default value.
        k_kernel (ndarray): Fourier space kernel of k vectors. Will be generated if
                            none is supplied.
        inverse_k_squared_kernel (ndarray): Regularised Fourier space kernel of
                                            inverse |k|^2 values. Will be generated
                                            if none is supplied.
        reg_param (float): Regularisation parameter to compensate for potential
                           singularities when dividing by image_in. The default
                           value of 0.1 is fine for most cases, assuming that
                           the intensities are normalised to unity.
        reg_param_tie (float): Regularisation parameter used in constructing the
                               inverse |k|^2 kernel, when none is provided. If this
                               is left as the default of None, a (typically)
                               reasonable value will be generated based on the
                               resolution and dimensions of the images.

    Returns:
        phase_retrieved (ndarray): The reconstructed phase map.
    """

    # If an in-focus image is supplied, make sure it is the same size as the intensity derivative.
    if image_in is not None:
        assert image_in.shape == intensity_derivative.shape

    # Set the resolution of the system.
    resolution = intensity_derivative.shape
    assert len(intensity_derivative.shape) == 2
    assert len(image_width) == 2

    # Construct any kernel that wasn't supplied.
    if k_kernel is None:
        k_kernel = construct_k_kernel(resolution, image_width)
    if inverse_k_squared_kernel is None:
        inverse_k_squared_kernel = construct_inverse_k_squared_kernel(resolution, image_width, reg_param_tie)

    # Use the in-focus image if it was supplied.
    if image_in is not None:
        regularised_inverse_intensity = image_in / (image_in * image_in + reg_param * reg_param)
        prefactor = (1. / wavelength) / (2. * PI)
        derivative_vec = np.zeros((resolution[0], resolution[1], 2), dtype=complex)
        derivative_vec[:, :, 0] = convolve(intensity_derivative, k_kernel[:, :, 0] *
                                           inverse_k_squared_kernel
                                           ) * regularised_inverse_intensity
        derivative_vec[:, :, 1] = convolve(intensity_derivative, k_kernel[:, :, 1] *
                                           inverse_k_squared_kernel
                                           ) * regularised_inverse_intensity

        derivative_vec[:, :, 0] = fft.fftshift(fft.fft2(derivative_vec[:, :, 0]))
        derivative_vec[:, :, 1] = fft.fftshift(fft.fft2(derivative_vec[:, :, 1]))
        derivative = dot_fields(k_kernel, derivative_vec) * inverse_k_squared_kernel
        phase_retrieved = prefactor * fft.ifft2(fft.ifftshift(derivative))
    # Otherwise assume a uniform in-focus intensity of `image_intensity`.
    else:
        prefactor = (1. / wavelength) / (2 * PI * image_intensity)
        filtered = convolve(intensity_derivative, inverse_k_squared_kernel)
        phase_retrieved = prefactor * filtered

    # Any imaginary components here are just be numerical errors, so get rid of them.
    phase_retrieved = np.real(phase_retrieved)
    return phase_retrieved


def project_electrostatic_phase(specimen,
                                accel_volt,
                                mean_inner_potential,
                                image_width):
    """
    Computes the image plane phase shift of the specimen.

    Args:
        specimen (ndarray): A specimen mask. Typically binary, but non-binary values can
                            be included for, say, softening edges using a low pass filter.
        accel_volt (float): The electron accelerating voltage (in Volts).
        mean_inner_potential (complex float): The mean inner potential of the specimen. An
                                              imaginary component can be used to simulate
                                              attenuation.
        image_width (sequence of floats): The dimensions of the specimen in units of length.

    Returns:
        phase (ndarray): The exit phase of the electron wavefunction.
    """
    # todo: alter this for optional x-ray imaging, or write a second version for that purpose.

    # Get the resolution of the specimen and ensure that it is three-dimensional.
    resolution = specimen.shape
    assert len(resolution) == 3

    # Determine wavelength and differential element.
    wavelength = accel_volt_to_lambda(accel_volt)
    dz = image_width[2] / resolution[2]

    # Project potential in the z-direction.
    return np.sum(specimen, axis=2) * PI/(accel_volt * wavelength) * mean_inner_potential * dz


def project_magnetic_phase(specimen,
                           magnetisation,
                           image_width,
                           mhat=None,
                           moment=None,
                           k_kernel=None,
                           inverse_k_squared_kernel=None):
    """
    Computes the magnetic phase from a moment distribution file or uniform magnetisation direction.

    Args:
        specimen (ndarray): Specimen mask file.
        magnetisation (float): Volume magnetisation strength.
        image_width (sequence of floats): Dimensions of image in units of length.
        mhat (sequence of floats|optional): Unit vector in direction of magnetisation (for
                                            uniform magnetisations. Ignored if moment is not
                                            None. Default is None.
        moment (ndarray): Three-dimentional array of vectors representing the magnetisation direction
                          as a function of space. Mask must have already been applied. Specimen array
                          has no effect if moment is not None. Values less than or greater than unity
                          in this array will affect the strength of the field in those voxels. Default
                          is None.
        k_kernel (ndarray): Two dimensional array of k-vectors. Will be constructed if not supplied
                            (default is None).
        inverse_k_squared_kernel (ndarray): Regularised inverse k-squared kernel. Will be constructed
                                            if not supplied (default is None).

    Returns:
        phase (ndarray): Phase map representing the magnetic phase shift.
    """
    # Set the resolution from the specimen.
    resolution = specimen.shape
    assert len(image_width) == 3

    # Make sure there is an argument that describes the magnetisation.
    assert mhat is not None or moment is not None

    # Construct the kernels if they are not supplied.
    if k_kernel is None:
        k_kernel = construct_k_kernel(resolution, image_width)
    if inverse_k_squared_kernel is None:
        inverse_k_squared_kernel = construct_inverse_k_squared_kernel(resolution[0:2], image_width[0:2])

    # Use the moment array, if it is given, to compute the phase.
    if moment is not None:
        if mhat is not None:
            warnings.warn("moment array is being used for magnetisation---ignoring mhat vector")
        moment_ = fft.fftn(moment, axes=[0, 1, 2]) * (image_width[2] / resolution[2])
        moment_ = fft.fftshift(moment_, axes=[0, 1, 2])
        moment_ = moment_[:, :, int(resolution[2] / 2), :]
        mhatcrossk_z = np.cross(moment_, k_kernel)[:, :, 2]

        # Set shape function to unity. Moment array must be pre-masked using specimen.mask_moment()
        D0 = 1
    else:
        # Use the magnetisation direction (mhat) if no moment array is given.
        mhat = mhat / np.linalg.norm(mhat)
        mhatcrossk_z = np.cross(mhat, k_kernel)[:, :, 2]
        D0 = fft.fftshift(fft.fftn(specimen))[:, :, int(resolution[2] / 2)] * (image_width[2] / resolution[2])

    # Finish computing the magnetic phase.
    phase = (1j * PI * mu0 * magnetisation / phi0) * D0 * inverse_k_squared_kernel * mhatcrossk_z
    phase = fft.ifftshift(phase)
    return np.real(fft.ifft2(phase))


def add_noise(image, i_in, sigma):
    """
    Incorporate Poisson noise into an image (see https://doi.org/10.1103/PhysRevA.90.023859).

    Args:
        image (ndarray): The image to add noise to.
        i_in (float): The incident intensity. Typically unity.
        sigma (float): The fractional noise level in the incident beam; e.g., for 5% noise,
                       set sigma=0.05.

    Returns:
        output (ndarray): The noisy image.
    """
    if sigma > 0:
        output = np.where(image >= 0, np.random.poisson(image /
                                                        (sigma*sigma * i_in)) * (sigma*sigma * i_in), 0)
    else:
        output = image
    return output


def apodise(image, rad_sup=0.5):
    """
    Apodise image with circular (or elliptical) window function.

    Args:
        image (ndarray): The image to apodise.
        rad_sup (float|optional): The fractional radius (or semi axes) of the window function. Default is
                                  0.5, which is the largest window function that remains circular
                                  (elliptical) within the image.

    Returns:
        output (ndarray): The apodised image.
    """

    # Create meshgrid of pixel locations centred on the image.
    x = np.arange(int(-image.shape[0] / 2), int(image.shape[0] / 2))
    y = np.arange(int(-image.shape[1] / 2), int(image.shape[1] / 2))
    X,Y = np.meshgrid(x / image.shape[0], y / image.shape[1])

    # Set output equal to the image inside the window, zero elsewhere.
    output = np.where(X**2 + Y**2 < rad_sup**2, image, 0)
    return output


def remove_offset(phase, rad_inf=0.45, rad_sup=0.5):
    """
    Adjust a phase map by a constant value. Assumes that the phase shift is zero over some annulus
    around the specimen.

    Args:
        phase (ndarray): The phase to be adjusted.
        rad_inf (float|optional): Inner radius of the annulus as a fraction of the width of the
                                  phase map. Default is 0.45.
        rad_sup (float|optional): Outer radius of the annulus as a fraction of the width of the
                                  phase map. Default is 0.5, which is the largest outer radius
                                  that can remain circular (or elliptical) within the phase map.

    Returns:
        output_phase (ndarray): The adjusted phase map.
    """

    # Create meshgrid of pixel locations centred on the image.
    x = np.arange(int(-phase.shape[0] / 2), int(phase.shape[0] / 2))
    y = np.arange(int(-phase.shape[1] / 2), int(phase.shape[1] / 2))
    X, Y = np.meshgrid(x / phase.shape[0], y / phase.shape[1])

    # Find the average value of the phase within the annulus.
    offset = phase[np.where(np.logical_and(rad_sup ** 2 > X ** 2 + Y ** 2, X ** 2 + Y ** 2 > rad_inf ** 2))]
    offset_avg = np.mean(offset)

    # Return the phase with the average offset removed.
    return phase - offset_avg


def separate_phase_components(forward, reverse):
    """
    Separate the magnetic and electrostatic phases using specimen flipping (time reversal).

    Args:
        forward (ndarray): The phase map in the forward direction (the direction whose phases
                           one is trying to extract).
        reverse (ndarray): The phase map in the reverse direction. Note that the reverse phase
                           map must be mirrored in the axis of rotation so that the forward
                           and reverse phases line up.

    Returns:
        elec (ndarray): The electrostatic phase map.
        mag (ndarray): The magnetic phase map.
    """

    elec = (forward + reverse) / 2
    mag = (forward - reverse) / 2
    return elec, mag


def solenoidal_moment(resolution, axis=2):
    # todo: I don't think I got this working yet. Test it, and fix it if needed.
    """
    Generate a solenoidal (rotational) moment distribution.

    Args:
        resolution (sequence of ints): The dimensions of the moment array in pixels.
        axis (int|optional): Axis about which the vectors circulate.

    Returns:
        moment (ndarray): The solenoidal moment distribution.
    """
    vhat = unit_vector_from_axis(axis)
    _resolution = list(resolution)
    del _resolution[axis]

    radial = construct_k_kernel(_resolution, (1, 1))

    moment = np.cross(vhat, radial)

    return moment
