import numpy as np
from .utils import *
from numpy import pi as PI
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

    kernel = np.zeros(resolution, dtype=complex)
    for i in range(resolution[0]):
        for j in range(resolution[1]):
            i0 = i - resolution[0] / 2
            j0 = j - resolution[1] / 2
            da = [1 / x for x in image_width]
            kernel[i, j] = (i0 * i0) * da[0] * da[0] + j0 * j0 * da[1] * da[1]

    return kernel

def _set_transfer_function(defocus, wavelength, resolution, image_width, k_squared_kernel=None):
    """
    Sets the imaging system's transfer function for a given
    defocus.
    """
    if k_squared_kernel is None:
        k_squared_kernel = construct_k_squared_kernel(resolution, image_width)
    return np.exp(1j * (PI * defocus * wavelength * k_squared_kernel))


def transfer_image(defocus, wavelength, image_width, phase, is_image=True):
    """
    Uses the defocus exact_phase (at the image plane) to produce an out of focus
    image or wavefield.
    """
    transfer_function = _set_transfer_function(defocus,
                                               wavelength,
                                               resolution=phase.shape,
                                               image_width=image_width)
    wavefunction = np.exp(1j*phase)
    wavefunction = convolve(wavefunction, transfer_function)
    if is_image:
        return np.absolute(wavefunction) * np.absolute(wavefunction)
    else:
        return wavefunction


def construct_inverse_k_squared_kernel(resolution, image_width, reg_param=None):
    if reg_param is None:
        reg_param = 0.1 / (np.mean(image_width) * np.mean(resolution))
    kernel = np.zeros(resolution, dtype=complex)
    da = tuple([1 / x for x in image_width])
    for i in range(resolution[0]):
        for j in range(resolution[1]):
            i0 = i - resolution[0] / 2
            j0 = j - resolution[1] / 2

            kernel[i, j] = (i0 * i0 * da[0] * da[0]) + (j0 * j0 * da[1] * da[1])
    kernel = kernel / (kernel * kernel + reg_param ** 4)
    return kernel


def construct_k_kernel(resolution, image_width):
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
    the micrographs.
    :return:
    """

    if image_in is not None:
        assert image_in.shape == intensity_derivative.shape
    resolution = intensity_derivative.shape
    assert len(intensity_derivative.shape) == 2
    assert len(image_width) == 2



    if k_kernel is None:
        k_kernel = construct_k_kernel(resolution, image_width)
    if inverse_k_squared_kernel is None:
        inverse_k_squared_kernel = construct_inverse_k_squared_kernel(resolution, image_width, reg_param_tie)


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

    else:
        prefactor = (1. / wavelength) / (2 * PI * image_intensity)

        filtered = convolve(intensity_derivative, inverse_k_squared_kernel)
        phase_retrieved = prefactor * filtered
    phase_retrieved = np.real(phase_retrieved)
    return phase_retrieved


def project_electrostatic_phase(specimen, accel_volt, mean_inner_potential, image_width):
    """
    Computes the image plane phase shift of the specimen.

    """
    resolution = specimen.shape
    assert len(resolution) == 3
    wavelength = accel_volt_to_lambda(accel_volt)
    dz = image_width[2] / resolution[2]
    return np.sum(specimen, axis=2) * PI/(accel_volt * wavelength) * mean_inner_potential * dz


def project_magnetic_phase(specimen,
                           mhat,
                           magnetisation,
                           image_width,
                           k_kernel=None,
                           inverse_k_squared_kernel=None):
    resolution = specimen.shape
    assert len(image_width) == 3
    if k_kernel is None:
        k_kernel = construct_k_kernel(resolution, image_width)
    if inverse_k_squared_kernel is None:
        inverse_k_squared_kernel = construct_inverse_k_squared_kernel(resolution[0:2], image_width[0:2])
    mhat = mhat / np.linalg.norm(mhat)

    mhatcrossk_z = np.cross(mhat, k_kernel)[:, :, 2]
    D0 = fft.fftshift(fft.fftn(specimen))[:, :, int(resolution[2]/2)] * (image_width[2] / resolution[2])

    phase = (1j * PI * mu0 * magnetisation / phi0) * D0 * inverse_k_squared_kernel * mhatcrossk_z
    phase = fft.ifftshift(phase)
    return np.real(fft.ifft2(phase))


def add_noise(image, i_in, sigma):

    return np.where(image >= 0, np.random.poisson(image / (sigma*sigma * i_in)) * (sigma*sigma * i_in), 0)


def apodise(image, rad_sup=0.5):

    x = np.arange(int(-image.shape[0] / 2), int(image.shape[0] / 2))
    y = np.arange(int(-image.shape[1] / 2), int(image.shape[1] / 2))
    X,Y = np.meshgrid(x / image.shape[0], y / image.shape[1])
    output = np.where(X**2 + Y**2 < rad_sup**2, image, 0)
    return output


def remove_offset(phase, rad_inf=0.45, rad_sup=0.5):
    x = np.arange(int(-phase.shape[0] / 2), int(phase.shape[0] / 2))
    y = np.arange(int(-phase.shape[1] / 2), int(phase.shape[1] / 2))
    X, Y = np.meshgrid(x / phase.shape[0], y / phase.shape[1])
    offset = phase[np.where(np.logical_and(rad_sup ** 2 > X ** 2 + Y ** 2, X ** 2 + Y ** 2 > rad_inf ** 2))]
    offset_avg = np.mean(offset)
    return phase - offset_avg


def separate_phase_components(forward, reverse):
    elec = (forward + reverse) / 2
    mag = (forward - reverse) / 2
    return elec, mag
