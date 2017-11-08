import numpy as np
from .utils import convolve
from numpy import pi as PI


def _construct_k_squared_kernel(resolution, image_width):

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
        k_squared_kernel = _construct_k_squared_kernel(resolution, image_width)
    return np.exp(1j * (PI * defocus * wavelength * k_squared_kernel))


def transfer_image(defocus, wavelength, image_width, phase, is_image):
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