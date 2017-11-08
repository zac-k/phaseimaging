import numpy as np
from numpy import fft
from .utils import convolve
PI = np.pi
from .plot import save_image
from .plot import plot_image
from .imaging import transfer_image


def intensity_derivative(image_under, image_over, defocus):
    return (image_under - image_over) / (2 * defocus)





def _construct_inverse_k_squared_kernel(resolution, image_width, reg_param_tie):

    kernel = np.zeros(resolution, dtype=complex)
    da = tuple([1 / x for x in image_width])
    for i in range(resolution[0]):
        for j in range(resolution[1]):
            i0 = i - resolution[0] / 2
            j0 = j - resolution[1] / 2

            kernel[i, j] = (i0 * i0 * da[0] * da[0]) + (j0 * j0 * da[1] * da[1])
    kernel = kernel / (kernel * kernel + reg_param_tie ** 4)
    return kernel


def _construct_k_kernel(resolution, image_width):
    kernel = np.zeros(list((resolution[0], resolution[1], 2)), dtype=complex)
    da = tuple([1 / x for x in image_width])
    for i in range(resolution[0]):
        for j in range(resolution[1]):
            i0 = i - resolution[0] / 2
            j0 = j - resolution[1] / 2
            kernel[i, j, 0] = i0 * da[0]
            kernel[i, j, 1] = j0 * da[1]
    return kernel

def dot_fields(field1, field2):
    """
    Pointwise dot-product of two two-dimentional vector fields.
    :param field1:
    :param field2:
    :return:
    """
    return field1[:, :, 0] * field2[:, :, 0] + field1[:, :, 1] * field2[:, :, 1]


def retrieve_phase_tie(defocus,
                       wavelength,
                       image_width,
                       image_under,
                       image_over,
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
    assert image_under.shape == image_over.shape
    if image_in is not None:
        assert image_in.shape == image_over.shape
    resolution = image_under.shape
    assert len(image_under.shape) == 2
    assert len(image_width) == 2

    if reg_param_tie is None:
        reg_param_tie = 0.1 / (np.mean(image_width) * np.mean(resolution))

    if k_kernel is None:
        k_kernel = _construct_k_kernel(resolution, image_width)
    if inverse_k_squared_kernel is None:
        inverse_k_squared_kernel = _construct_inverse_k_squared_kernel(resolution, image_width, reg_param_tie)

    derivative = intensity_derivative(image_under, image_over, defocus)
    if image_in is not None:
        regularised_inverse_intensity = image_in / (image_in * image_in + reg_param * reg_param)
        prefactor = (1. / wavelength) / (2. * PI)
        derivative_vec = np.zeros((resolution[0], resolution[1], 2), dtype=complex)
        derivative_vec[:, :, 0] = convolve(derivative, k_kernel[:, :, 0] *
                                           inverse_k_squared_kernel
                                           ) * regularised_inverse_intensity
        derivative_vec[:, :, 1] = convolve(derivative, k_kernel[:, :, 1] *
                                           inverse_k_squared_kernel
                                           ) * regularised_inverse_intensity

        derivative_vec[:, :, 0] = fft.fftshift(fft.fft2(derivative_vec[:, :, 0]))
        derivative_vec[:, :, 1] = fft.fftshift(fft.fft2(derivative_vec[:, :, 1]))
        derivative = dot_fields(k_kernel, derivative_vec) * inverse_k_squared_kernel
        phase_retrieved = prefactor * fft.ifft2(fft.ifftshift(derivative))

    else:
        prefactor = (1. / wavelength) / (2 * PI * image_intensity)

        filtered = convolve(derivative, inverse_k_squared_kernel)
        phase_retrieved = prefactor * filtered
    phase_retrieved = np.real(phase_retrieved)
    return phase_retrieved


def project_electrostatic_phase(specimen, accel_volt, mean_inner_potential, image_width):
    """
    Computes the image plane phase shift of the specimen.

    """
    resolution = specimen.shape
    assert len(resolution) == 3


    if len(image_width) == 3:
        image_width = image_width[0]
    wavelength = accel_volt_to_lambda(accel_volt)
    dz = image_width[0] / resolution[0]
    return np.sum(specimen, axis=0) * PI/(accel_volt * wavelength) * mean_inner_potential * dz


def accel_volt_to_lambda(accel_volt):
    """
    Calculates wavelength in metres from accelerating voltage in volts.

    """
    return 12.25e-10/np.sqrt(accel_volt) * 1/np.sqrt(1+(1.6e-19 * accel_volt/(2*9.11e-31*9e16)))


def import_specimen(specimen_file):

    with open(specimen_file, 'r') as f:
        # todo: modify the next two lines to allow reading rectangular specimens
        specimen_size = int(len(f.readline().split()))
        resolution = (specimen_size, specimen_size, specimen_size)
        f.seek(0)
        specimen = np.zeros((specimen_size, specimen_size, specimen_size))
        for k in range(resolution[0]):
            for i in range(resolution[1]):
                specimen[i, :, k] = f.readline().split()
            f.readline()
    return specimen
