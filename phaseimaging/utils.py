from numpy import fft
import numpy as np

m0 = 9.11e-31  # Electron rest mass
e = -1.6e-19  # Electron charge
c = 3e8  # Speed of light
PI = np.pi

def convolve(image, kernel):
    """
    Filters an image in Fourier space. The kernel here is already in
    Fourier space. Use scipy.signal.fftconvolve() for convolving
    two real space images.
    """
    assert len(image) == len(kernel) and len(image[0]) == len(kernel[0])
    image = fft.fft2(image)
    image = fft.fftshift(image)
    image *= kernel
    image = fft.ifftshift(image)
    image = fft.ifft2(image)
    return image

def intensity_derivative(image_under, image_over, defocus):
    return (image_under - image_over) / (2 * defocus)


def dot_fields(field1, field2):
    """
    Pointwise dot-product of two two-dimentional vector fields.
    :param field1:
    :param field2:
    :return:
    """
    return field1[:, :, 0] * field2[:, :, 0] + field1[:, :, 1] * field2[:, :, 1]


def accel_volt_to_lambda(accel_volt):
    """
    Calculates wavelength in metres from accelerating voltage in volts.

    """
    return 12.25e-10/np.sqrt(accel_volt) * 1/np.sqrt(1 + (-e * accel_volt / (2 * m0 * c*c)))

def lambda_to_accel_volt(wavlen):
    """
    Calculates electron accelerating voltage in volts from wavelength in metres.
    """
    return m0*c*c*(np.sqrt(1 + 2 * (12.25e-10*12.25e-10 * -e)/(wavlen*wavlen * m0 * c*c)) - 1) / -e


def import_specimen(specimen_file):
    """

    Args:
        specimen_file (str): Path of specimen file.
    returns:
        ndarray: Three-dimensional binary specimen
        mask. A value of 1 indicates a voxel where
        the specimen exists, and a value of 0
        indicates a voxel where it does not.
    """
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


def normalised_rms_error(exact, reconstructed, display=False, name=None):

    """
    Calculate normalised rms error between exact and reconstructed
    arrays of arbitrary (but equal) dimensions
    :param exact:
    :param reconstructed:
    :return:
    """

    assert np.shape(exact) == np.shape(reconstructed)
    total = 0
    norm = 0

    len_flat = np.prod(np.shape(exact), dtype=int)
    exact = exact.reshape(len_flat)
    reconstructed = reconstructed.reshape(len_flat)

    for i in range(len_flat):
        total += (exact[i] - reconstructed[i]) * (exact[i] - reconstructed[i])
        norm += exact[i] * exact[i]
    error = float(np.real(np.sqrt(total / norm)))
    if name is not None:
        text = "in " + name + " "
    else:
        text = ""
    if display:
        print("Normalised RMS error " + text + "= {0: .1%}".format(error))
    return error


def rotate_vector(v, angle, axis=0):
    cos = np.cos(angle)
    sin = np.sin(angle)
    if axis == 0:
        R = np.array([[1, 0, 0], [0, cos, -sin], [0, sin, cos]])
    elif axis == 1:
        R = np.array([[cos, 0, sin], [0, 1, 0], [-sin, 0, cos]])
    elif axis == 2:
        R = np.array([[cos, -sin, 0], [sin, cos, 0], [0, 0, 1]])

    return np.dot(v, R)
