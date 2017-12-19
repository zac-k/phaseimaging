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

    Args:
        image (ndarray): Real space image to be filtered.
        kernel (ndarray): Fourier-space convolution kernel; i.e., a transfer function.
    Returns:
        image (ndarray): The filtered image.
    """
    assert len(image) == len(kernel) and len(image[0]) == len(kernel[0])
    image = fft.fft2(image)
    image = fft.fftshift(image)
    image *= kernel
    image = fft.ifftshift(image)
    image = fft.ifft2(image)
    return image

def intensity_derivative(image_under, image_over, defocus):
    """
    Computes the intensity derivative from out-of-focus images using a CDA.

    Args:
        image_under (ndarray): The underfocussed image
        image_over (ndarray): The overfocussed image
        defocus (float): The defocus distance, whose magnitude is the same for both images.
    Returns:
        derivative (ndarray): The longitudinal intensity derivative.
    """
    return (image_under - image_over) / (2 * defocus)


def dot_fields(field1, field2):
    """
    Pointwise dot-product of two two-dimentional vector fields.

    Args:
        field1 (ndarray): Two-dimensional vector field.
        field2 (ndarray): A second two-dimensional vector field.
    Returns:
        dot_product (ndarray): Two-dimensional scalar field.
    """
    return field1[:, :, 0] * field2[:, :, 0] + field1[:, :, 1] * field2[:, :, 1]


def accel_volt_to_lambda(accel_volt):
    """
    Calculates wavelength in metres from accelerating voltage in volts.

    Args:
        accel_volt (float): The electron accelerating voltage in volts
    Returns:
        lambda (float): The relativistic wavelength of the electrons in metres.
    """
    return 12.25e-10/np.sqrt(accel_volt) * 1/np.sqrt(1 + (-e * accel_volt / (2 * m0 * c*c)))

def lambda_to_accel_volt(wavlen):
    """
    Calculates electron accelerating voltage in volts from wavelength in metres.

    Args:
        wavlen (float): The relativistic wavelength (lambda) of the electrons.
    Returns:
        accel_volt (float): The electron accelerating voltage in Volts.

    """
    return m0*c*c*(np.sqrt(1 + 2 * (12.25e-10*12.25e-10 * -e)/(wavlen*wavlen * m0 * c*c)) - 1) / -e


def import_specimen(specimen_file):
    """
    Load a specimen array from file.

    Args:
        specimen_file (str): Path of specimen file.
    returns:
        specimen (ndarray): Three-dimensional binary specimen
        mask. A value of 1 indicates a voxel where
        the specimen exists, and a value of 0
        indicates a voxel where it does not.
    """

    # Check if the file has a .npy extension.
    if specimen_file[-4:] == '.npy':
        # Load specimen from numpy file and convert boolean to float.
        specimen = np.load(specimen_file).astype('float64')
    else:
        # Assume the file is in text format. Open the file and read
        # the contents into the output array. Text format specimens must be cubic
        with open(specimen_file, 'r') as f:
            specimen_size = int(len(f.readline().split()))
            resolution = (specimen_size, specimen_size, specimen_size)
            f.seek(0)
            specimen = np.zeros((specimen_size, specimen_size, specimen_size))
            for k in range(resolution[0]):
                for i in range(resolution[1]):
                    specimen[i, :, k] = f.readline().split()
                f.readline()
    return specimen


def import_moment(moment_file):
    """
    Load a moment distribution array from file.

    Args:
        moment_file (str): Path of moment file.
    Returns:
        moment (ndarray): Three-dimensional vector field describing the magnetisation direction
                          and magnitude at each point in space.
    """

    if moment_file[-4:] == '.npy':
        moment = np.load(moment_file)
    else:
        Exception("Only .npy moment files allowed!")
    return moment


def normalised_rms_error(exact, reconstructed, display=False, name=None):
    """
    Calculate normalised rms error between exact and reconstructed
    arrays of arbitrary (but equal) dimensions

    Args:
        exact (ndarray): The exact field.
        reconstructed (ndarray): The reconstructed field.
        display (bool|optional): If `True`, result will also be printed
                                 to the console. Default is False.
        name (string|optional): The name will be printed along with the
                                error if `display` is true.
    Returns:
        error (float): Normalised rms error in the reconstructed field
                       relative to the exact field.
    """

    # todo: refactor this code into a more numpy-friendly (vectorised) form
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

    if display:
        if name is not None:
            text = "in " + name + " "
        else:
            text = ""
        print("Normalised RMS error " + text + "= {0: .1%}".format(error))
    return error


def rotate_vector(v, angle, axis=0):
    """
    Rotates a vector about a given Cartesian axis, by a given angle.

    Args:
        v (ndarray|sequence of floats): The vector to be rotated.
        angle (float): The angle of rotation in radians.
        axis (int|optional): The axis to be rotated about. Default is 0.
    Returns:
        output (ndarray|sequence of floats): The rotated vector.
    """
    cos = np.cos(angle)
    sin = np.sin(angle)
    if axis == 0:
        R = np.array([[1, 0, 0], [0, cos, -sin], [0, sin, cos]])
    elif axis == 1:
        R = np.array([[cos, 0, sin], [0, 1, 0], [-sin, 0, cos]])
    elif axis == 2:
        R = np.array([[cos, -sin, 0], [sin, cos, 0], [0, 0, 1]])

    return np.dot(v, R)

def unit_vector_from_axis(axis):
    """
    Generates a unit vector pointing in the positive direction along a given axis.

    Args:
        axis (int): 0, 1, or 2, representing the axis direction.

    Returns:
        unit_vector (tuple of ints): The unit vector pointing in the positive direction
                                     of the axis.
    """
    assert axis in range(3)
    unit_vector = [0, 0, 0]
    unit_vector[axis] = 1
    return unit_vector

def vec_isin(a, b):
    """
    Given two NxN arrays, determine which rows of a exist somewhere in b.

    Args:
        a (ndarray): Array of elements.
        b (ndarray): Test array; the vectors against which to test each row of a.

    Returns:
        vec_isin (ndarray of bools): Vector of length equal to the number of rows in a,
                                     equal to True if that row is in b, or False if it
                                     is not.
    """
    # Convert each row of both arrays to strings
    a_strings = [','.join(x.astype(str)) for x in a]
    b_strings = [','.join(x.astype(str)) for x in b]

    # Return boolean vector of which strings from a are in b.
    return np.isin(a_strings, b_strings)

