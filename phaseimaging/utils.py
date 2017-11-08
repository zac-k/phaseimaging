from numpy import fft

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