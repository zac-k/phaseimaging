from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import warnings


def plot_image(image, limits=None, title=None):
    """
    Plot an image array on the screen.

    Args:
        image (ndarray): The image data to plot.
        limits (sequence of floats|optional): Two-element sequence providing the minimum
                                              and maximum values of the plot. Values outside
                                              this range will be clipped at the min/max. If
                                              limits are not provided, the color axis will be
                                              scaled according to the max/min values in the
                                              image array.
        title (string|optional): The title to display on the figure. Default is None.

    Returns:
        None
    """

    # If the image is complex, warn the user and plot only the real part.
    if np.iscomplexobj(image):
        warnings.warn('Image is complex; only displaying real part')
        image = np.real(image)

    # Plot the scaled image if no limits are provided.
    if limits is None:
        plt.imshow(image, cmap='gray')
    # Otherwise, plot the image using the limits.
    else:
        vmin = limits[0]
        vmax = limits[1]
        plt.imshow(image, cmap='gray', vmin=vmin, vmax=vmax)
    # Add the title if one is provided.
    if title is not None:
        plt.title(title)
    plt.show(block=True)


def plot_quiver(vec_field, width=(1, 1, 1), title=None):
    """
    Plot a three-dimensional quiver plot on screen.

    Args:
        vec_field (ndarray): The vector field to plot.
        width (sequence of floats|optional): The dimensions of the field in units of length.
                                             Default is the unit cube.
        title (string|optional): Title to display on the plot. Default is None.

    Returns:
        None
    """

    # Generate a figure with 3D axes.
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Get resolution
    res = np.shape(vec_field)[0:3]

    # Make the voxel grid
    x, y, z = np.meshgrid(np.arange(-width[0] / 2, width[0] / 2, width[0] / res[0]),
                          np.arange(-width[1] / 2, width[1] / 2, width[1] / res[1]),
                          np.arange(-width[2] / 2, width[2] / 2, width[2] / res[2]))
    # todo: fix scaling in quiver
    ax.quiver(x, y, z, vec_field[:, :, :, 0], vec_field[:, :, :, 1], vec_field[:, :, :, 2],
              length=0.01 * np.linalg.norm(width), normalize=True)

    # Add the title if one is provided.
    if title is not None:
        plt.title(title)
    plt.show(block=True)


def save_image(image, output_path, limits=None):
    """
    Save an image array in an image format. This produces a clean image with no axes or border.
    Args:
        image (ndarray): The image data to save.
        output_path (string): The path of file to save.
        limits (sequence of floats|optional): Two-element sequence providing the minimum
                                              and maximum values of the image. Values outside
                                              this range will be clipped at the min/max. If
                                              limits are not provided, the color axis will be
                                              scaled according to the max/min values in the
                                              image array.

    Returns:
        None
    """

    # Generate a figure with no axes, border, etc.
    fig = plt.figure()
    fig.set_size_inches(1, 1)
    ax = plt.Axes(fig, [0, 0, 1, 1])
    ax.set_axis_off()
    fig.add_axes(ax)

    # If the image is complex, warn the user and discard the imaginary part.
    if np.iscomplexobj(image):
        warnings.warn('Image is complex; only saving real part')
        image = np.real(image)

    # Plot the image, scaled according to the limits if they are provided.
    if limits is None:
        ax.imshow(image, cmap='gray', aspect='auto', interpolation='none')
    else:
        vmin = limits[0]
        vmax = limits[1]
        ax.imshow(image, cmap='gray', vmin=vmin, vmax=vmax, aspect='auto', interpolation='none')

    # Save the figure with a high resolution, and close it.
    plt.savefig(output_path, dpi=800)
    plt.close()