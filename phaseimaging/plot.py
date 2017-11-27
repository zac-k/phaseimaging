from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
import warnings


def plot_image(image, limits=None, title=None):

    if np.iscomplexobj(image):
        warnings.warn('Image is complex; only displaying real part')
        image = np.real(image)
    if limits is None:
        plt.imshow(image, cmap='gray')
    else:
        vmin = limits[0]
        vmax = limits[1]
        plt.imshow(image, cmap='gray', vmin=vmin, vmax=vmax)
    if title is None:
        title = ""
    plt.title(title)
    plt.show(block=True)


def plot_quiver(vec_field, width=(1, 1, 1), title=None):

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

    if title is None:
        title = ""
    plt.title(title)
    plt.show(block=True)


def save_image(image, output_path, limits=None):

    fig = plt.figure()
    fig.set_size_inches(1, 1)
    ax = plt.Axes(fig, [0, 0, 1, 1])
    ax.set_axis_off()
    fig.add_axes(ax)
    if np.iscomplexobj(image):
        warnings.warn('Image is complex; only saving real part')
        image = np.real(image)
    if limits is None:
        ax.imshow(image, cmap='gray', aspect='auto', interpolation='none')
    else:
        vmin = limits[0]
        vmax = limits[1]
        ax.imshow(image, cmap='gray', vmin=vmin, vmax=vmax, aspect='auto', interpolation='none')
    # Remove tick marks from plot

    plt.savefig(output_path, dpi=800)
    plt.close()