import matplotlib.pyplot as plt
import numpy as np
import warnings


def plot_image(image, limits=None):

    if np.iscomplexobj(image):
        warnings.warn('Image is complex; only displaying real part')
        image = np.real(image)
    if limits is None:
        plt.imshow(image, cmap='gray')
    else:
        vmin = limits[0]
        vmax = limits[1]
        plt.imshow(image, cmap='gray', vmin=vmin, vmax=vmax)
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