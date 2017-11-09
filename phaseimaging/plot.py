import matplotlib.pyplot as plt
import numpy as np


def plot_image(image, type=None):

    if np.iscomplexobj(image):
        print('Image is complex; only displaying real part')
        image = np.real(image)
    if type == 'image':
        vmin = 0
        vmax = 2
        plt.imshow(image, cmap='gray', vmin=vmin, vmax=vmax)
    elif type == 'phase':
        vmin = -3
        vmax = 3
        plt.imshow(image, cmap='gray', vmin=vmin, vmax=vmax)
    else:
        plt.imshow(image, cmap='gray')
    plt.show(block=True)


def save_image(image, output_path, limits=None):

    fig = plt.figure()
    fig.set_size_inches(1, 1)
    ax = plt.Axes(fig, [0, 0, 1, 1])
    ax.set_axis_off()
    fig.add_axes(ax)
    if np.iscomplexobj(image):
        print('Image is complex; only saving real part')
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