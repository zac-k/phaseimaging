### phaseimaging.plot_image

Plots a numpy array as an image, for visualising phase, intensity, etc.

    phaseimaging.plot_image(image, limits=None)

|  |  |  |
|---|---|---|
| Parameters: | **image** : *ndarray* |   |
|  |  | The array to be displayed as an image. For complex arrays, only the real part will be displayed, with a warning displayed on the command line.  |
|  | **limits** : *tuple*, *list*, *optional* |  |
|  |  | Two element tuple or list containing the minimum and maximum values to be displayed. Intensities outside these limits will be clipped. If omitted, grayscale will be scaled according to the maximum and minimum values of the image array.

