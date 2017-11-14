### phaseimaging.save_image

Similar to `plot_image`, but saves the visualisation in an image format rather than plotting it on-screen.

    phaseimaging.save_image(image, output_path, limits=None)

|  |  |  |
|---|---|---|
| Parameters: | **image** : *ndarray* |   |
|  |  | The array to be saved as an image. For complex arrays, only the real part will be kept, with a warning displayed on the command line.  |
|   | **output_path** : *string* |   |
|  |  | The output path of the file to be saved. The extension will determine the format (e.g., 'png', 'eps', etc.).  |
|  | **limits** : *tuple*, *list*, *optional* |  |
|  |  | Two element tuple or list containing the minimum and maximum values to be displayed. Intensities outside these limits will be clipped. If omitted, grayscale will be scaled according to the maximum and minimum values of the image array.

[Back to Functions](functions.md)