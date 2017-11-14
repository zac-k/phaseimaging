### phaseimaging.intensity_derivative

Compute longitudinal derivative of the intensity.

    phaseimaging.intensity_derivative(image_under, image_over, defocus)

|  |  |  |
|---|---|---|
| Parameters: | **image_under** : *array-like* |   |
|  |  | Intensity at the under-focus plane |
|  | **image_over** : *array-like* |  |
|  |  | Intensity at the over-focus plane |
|  | **defocus** : *float* |  |
|  |  | Defocus distance. This is the size of the defocus, which much be equal for both defocussed intensities |
| Returns: | **derivative** : *ndarray* |  |
|  |  | Computed derivative |

[Back to Functions](functions.md)