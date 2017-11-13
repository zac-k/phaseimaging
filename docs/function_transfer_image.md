### phaseimaging.transfer_image

Generate in- or out-of-focus intensity or wavefield from exit phase.

    phaseimaging.transfer_image(defocus, wavelength, image_width, phase, is_image=True)

|  |  |  |
|---|---|---|
| Parameters: | **defocus** : *float* |   |
|  |  | Microscope defocus in metres. Positive for over-focus, negative for under-focus, and zero for in-focus. |
|  | **wavelength** : *float* |  |
|  |  | Wavelength of incident beam in metres |
|  | **image_width** : *tuple*, *list* |  |
|  |  | Three element tuple or list containing the width of the specimen array in metres, in the x-, y-, and z-direction, respectively. |
|  | **phase** : *ndarray* |  |
|  |  | Exit phase of the wavefield |
|  | **is_image** : *bool*, *optional* |  |
|  |  | If true, return will be the intensity at the in- or out-of-focus plane. If false, return will be the complex wavefield at the in- or out-of-focus plane. Default is `True`. |
| Returns: | **image** : *ndarray* |  |
|  |  | Computed image. The intensity or wavefield depending on the value of `is_image`. Default is intensity. |
