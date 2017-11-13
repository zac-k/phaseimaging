### *class* phaseimaging.Intensity

| Attributes |  |  
|---|---|
| `resolution` | Dimensions of the image in pixels. |
| `shape` | Alias of `resolution`. |
| `width` | Dimensions of the image in units of length. |
| `image` | Numpy array containing the image data. |
| `defocus` | Distance corresponding to the defocus of the intensity measurement. |
| `incident` | Uniform incident intensity of the beam. |

| Methods |  |  
|---|---|
| `plot([limits])` | Generates an on-screen image using matplotlib. |
| `apodise([rad_sup])` | Apodises the image with a circular window function of radius rad_sup (as a fraction of the image width). |
| `add_noise(sigma)` | Adds Poisson noise to the images with a fractional noise level of sigma. |
| `transfer(phase, beam)` | Generates an image by transferring the phase of the beam to an in-, or out-of-focus, plane. |

[Back to Classes](classes.md)   