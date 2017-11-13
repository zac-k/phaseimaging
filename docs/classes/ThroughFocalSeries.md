### *class* phaseimaging.ThroughFocalSeries
   
| Attributes |  |  
|---|---|
| `intensities` | Array of `Intensity` objects. |
| `derivative` | Longitudinal derivate of the intensity, computed from the through-focal series. |
| `defoci` | Defoci corresponding to the respective intensities. |
| `len` | Number of images in the series. |

| Methods |  |  
|---|---|
| `transfer_images(phase, beam)` | Generates the through focal series corresponding to the phase, using the distances given by `defoci`. |
| `add_noise(sigma)` | Adds Poisson noise to every image in the series, using the fractional noise level sigma. |
| `compute_derivative()` | Computes the longitudinal intensity derivative from the through-focal series. |

[Back to Classes](classes.md)