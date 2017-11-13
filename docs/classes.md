## Classes

### Image

| Attributes |  |  
|---|---|
| `resolution` | Dimensions of the image in pixels. |
| `width` | Dimensions of the image in units of length. |
| `image` | Numpy array containing the image data. |

| Methods |  |  
|---|---|
| `plot([limits])` | Generates an on-screen image using matplotlib. |
| `apodise([rad_sup])` | Apodises the image with a circular window function of radius rad_sup (as a fraction of the image width). |

   

### Intensity

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

   

### Phase

| Attributes |  |  
|---|---|
| `resolution` | Dimensions of the phase map in pixels. |
| `width` | Dimensions of the phase map in units of length. |
| `image` | Numpy array containing the phase. |
| `defocus` | Distance corresponding to the defocus of the phase. Typically zero. |
| `k_kernel` | Kernel of k-vectors for phase retrieval. Generated during the phase retrieval step if it does not already exist. |
| `inverse_k_squared_kernel` | Regularised k^-2 kernel for phase retrieval. Generated during the phase retrieval step if it does not already exist. |

| Methods |  |  
|---|---|
| `plot([limits])` | Generates an on-screen image using matplotlib. |
| `apodise([rad_sup])` | Apodises the image with a circular window function of radius rad_sup (as a fraction of the image width). |
| `construct_kernels([reg_param])` | Generated `k_kernel` and `inverse_k_squared_kernel`. |
| `project_electrostatic(specimen, beam)` | Generates a phase map by projecting the electrostatic potential of the specimen. This is added to the existing phase.|
| `project_magnetic(specimen, beam)` | Generates a phase map by projecting the magnetic vector potential of the specimen. This is added to the existing phase.|
| `retrieve_phase_tie(tfs, beam)` | Retrieves the phase from a through-focal series of intensities. This replaces any pre-existing phase map.|

  
### Wavefield
    
| Attributes |  |  
|---|---|
| `resolution` | Dimensions of the image in pixels. |
| `width` | Dimensions of the image in units of length. |
| `image` | Numpy array containing the wavefield. |
| `defocus` | Distance corresponding to the defocus of the wavefield. Typically zero. |

| Methods |  |  
|---|---|
| `plot([limits])` | Generates an on-screen image using matplotlib. |
| `apodise([rad_sup])` | Apodises the image with a circular window function of radius rad_sup (as a fraction of the image width). |

### Specimen
   
| Attributes |  |  
|---|---|
| `resolution` | Dimensions of the specimen array in pixels. |
| `width` | Dimensions of the specimen array in units of length. |
| `image` | Numpy array containing the specimen. |
| `mean_inner_potential` | Mean inner potential of the specimen. Complex values can be used to simulate attenuation. |
| `magnetisation` | Volume magnetisation. |
| `mhat` | Vector pointing in the magnetisation direction. The length of this vector has no effect. |

| Methods |  |  
|---|---|
| `plot([limits])` | Generates an on-screen image, from a central slice, using matplotlib. |
| `apodise([rad_sup])` | Apodises the with a cylindrical window function of radius rad_sup (as a fraction of the array width). |

### Beam

| Attributes |  |  
|---|---|
| `wavelength` | Wavelength of the beam. |
| `accel_volt` | Accelerating voltage for an electron beam. Automatically generated from the wavelength on initialisation. |


### ThroughFocalSeries
   
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
