# phaseimaging

# Index

* [Introduction](#introduction)
* [Functions](functions/functions)
* [Classes](classes/classes)
* [Sample Implementations](sample_implementations/sample_implementations.md)

## Introduction

This library provides a variety of functions related to phase imaging. Additionally, classes and their associated methods are included for if an object-oriented approach is preferred. Some of the included functions and methods can be used to process experimental images, but I have also included functions to simulate these images. For experimental data, pre-processing will be required for meaningful results. In the future, I will add image registration and other preprocessing functions to enable perfoming all steps after image acquisition using only this library.

Note that although I have coded most of the functions to allow for rectangular images and specimen arrays, most of my testing has used square and cubic arrays, so there is a potential for bugs if you are using rectangular arrays. Please let me know if you experience any issues with any of the functions.

This documentation is a work in progress, so newer parts of the code may not always be documented here. Check the documentation in the source code if any explanation is missing from here, and if it's not there either, feel free to contact me.




### Procedural style
```python
import phaseimaging as phim

# Import the specimen from file
specimen = phim.import_specimen('../specimen')

# Project the electrostatic phase
phase_elec = phim.project_electrostatic_phase(specimen, 300e3, -17 + 1j, (100e-9, 100e-9, 100e-9))

# Set magnetisation strength
mass_mag = 80  # emu/g
density = 5.18  # g/cm^3
magnetisation = mass_mag * density * 1000  # A/m

# Project magnetic phase using a magnetisation vector in the x-direction
phase_mag = phim.project_magnetic_phase(specimen,
                                    (1,0,0),
                           magnetisation,
                                    (100e-9, 100e-9, 100e-9))

# Combine the two components of the phase
phase = phase_mag + phase_elec

# Transfer to under- and over-focus planes
image_under = phim.transfer_image(-8e-6, 1.96e-12, (100e-9, 100e-9), phase)
image_in = phim.transfer_image(0, 1.96e-12, (100e-9, 100e-9), phase)
image_over = phim.transfer_image(8e-6, 1.96e-12, (100e-9, 100e-9), phase)

# Add noise to images and apply apodisation
image_under = phim.apodise(phim.add_noise(image_under, 1, 0.15))
image_in = phim.apodise(phim.add_noise(image_in, 1, 0.15))
image_over = phim.apodise(phim.add_noise(image_over, 1, 0.15))

# Compute the longitudinal derivative of the intensity
derivative = phim.intensity_derivative(image_under, image_over, 8e-6)

# Retrieve the phase
phase_ret = phim.retrieve_phase_tie(1.96e-12, (100e-9, 100e-9), derivative, image_in)
phase_ret = phim.remove_offset(phase_ret)
phim.normalised_rms_error(phase, phase_ret, display=True)

# Plot the phases and intensity images
phim.plot_image(phase, limits=[-3, 3])
phim.plot_image(image_under, limits=[0, 2])
phim.plot_image(image_in, limits=[0, 2])
phim.plot_image(image_over, limits=[0, 2])
phim.plot_image(phase_ret, limits=[-3, 3])
```
### Object-oriented style

```python
import phaseimaging as phim

# Set magnetisation strength
mass_mag = 80  # emu/g
density = 5.18  # g/cm^3
magnetisation = mass_mag * density * 1000  # A/m

# Build specimen
specimen = phim.Specimen(width = (100e-9, 100e-9, 100e-9),
                         mean_inner_potential=-17+1j,
                         magnetisation=magnetisation,
                         mhat=(1, 0, 0),
                         specimen_file='../specimen')

# Initialise phase and beam
phase = phim.Phase(resolution=specimen.resolution[0:2], width=specimen.width[0:2])
beam = phim.Beam(phim.accel_volt_to_lambda(300e3))

# Project the electrostatic and magnetic phases
phase.project_electrostatic(specimen, beam)
phase.project_magnetic(specimen, beam)

# Generate through-focal series of intensities and compute derivative
through_focal_series = phim.ThroughFocalSeries(phase.resolution, phase.width, [-8e-6, 0, 8e-6])
through_focal_series.transfer_images(phase, beam)
through_focal_series.add_noise(0.05)
through_focal_series.compute_derivative()

# Apodise images in through-focal series
through_focal_series.apodise()

# Use aliases for intensities, for brevity
image_under = through_focal_series.intensities[0]
image_in = through_focal_series.intensities[1]
image_over = through_focal_series.intensities[-1]

# Retrieve phase
phase_ret = phim.Phase(resolution=specimen.resolution[0:2], width=specimen.width[0:2])
phase_ret.retrieve_phase_tie(through_focal_series, beam)
phase.remove_offset()
phase_ret.normalised_rms_error(phase, display=True)

# Plot intensities and phases
phase.plot(limits=[-3, 3])
image_under.plot(limits=[0, 2])
image_in.plot(limits=[0, 2])
image_over.plot(limits=[0, 2])
phase_ret.plot(limits=[-3, 3])
```

    
### Output

    Normalised RMS error =  28.1%
    
![projected phase](README/phase.png) ![under-focus image](README/image_under.png) ![in-focus image](README/image_in.png) ![over-focus image](README/image_over.png) ![retrieved phase](README/phase_ret.png)