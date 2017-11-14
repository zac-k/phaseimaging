import phaseimaging as phim
# Import the specimen from file
specimen = phim.import_specimen('../../specimen')

# Set the physical conditions
accelerating_voltage = 300e3
mean_inner_potential = -17 + 1j
specimen_array_dimensions = (100e-9, 100e-9, 100e-9)
mhat = (1, 0, 0)

# Project the electrostatic phase
phase_elec = phim.project_electrostatic_phase(specimen,
                                              accelerating_voltage,
                                              mean_inner_potential,
                                              specimen_array_dimensions)

# Set magnetisation strength
mass_mag = 80  # emu/g
density = 5.18  # g/cm^3
magnetisation = mass_mag * density * 1000  # A/m

# Project magnetic phase using a magnetisation vector in the x-direction
phase_mag = phim.project_magnetic_phase(specimen, mhat, magnetisation, specimen_array_dimensions)

# Combine the two components of the phase
phase = phase_mag + phase_elec

# Set the imaging parameters
wavlen = phim.accel_volt_to_lambda(accelerating_voltage)
image_dimensions = specimen_array_dimensions[0:2]
defocus = 8e-6
noise_level = 0.15

# Transfer to under- and over-focus planes
image_under = phim.transfer_image(-defocus, wavlen, image_dimensions, phase)
image_in = phim.transfer_image(0, wavlen, image_dimensions, phase)
image_over = phim.transfer_image(defocus, wavlen, image_dimensions, phase)

# Add noise to images and apply apodisation
image_under = phim.apodise(phim.add_noise(image_under, 1, noise_level))
image_in = phim.apodise(phim.add_noise(image_in, 1, noise_level))
image_over = phim.apodise(phim.add_noise(image_over, 1, noise_level))

# Compute the longitudinal derivative of the intensity
derivative = phim.intensity_derivative(image_under, image_over, defocus)

# Retrieve the phase
phase_ret = phim.retrieve_phase_tie(wavlen, image_dimensions, derivative, image_in)
phase_ret = phim.remove_offset(phase_ret)
phim.normalised_rms_error(phase, phase_ret, display=True)

# Plot the phases and intensity images
phim.plot_image(phase, limits=[-3, 3])
phim.plot_image(image_under, limits=[0, 2])
phim.plot_image(image_in, limits=[0, 2])
phim.plot_image(image_over, limits=[0, 2])
phim.plot_image(phase_ret, limits=[-3, 3])
