import phaseimaging as phim

# Import the specimen from file
specimen = phim.import_specimen('C:/Users/zac/PycharmProjects/phaseimaging/specimen')

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
image_over = phim.transfer_image(8e-6, 1.96e-12, (100e-9, 100e-9), phase)

# Add noise to images
image_under = phim.add_noise(image_under, 1, 0.15)
image_over = phim.add_noise(image_over, 1, 0.15)

# Compute the longitudinal derivative of the intensity
derivative = phim.intensity_derivative(image_under, image_over, 8e-6)

# Retrieve the phase
phase_ret = phim.retrieve_phase_tie(1.96e-12, (100e-9, 100e-9), derivative)

# Plot the phases and intensity images
phim.plot_image(phase, limits=[-3, 3])
phim.plot_image(image_under, limits=[0, 2])
phim.plot_image(image_over, limits=[0, 2])
phim.plot_image(phase_ret, limits=[-3, 3])
