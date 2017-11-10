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
                         specimen_file='C:/Users/zac/PycharmProjects/phaseimaging/specimen')

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

# Use aliases for intensities, for brevity
image_under = through_focal_series.intensities[0]
image_over = through_focal_series.intensities[-1]

# Retrieve phase
phase_ret = phim.Phase(resolution=specimen.resolution[0:2], width=specimen.width[0:2])
phase_ret.retrieve_phase_tie(through_focal_series, beam)

# Plot intensities and phases
phase.plot(limits=[-3, 3])
image_under.plot(limits=[0, 2])
image_over.plot(limits=[0, 2])
phase_ret.plot(limits=[-3, 3])
