import phaseimaging as phim
import copy

# Set magnetisation strength
mass_mag = 80  # emu/g
density = 5.18  # g/cm^3
magnetisation = mass_mag * density * 1000  # A/m

# Build specimen
specimen = phim.Specimen(width=(100e-9, 100e-9, 100e-9),
                         mean_inner_potential=-17+1j,
                         magnetisation=magnetisation,
                         mhat=(0.707, 0.707, 0),
                         specimen_file='../../specimen')

# Initialise phase and beam
phase = phim.Phase(resolution=specimen.resolution[0:2], width=specimen.width[0:2], name="Total Phase")
beam = phim.Beam(phim.accel_volt_to_lambda(300e3))

# Project the electrostatic and magnetic phases. Store a deep copy
# of the magnetic component for computing the error later.
phase.project_magnetic(specimen, beam)
phase_mag = copy.deepcopy(phase)
phase_mag.name = "Exact Magnetic Phase"
phase.project_electrostatic(specimen, beam)

# Set noise level and defoci.
noise_level = 0.01
defoci = [-8e-6, 0, 8e-6]

# Generate through-focal series of intensities and compute derivative.
through_focal_series = phim.ThroughFocalSeries(phase.resolution, phase.width, defoci=defoci)
through_focal_series.transfer_images(phase, beam)
through_focal_series.add_noise(noise_level)
through_focal_series.compute_derivative()

# Apodise images in through-focal series
through_focal_series.apodise()

# Retrieve the phase and compute the total error.
phase_ret = phim.Phase(resolution=specimen.resolution[0:2], width=specimen.width[0:2], name="Retrieved Phase")
phase_ret.retrieve_phase_tie(through_focal_series, beam)
phase_ret.remove_offset()

# Rotate the specimen for extracting the magnetic phase. Include
# a small error in the x and y rotations. The z rotation can
# be adjusted after projection if desired.
specimen.rotate(angle=178, axis=0)
specimen.rotate(angle=1, axis=1)

# Project the phases in the reverse direction.
phase_reverse = phim.Phase(resolution=specimen.resolution[0:2], width=specimen.width[0:2], name="Total Reverse Phase")
phase_reverse.project_electrostatic(specimen, beam)
phase_reverse.project_magnetic(specimen, beam)

# Generate the through-focal series in the reverse direction and compute the derivative.
through_focal_series_reverse = phim.ThroughFocalSeries(phase.resolution, phase.width, defoci=defoci)
through_focal_series_reverse.transfer_images(phase_reverse, beam)
through_focal_series_reverse.add_noise(noise_level)
through_focal_series_reverse.compute_derivative()

# Retrieve the reverse direction total phase.
phase_ret_reverse = phim.Phase(resolution=specimen.resolution[0:2], width=specimen.width[0:2],
                               name="Retrieved Reverse Phase")
phase_ret_reverse.retrieve_phase_tie(through_focal_series_reverse, beam)
phase_ret_reverse.remove_offset()

# Mirror the reverse phase and use the forward and reverse total phases to obtain
# the electrostatic and magnetic components.
phase_ret_reverse.flip(axis=1)
phase_ret_elec = phim.Phase(resolution=specimen.resolution[0:2], width=specimen.width[0:2],
                            name="Retrieved Electrostatic Phase")
phase_ret_mag = phim.Phase(resolution=specimen.resolution[0:2], width=specimen.width[0:2],
                            name="Retrieved Magnetic Phase")
phase_ret_elec.image, phase_ret_mag.image = phim.separate_phase_components(phase_ret.image, phase_ret_reverse.image)

# Display errors in magnetic component and total phase.
phase_ret_mag.normalised_rms_error(phase_mag, display=True)
phase_ret.normalised_rms_error(phase, display=True)

# Use aliases for intensities, for brevity. Only "forward" intensity
# measurements will be plotted here.
image_under = through_focal_series.intensities[0]
image_in = through_focal_series.intensities[1]
image_over = through_focal_series.intensities[-1]

# Plot intensities and phases

image_under.plot(limits=[0, 2])
image_in.plot(limits=[0, 2])
image_over.plot(limits=[0, 2])
phase.plot(limits=[-3, 3])
phase_ret.plot(limits=[-3, 3])
phase_mag.plot(limits=[-3, 3])
phase_ret_mag.plot(limits=[-3, 3])
