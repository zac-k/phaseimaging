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
phase = phim.Phase(resolution=specimen.resolution[0:2], width=specimen.width[0:2])
beam = phim.Beam(phim.accel_volt_to_lambda(300e3))

# Project the electrostatic and magnetic phases
phase.project_magnetic(specimen, beam)
phase_mag = copy.deepcopy(phase)
phase.project_electrostatic(specimen, beam)

# Set noise level and defoci
noise_level = 0.01
defoci = [-8e-6, 0, 8e-6]

# Generate through-focal series of intensities and compute derivative
through_focal_series = phim.ThroughFocalSeries(phase.resolution, phase.width, defoci=defoci)
through_focal_series.transfer_images(phase, beam)
through_focal_series.add_noise(noise_level)
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
phase_ret.remove_offset()
phase_ret.normalised_rms_error(phase, display=True)


specimen.rotate(angle=180, axis=0)

phase_reverse = phim.Phase(resolution=specimen.resolution[0:2], width=specimen.width[0:2])
phase_reverse.project_electrostatic(specimen, beam)
phase_reverse.project_magnetic(specimen, beam)
through_focal_series_reverse = phim.ThroughFocalSeries(phase.resolution, phase.width, defoci=defoci)
through_focal_series_reverse.transfer_images(phase_reverse, beam)
through_focal_series_reverse.add_noise(noise_level)
through_focal_series_reverse.compute_derivative()

phase_ret_reverse = phim.Phase(resolution=specimen.resolution[0:2], width=specimen.width[0:2])
phase_ret_reverse.retrieve_phase_tie(through_focal_series_reverse, beam)
phase_ret_reverse.remove_offset()


phase_ret_reverse.flip(axis=1)
phase_ret_elec = phim.Phase(resolution=specimen.resolution[0:2], width=specimen.width[0:2])
phase_ret_mag = phim.Phase(resolution=specimen.resolution[0:2], width=specimen.width[0:2])
phase_ret_elec.image, phase_ret_mag.image = phim.separate_phase_components(phase_ret.image, phase_ret_reverse.image)

phase_ret_mag.normalised_rms_error(phase_mag, display=True)
phase_ret.normalised_rms_error(phase, display=True)

# Plot intensities and phases
phase.plot(limits=[-3, 3])
phase_reverse.plot(limits=[-3, 3])
image_under.plot(limits=[0, 2])
image_in.plot(limits=[0, 2])
image_over.plot(limits=[0, 2])
phase_ret.plot(limits=[-3, 3])
phase_ret_reverse.plot(limits=[-3, 3])
phase_mag.plot(limits=[-3, 3])
phase_ret_mag.plot(limits=[-3, 3])
