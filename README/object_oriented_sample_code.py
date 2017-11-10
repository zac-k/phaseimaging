import phaseimaging as phim

mass_mag = 80  # emu/g
density = 5.18  # g/cm^3
magnetisation = mass_mag * density * 1000  # A/m

specimen = phim.Specimen(width = (100e-9, 100e-9, 100e-9),
                         mean_inner_potential=-17+1j,
                         magnetisation=magnetisation,
                         mhat=(1, 0, 0),
                         specimen_file='C:/Users/zac/PycharmProjects/phaseimaging/specimen')
phase = phim.Phase(resolution=specimen.resolution[0:2], width=specimen.width[0:2])
beam = phim.Beam(phim.accel_volt_to_lambda(300e3))
phase.project_electrostatic(specimen, beam)
phase.project_magnetic(specimen, beam)

through_focal_series = phim.ThroughFocalSeries(phase.resolution, phase.width, [-8e-6, 0, 8e-6])
through_focal_series.transfer_images(phase, beam)
through_focal_series.add_noise(0.05)
through_focal_series.compute_derivative()

image_under = through_focal_series.intensities[0]
image_over = through_focal_series.intensities[-1]

phase_ret = phim.Phase(resolution=specimen.resolution[0:2], width=specimen.width[0:2])
phase_ret.retrieve_phase_tie(through_focal_series, beam)

phase.plot(limits=[-3, 3])
image_under.plot(limits=[0, 2])
image_over.plot(limits=[0, 2])
phase_ret.plot(limits=[-3, 3])
