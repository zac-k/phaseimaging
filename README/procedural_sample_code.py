import phaseimaging as phim

specimen = phim.import_specimen('C:/Users/zac/PycharmProjects/phaseimaging/specimen')
phase_elec = phim.project_electrostatic_phase(specimen, 300e3, -17 + 1j, (100e-9, 100e-9, 100e-9))
mass_mag = 80  # emu/g
density = 5.18  # g/cm^3
magnetisation = mass_mag * density * 1000  # A/m
phase_mag = phim.project_magnetic_phase(specimen,
                                    (1,0,0),
                           magnetisation,
                                    (100e-9, 100e-9, 100e-9))

phase = phase_mag + phase_elec
phim.plot_image(phase, limits=[-3, 3])
image_under = phim.transfer_image(-8e-6, 1.96e-12, (100e-9, 100e-9), phase)
image_over = phim.transfer_image(8e-6, 1.96e-12, (100e-9, 100e-9), phase)

image_under = phim.add_noise(image_under, 1, 0.15)
image_over = phim.add_noise(image_over, 1, 0.15)

derivative = phim.intensity_derivative(image_under, image_over, 8e-6)

phase_ret = phim.retrieve_phase_tie(1.96e-12, (100e-9, 100e-9), derivative)
phim.plot_image(image_under, limits=[0, 2])
phim.plot_image(image_over, limits=[0, 2])
phim.plot_image(phase_ret, limits=[-3, 3])