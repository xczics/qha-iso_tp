from qha.api import qha

# Create a qha object
Job={
    'calculation': 'iso',
    'input_from_file': True,
    'T_MIN': 5,
    'NT': 405,
    'DT': 5,
    'NTV': 220,
    'DELTA_P': 1.0,
    'input': 'input01_Liso',
    'input_h': 'input01_Hiso',
    'P_MIN': -2.1,
    'high_verbosity': True,
    'thermodynamic_properties':['F','U','H','V','iso'],
    'energy_unit': 'ry',
    'output_directory': './results/',
    'freq2THz': 0.029979245368431,
    'an': 2,
}
QHAJob=qha(**Job)
print(QHAJob.check())
QHAJob.run()
print(QHAJob.thermodict()['iso'].shape)
ISO=QHAJob.propertiesTP('iso')
print(ISO(997.5,10.0))
V=QHAJob.propertiesTP('V')
print(V(1000,10.0))
print(ISO(997.5,10.0))