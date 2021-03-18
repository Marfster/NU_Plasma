# NU Data Reduction

This program was developed during my PhD. It performs the data reduction and interference correction of data produced by a Multi-Collector Inductively-Coupled Mass Spectrometer (NU Plasma II).

## Code for data reduction: `nu_data_reduction.py`
Classes:
* `NU_data_read` - for Nu Plasma II: all methods to read the raw data, perform baseline correction and store the data into a dictionary
* `Neptune_data_read` - for Thermo Neptune Plus: all methods to read the raw data, perform baseline correction and store the data into a dictionary
* `normalisation` - Methods used for mass fractionation and interference correction for a single cycle of a sample measurement
* `evaluation` - Methods for background correction, mass fractionation correction, outlier rejection for a whole sample measurement

## Code for Isotope configuration dictionaries: `iso_properties.py`
Classes:



### TODO:
  - finish Readme file
  - move from python dictionaries and pandas Dataframes to Numpy arrays to improve speed
  - move data evaluation of a measurement day from Jupyter notebooks to python file
  - safe data into PostgresSQL databases instead of HDF5 files
