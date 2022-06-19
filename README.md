![Testing](https://github.com/jamesrhester/ImgCIFHandler.jl/workflows/Tests/badge.svg)

# ImgCIFHandler

This Julia package provides functions for working with ImgCIF files. imgCIF is a raw data standard developed
by the International Union of Crystallography, and recently enhanced with data
pointers so that the raw data can be preserved outside of the imgCIF file 
itself.

# Installation

## Prerequisites

For some functions, the programs `tar`, `curl` and `awk` must be available in the environment, 
as well as `gzip` and `bunzip2`. This requirement may be removed in the future.

## Installation

To install this package, at the Julia package prompt (reached by typing `]` at the Julia prompt) type `add ImgCIFHandler`.
See also additional instructions in the `tools` directory for installing and updating the imgCIF checking tool.

# Use

After installation, typing 'using ImgCIFHandler' will import the package. The following functions will be in scope:

1. `imgload`. This function can be used to
access raw data stored using pointers within imgCIF files. The example file [`b4_master.cif`](test/testfiles/b4_master.cif) contains such pointers.
For example, to read the raw data array labelled `ext1` in `b4_master.cif` the following call is used:

```julia
x = imgload("b4_master.cif","ext1")
```

Note that this command will only work if the images referred to in `b4_master.cif` actually exist. To create those files locally, unpack the 
`b4_mini.tar.bz2` contents into location `test/testfiles/test_cbf_unzipped`.

2. `get_beam_centre(filename)`. This will calculate the beam centre in pixel and detector coordinates with all axes at zero. 
`get_beam_centre(filename,scan_name,frameno)` will calculate the beam centre for the given frame number of the scan.

# Supported external formats

This release supports `HDF5`, `CBF`, and `SMV` (ADSC) formats, as well as Tar archives (optionally compressed using Bzip2 or Gzip) or 
Zip archives containing files in those formats. In all cases the format stated by the imgCIF file is used.

# Contributions

Contributions are welcome. You may raise an issue using the tab above, or clone the repository and create a pull request with your
changes.
