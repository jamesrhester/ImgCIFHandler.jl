# imgCIF Tools

## `image_test.jl`

This program runs a series of tests on an imgCIF file. 

### Installation

1. [Install Julia](https://julialang.org/downloads) if you don't already have it.
2. Copy **all** `.jl` files and `Project.toml` from here to a convenient directory.

### Updating

Overwrite all `.jl` files and `Project.toml` from step 2 above with
the latest copies from here.  Delete the file `Manifest.toml` from the
same directory if present as otherwise the latest version of
`ImgCIFHandler.jl` may not be used.

### Usage

For help, run `julia image_test.jl --help` after installation. 

The first time `image_test.jl` is run, several minutes will be occupied with downloading and 
installing all supporting Julia packages. Subsequent runs should be much faster.

Note the `--sub <original_url> <local_file>` option (which may be repeated for multiple
urls) which links a local file with a remote URL that may be present in the imgCIF file
being checked. This
allows interactive preparation and checking of imgCIF descriptions and archive files without 
needing to download the whole archive each time the program is run.

### Examples

```
julia image_tests.jl tests/all_scans.cif
```

Run checks on `tests/all_scans.cif`, testing for the presence of any remote
archives but not downloading any images and not testing any images.

```
julia image_tests.jl -i tests/all_scans.cif
```

Run checks on `tests/all_scans.cif`, including checks on the first image
found in the remote repository specified in `tests/all_scans.cif`,
with low-resolution printout to the terminal.

```
julia image_tests.jl -i -o -p 3 tests/all_scans.cif
```

As above, but use the 3rd image (`-p 3`) found in the archive instead of
the first image and output (`-o`) a test image to file instead of terminal.

```
julia image_tests.jl --skip -o -s https://zenodo.org/record/5886687/files/cbf_b4_1.tar.bz2 /home/myself/downloads/cbf_b4_1.tar.bz2 b4_master_remote.cif
```

Do not verify imgCIF metadata (`--skip`), running tests only on the images
themselves, using 
local file `/home/myself/downloads/cbf_b4_1.tar.bz2` in place of 
`https://zenodo.org/record/5886687/files/cbf_b4_1.tar.bz2` whenever 
encountered in `b4_master_remote.cif` (`-s` option). Store test image in file
`b4_master_remote.cif.png` (`-o` option).

```
JULIA_DEBUG=Main julia image_tests.jl <filename>
```
or
```
JULIA_DEBUG=ImgCIFHandler julia image_tests.jl <filename>
```

Provide extra debugging output during the run that may or may not help.
