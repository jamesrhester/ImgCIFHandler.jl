# imgCIF Tools

## `image_test.jl`

This program runs a series of tests on an imgCIF file. 

### Installation

1. [Install Julia](https://julialang.org/downloads) if you don't already have it.
2. Copy **all** `.jl` files and `Project.toml` from here to a convenient directory.

### Usage

For help, run `julia image_test.jl --help` after installation. 

The first time `image_test.jl` is run, several minutes will be occupied with downloading and 
installing all supporting Julia packages. Subsequent runs should be much faster.

Note the `--sub <original_url> <local_file>` option (which may be repeated for multiple
urls) which links a local file with a remote URL that may be present in the imgCIF file
being checked. This
allows interactive preparation and checking of imgCIF descriptions and archive files without 
needing to download the whole archive each time the program is run.

### Updating

Overwrite all `.jl` files and `Project.toml` from step 2 above with
the latest copies from here.  Delete the file `Manifest.toml` from the
same directory if present as otherwise the latest version of
`ImgCIFHandler.jl` may not be used.
