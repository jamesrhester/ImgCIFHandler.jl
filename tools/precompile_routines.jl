# Instructions

"""
This file contains an exemplar run of the image checking routines.
By running the following commands at the Julia command prompt, you
can create a compiled image for faster processing, which is provided
using the -J julia option. You must provide an imgCIF file containing
external data pointers, as well as a local copy of the external data
pointed to. The script below will not work unless you change the
contents of mycif and subs accordingly!
"""

# julia> PackageCompiler.create_sysimage(sysimage_path = <output>,
# precompile_execution_file = "precompile_routines.jl")
#
# Then execute julia: julia -J<path_to_output> for no compilation time! 
#
include("image_tests.jl")

# Only useful if your CIF and data are in the ImgCIFHandler package
# test/testfiles directory

homedir = FilePaths.Path(joinpath(pkgdir(ImgCIFHandler),"test","testfiles"))

#
# Change this to a CIF file on your local machine
# e.g.
# mycif = first(Cif(p"/home/myself/stuff/c3h8.cif")).second
#
mycif = first(Cif(joinpath(homedir,p"all_scans.cif"))).second

fix_loops!(mycif)

# Change this to correspond to the external data and local copy on your
# local machine!!!
#
# e.g.
# subs = Dict("https://zenodo.org/record/1234/files/abc.tgz"=>"/home/myself/stuff/abc.tgz")
#
subs =  Dict("https://zenodo.org/record/6365376/files/cbf_m0220c.tar.bz2"=>"$(joinpath(homedir,"cbf_m0220c.tar.bz2"))")

result,img = run_img_checks(mycif,
                            images=true,
                            always=true,
                            subs = subs,
                            connected = true,
                            full = true,
                            savepng = true
                            )
