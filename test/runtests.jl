# Test
using ImgCIFHandler
using CrystalInfoFramework
using DataFrames
using Test
using URIs

const b4master = joinpath(@__DIR__,"testfiles/b4_master.cif")
const b4master_rem = joinpath(@__DIR__,"testfiles/b4_master_remote.cif")
const b4master_arch = joinpath(@__DIR__,"testfiles/b4_master_archive.cif")
const multi_scan = joinpath(@__DIR__,"testfiles/all_scans.cif")

extract_files() = begin
    clean_up() # in case we failed last time
    # Uncompress archive
    archfile = joinpath(@__DIR__,"testfiles/b4_mini.tar.bz2")
    run(`bunzip2 -k $archfile`)
    # extract files into directory
    detar_dir = joinpath(@__DIR__,"testfiles/test_cbf_unzipped")
    mkpath(detar_dir)
    cd(detar_dir)
    run(`tar -xvf ../b4_mini.tar`)
end

clean_up() = begin
    rm(joinpath(@__DIR__,"testfiles/test_cbf_unzipped"),force=true,recursive=true)
    rm(joinpath(@__DIR__,"testfiles/b4_mini.tar"),force=true)
end

extract_files()

include("format_tests.jl")

# Test our own imgCIF routines

include("imgcif_tests.jl")
include("server_tests.jl")

clean_up()

