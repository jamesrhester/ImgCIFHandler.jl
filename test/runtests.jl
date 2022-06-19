# Test
using ImgCIFHandler
using Test
using URIs

const b4master = joinpath(@__DIR__,"testfiles/b4_master.cif")

extract_files() = begin
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
    rm(joinpath(@__DIR__,"testfiles/test_cbf_unzipped"),recursive=true)
    rm(joinpath(@__DIR__,"testfiles/b4_mini.tar"))
end

@testset "Test HDF5 file loading" begin
    q = joinpath(@__DIR__,"testfiles/simple3D.h5")
    x = imgload(q,Val(:HDF),path="/entry/data/test",frame=1)
    @test size(x) == (4,3)
end

@testset "Test ADSC file loading" begin
    q = joinpath(@__DIR__,"testfiles/tartaric_2_003.img")
    x = imgload(q,Val(:SMV))
    @test size(x) == (2048,2048)
end

@testset "Test variants of imgload" begin
    extract_files()
    x = imgload(b4master)
    @test size(x) == (4148,4362)
end

@testset "Test CBF file loading" begin
    q = joinpath(@__DIR__,"testfiles/test_cbf_unzipped/s01f0002.cbf")
    x = imgload(q,Val(:CBF))
    @test size(x) == (4148,4362)
end

@testset "Test extraction from archive" begin
    loc = unescapeuri(joinpath(@__DIR__,"testfiles/b4_mini.tar"))
    x = imgload(URI(scheme="file",path=loc),Val(:CBF),arch_type="TAR",arch_path="s01f0003.cbf")
    @test size(x) == (4148,4362)
    loc = unescapeuri(joinpath(@__DIR__,"testfiles/b4_mini.tar.bz2"))
    x = imgload(URI(scheme="file",path=loc),Val(:CBF),arch_type="TBZ",arch_path="s01f0003.cbf")
    @test size(x) == (4148,4362)
end

@testset "Test axis getting and setting" begin
    a,t,p = get_detector_axis_settings(b4master)
    @test issetequal(a,["trans","two_theta"])
    a,t,p = get_detector_axis_settings(b4master,"SCAN1",3)
    trans = indexin(["trans"],a)[]
    @test p[trans] == 287.22
    handle = ImgCIFHandler.cbf_read_file(b4master)
    twotheta = indexin(["two_theta"],a)[]
    p[twotheta] = 21.2
    ImgCIFHandler.cbf_set_axis_positions(handle,a,t,p)

    # Read it back in

    ImgCIFHandler.cbf_find_category(handle,"diffrn_scan_frame_axis")
    ImgCIFHandler.cbf_find_column(handle,"axis_id")
    ImgCIFHandler.cbf_find_row(handle,"two_theta")
    ImgCIFHandler.cbf_find_column(handle,"angle")
    tt = ImgCIFHandler.cbf_get_string_value(handle)
    @test tt == "$(p[twotheta])"
end

@testset "Test beam centre calculation" begin
    c1,c2,i1,i2 = get_beam_centre(b4master)
    @test isapprox(c1, 172.4595, atol= 0.0001)
end


clean_up()
