@testset "Test CBF file loading" begin
    q = joinpath(@__DIR__,"testfiles/test_cbf_unzipped/s01f0002.cbf")
    x = imgload(q,Val(:CBF))
    @test size(x) == (4148,4362)
end

@testset "Test HDF5 file loading" begin
    q = joinpath(@__DIR__,"testfiles/_b4_three_frames.h5")
    x = imgload(q,Val(:HDF5),path="/data",frame=1)
    @test size(x) == (101,202)
end

@testset "Test ADSC file loading" begin
    q = joinpath(@__DIR__,"testfiles/tartaric_2_003.img")
    x = imgload(q,Val(:SMV))
    @test size(x) == (2048,2048)
end

@testset "Test TIFF file loading" begin
    q = joinpath(@__DIR__, "testfiles/test_XRD285_3frames.tiff")
    x = imgload(q, Val(:TIFF))
    @test size(x) == (2048, 2048)
end

@testset "Test variants of imgload" begin
    x = imgload(b4master)
    @test size(x) == (4148,4362)
end

@testset "Test extraction from archive" begin
    
    loc = unescapeuri(joinpath(@__DIR__,"testfiles/b4_mini.tar"))
    dl_info = DataFrame(:uri=>URI(scheme="file",path=loc),
                        :format=>"CBF",
                        :archive_format=>"TAR",
                        :archive_path=>"s01f0003.cbf")
    a = create_archives(dl_info.uri[], arch_type = "TAR")
    x = imgload(first(a), dl_info)
    @test size(x) == (4148,4362)
    
    loc = unescapeuri(joinpath(@__DIR__,"testfiles/b4_mini.tar.bz2"))
    dl_info.archive_format=["TBZ"]
    dl_info.uri=[URI(scheme="file",path=loc)]
    a = create_archives(dl_info.uri[], arch_type = "TBZ")
    x = imgload(first(a), dl_info)
    @test size(x) == (4148,4362)
end

@testset "Test extraction of multiple images from archive" begin
    cf_string = read(b4master_arch, String)
    # put in absolute path
    cf_string = replace(cf_string, "b4_mini.tar.bz2" => "file:" * joinpath(@__DIR__, "testfiles", "b4_mini.tar.bz2"))
    cf = first(cif_from_string(cf_string)).second
    a = first(create_archives(cf))
    imgsum = imgload(cf,["1","2","3"],a)
    @test size(imgsum) == (4148,4362)
    # check that we have a sum
    imgs = imgload.(Ref(cf),["1","2","3"], Ref(a))
    @test imgs[1]+imgs[2]+imgs[3] == imgsum
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
