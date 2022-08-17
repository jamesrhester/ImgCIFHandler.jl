# Test our imgCIF handling routines

@testset "Test beam centre calculation" begin
    (c1,c2),(i1,i2) = get_beam_centre(b4master)
    @test isapprox(c1, 172.4595, atol= 0.0001)
end

@testset "Test scan, frame no extraction" begin
    bb = first(Cif(Path(b4master))).second
    s,f = ImgCIFHandler.scan_frame_from_img_name("test_cbf_unzipped/s01f0003.cbf",nothing,bb)
    println("$s,$f")
    @test s == "SCAN1" && f == 3
    bb = first(Cif(Path(b4master_rem))).second
    uri = "https://zenodo.org/record/5886687/files/cbf_b4_1.tar.bz2"
    s,f = ImgCIFHandler.scan_frame_from_img_name(uri,"s01f0014.cbf",bb)
    @test s == "SCAN1" && f == 14
    bb = first(Cif(Path(multi_scan))).second
    uri = "https://zenodo.org/record/6365376/files/cbf_m0220c.tar.bz2"
    s,f = ImgCIFHandler.scan_frame_from_img_name(uri,"m0220c_02_0171.cbf",bb)
    @test s == "SCAN02" && f == 171
    ext_info = ImgCIFHandler.external_specs_from_bin_ids(["1410","1411","1412"],bb)
    @test ext_info.uri[1] == "https://zenodo.org/record/6365376/files/cbf_m0220c.tar.bz2"
    @test size(ext_info,1) == 3
    @test ext_info.archive_path[3] == "m0220c_02_0212.cbf"
end

@testset "Test gonio axis handling" begin
    cc = first(Cif(Path(multi_scan),native=true)).second
    gonio_axes = get_gonio_axes(cc)
    @test gonio_axes[1] == ["OMEGA","KAPPA","PHI"]
end

@testset "Test axis calculations" begin
    cc = first(Cif(Path(multi_scan),native=true)).second
    @test get_dependency_chain(cc,"PHI")==["PHI","KAPPA","OMEGA"]
end

@testset "Test pixel position calculations" begin

    for (testfile,sname) in zip((b4master,multi_scan),("1","02"))

        cbf_pc = get_pixel_coordinates(testfile,5,9,"SCAN$sname",6)

        tfc = first(Cif(Path(testfile))).second
        
        ich_pc = get_pixel_coordinates(tfc,5,9,"SCAN$sname",6)

        @test isapprox(ich_pc, cbf_pc, atol=1e-6)
    end
end


@testset "Test beam centre calculations" begin

    # compare to cbf calculations

    for (testfile,sname) in zip((b4master,multi_scan),("1","02"))
        cbf_centre = get_beam_centre(testfile)
        cbf_centre_scan = get_beam_centre(testfile,"SCAN$sname",2)
        
        tfc = first(Cif(Path(testfile))).second
    
        ich_centre = get_beam_centre(tfc)
        ich_centre_scan = get_beam_centre(tfc,"SCAN$sname",2)

        @test isapprox(ich_centre[1], cbf_centre[1], atol=1e-6)
        @test isapprox(ich_centre_scan[1], cbf_centre_scan[1], atol=1e-6)
        @test isapprox(ich_centre[2], cbf_centre[2], atol=1e-6)
        @test isapprox(ich_centre_scan[2], cbf_centre_scan[2], atol=1e-6)
    end
    
end


