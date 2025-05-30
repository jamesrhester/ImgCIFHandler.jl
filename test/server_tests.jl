# Test interactions with simple web servers

using LiveServer

fire_up_server() = begin
    serve_dir = joinpath(@__DIR__, "testfiles")
    @async LiveServer.serve(dir = serve_dir, port=8001, verbose=true, debug=true)
end

prepare_cf(imgcif) = begin
    cf = first(Cif(imgcif)).second
    la = create_archives(cf)
    return cf, la[]
end

fire_up_server()
sleep(1)   #give server time to start

@testset "http + TBZ" begin
    cf, la = prepare_cf(b4master_rem)
    a = peek_image(la, cf)
    @test !isnothing(a)
    @test isfile(a)
    @test ImgCIFHandler.get_any_local(la, cf) == a
end

@testset "http + TIFF" begin
    cf, la = prepare_cf(joinpath(@__DIR__, "testfiles/X285_tiff_test.cif"))
    a = peek_image(la, cf)
    @test !isnothing(a)
    @test isfile(a)
    @test ImgCIFHandler.get_any_local(la, cf) == a
end
