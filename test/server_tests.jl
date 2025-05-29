# Test interactions with simple web servers

using LiveServer

fire_up_server() = begin
    serve_dir = joinpath(@__DIR__, "testfiles/test_cbf_unzipped")
    @async LiveServer.serve(dir = serve_dir, port=8001, verbose=true)
end

prepare_cf() = begin
    cf = first(Cif(b4master_rem)).second
    la = create_archives(cf)
    return cf, la[]
end

fire_up_server()
sleep(2)   #give server time to start

@testset "http + TBZ" begin
    cf, la = prepare_cf()
    a = peek_image(la, cf)
    @test !isnothing(a)
    @test isfile(a)
    @test ImgCIFHandler.get_any_local(la, cf) == a
end
