# Test interactions with simple web servers

using LiveServer

fire_up_server() = begin
    serve_dir = joinpath(@__DIR__, "testfiles")
    @async LiveServer.serve(dir = serve_dir, port=8001, verbose=true, debug=true)
end

fire_up_rsync() = begin
    test_dir = joinpath(@__DIR__, "testfiles")
    rsync_loc = joinpath(@__DIR__, "rsyncd.conf")
    rsync_config = """
use chroot = no
port = 10873
lock file = $(joinpath(@__DIR__, "rsyncd.lock"))
log file = $(joinpath(@__DIR__, "rsyncd.log"))

[files]
    path = $test_dir
"""
    f = open(rsync_loc, "w")
    write(f, rsync_config)
    close(f)
    rscmd = `rsync --daemon --config=$rsync_loc --no-detach`
    pr = run(rscmd, wait=false)
    return pr
end

prepare_cf(imgcif) = begin
    cf = first(Cif(imgcif)).second
    la = create_archives(cf)
    return cf, la[]
end

@testset "rsync + SMV" begin

    pr = fire_up_rsync()

    sleep(1) #make sure is running
    print(run(`ss -tunapls`))

    cf, la = prepare_cf(joinpath(@__DIR__, "testfiles/rsync_test.imgcif"))
    a = peek_image(la, cf)

    @test !isnothing(a)
    @test isfile(a)
    @test ImgCIFHandler.get_any_local(la, cf) == a

    img = imgload(cf, la)
    @test size(img) == (3072, 3072)
    kill(pr)
end

fire_up_server()
sleep(1)   #give server time to start

@testset "http + ZIP + CBF" begin
    cf, la = prepare_cf(joinpath(@__DIR__, "testfiles/zip_test.imgcif"))
    a = peek_image(la, cf)
    @test !isnothing(a)
    @test isfile(a)
    @test ImgCIFHandler.get_any_local(la, cf) == a

    img = imgload(cf, la)
    @test size(img) == (1475,1679)

end

@testset "http + HDF5" begin
    cf, la = prepare_cf(joinpath(@__DIR__, "testfiles/hdf5_test.imgcif"))
    a = peek_image(la, cf)
    @test !isnothing(a)
    @test isfile(a)
    @test ImgCIFHandler.get_any_local(la, cf) == a

    img = imgload(cf, la)
    @test size(img) == (101, 202)
end

@testset "http + TBZ" begin
    cf, la = prepare_cf(b4master_rem)
    a = peek_image(la, cf)
    @test !isnothing(a)
    @test isfile(a)
    @test ImgCIFHandler.get_any_local(la, cf) == a

    img = imgload(cf, "1", la)
    @test size(img) == (4148, 4362)
end

@testset "http + TIFF" begin
    cf, la = prepare_cf(joinpath(@__DIR__, "testfiles/X285_tiff_test.cif"))
    a = peek_image(la, cf)
    @test !isnothing(a)
    @test isfile(a)
    @test ImgCIFHandler.get_any_local(la, cf) == a

    img = imgload(cf, "1", la)
    @test size(img) == (2048, 2048)
end

