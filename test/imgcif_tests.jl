# Test our imgCIF handling routines

@testset "Test gonio axis handling" begin
    cc = first(Cif(Path(multi_scan),native=true)).second
    gonio_axes = get_gonio_axes(cc)
    @test gonio_axes[1] == ["OMEGA","KAPPA","PHI"]
end

@testset "Test axis calculations" begin
    cc = first(Cif(Path(multi_scan),native=true)).second
    @test get_dependency_chain(cc,"PHI")==["PHI","KAPPA","OMEGA"]
end

