using Test, OptQit, Yao

@testset "Positive Semidefiniteness with tolerance" begin
    n = 3
    U = rand_unitary(2^n)
    A_npsd = rand(Float64,2^n)
    A_psd = rand(Float64,2^n)
    
    @test ispsd(U*A_npsd*U') == false
    @test ispsd(U*A_psd*U') == true

end