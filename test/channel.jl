using Test, OptQit, Yao,LinearAlgebra

@testset "Positive Semidefiniteness with tolerance" begin
    n = 3
    U = rand_unitary(2^n)
    eval_npsd = rand(Float64,2^n)
    eval_npsd[1] = -1e-6
    eval_npsd = normalize(eval_npsd)
    tol = -eval_npsd[1]
    A_npsd = round.(U*Diagonal(eval_npsd)*U',digits=14)
    eval_psd = rand(Float64,2^n)
    A_psd = round.(U*Diagonal(normalize(eval_psd))*U',digits=14)

    @test ispsd(A_npsd,tol=tol) ==  false
    @test ispsd(A_npsd,tol=2*tol) == true

    res, wit = ispsd_wit(A_npsd) 
    @test res == false
    @test isapprox(real(wit' * A_npsd * wit), -tol, atol=1e-14)

    res, wit = ispsd_wit(A_npsd, tol=2*tol) 
    @test res == true 

    @test ispsd(A_psd) == true
    res, wit = ispsd_wit(A_psd) 
    @test res == true 
    @test real(wit' * A_psd * wit) >= 0.0 
    @test isapprox(imag(wit' * A_psd * wit) , 0.0, atol=1e-14)

end