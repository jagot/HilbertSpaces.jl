using HilbertSpaces
using Base.Test

S = AngularSpace

i = bra(S)(1, 0)
j = ket(S)(0, 0)
C = op(S)(1, 0)

@testset "Comparisons" begin
    @test i == HilbertSpaces.AngKet(1, 0)'
    @test j == HilbertSpaces.AngBra(0, 0)'
    @test C == HilbertSpaces.SphericalTensor(1, 0)
end

@testset "Spaces" begin
    @test space(i) == AngularSpace
    @test space(j) == AngularSpace
    @test space(C) == AngularSpace
    @test op(i) == HilbertSpaces.SphericalTensor
    @test op(j) == HilbertSpaces.SphericalTensor
    @test op(HilbertSpaces.AngularSpace) == HilbertSpaces.SphericalTensor
    @test identity(op(i)) == op(S)(0, 0)
    @test bra(HilbertSpaces.AngularSpace) == HilbertSpaces.AngBra
    @test ket(HilbertSpaces.AngularSpace) == HilbertSpaces.AngKet
    @test scalar(i) == HilbertSpaces.AngScalar
    @test scalar(j) == HilbertSpaces.AngScalar
    @test scalar(HilbertSpaces.AngularSpace) == HilbertSpaces.AngScalar
end

@testset "Products" begin
    @test typeof(i|j) == scalar(S)
    @test typeof((i|C)|j) == scalar(S)
    @test typeof(i|(C|j)) == scalar(S)
    @test float(i|j) == 0
    @test float(i|i') ≈ 1
    @test float(j'|j) ≈ 1
    
    @test_throws ErrorException C|i
    @test_throws ErrorException j|C
end

@testset "Scalars" begin
    v = (i|C|j)
    @test float(v) ≈ 1/√3
end

@testset "Strings" begin
    @test string(bra(S)(0,0)) == "〈s|"
    @test string(bra(S)(1,0)) == "〈p|"
    @test string(bra(S)(1,1)) == "〈p⁺|"
    @test string(bra(S)(1,-1)) == "〈p⁻|"
    @test string(bra(S)(3,-2)) == "〈f;-2|"
    @test string(bra(S)(3//2, -1//2)) == "〈³⁄₂;⁻¹⁄₂|"
    @test string(ket(S)(1,0)) == "|p〉"

    @test string(op(S)(2,1)) == "C₁⁽²⁾"

    @test string(C|j) == "C₀⁽¹⁾|s〉"
    @test string(i|C) == "〈p|C₀⁽¹⁾"
    @test string(i|C|j) == @sprintf("〈p|C₀⁽¹⁾|s〉 ≈ %0.5g", float(i|C|j))
    @test string(i|j) == "〈p|s〉 ≈ 0"
    @test string(i|i') == "〈p|p〉 ≈ 1"
end
