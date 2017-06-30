using HilbertSpaces
using Base.Test

s = AngularSpace

i = bra(s)(1, 0)
j = ket(s)(0, 0)
C = op(s)(1,0)

@test space(i) == AngularSpace
@test space(j) == AngularSpace
@test space(C) == AngularSpace
@test op(i) == HilbertSpaces.SphericalTensor
@test identity(op(i)) == op(s)(0, 0)

@test float(i|j) == 0
@test float(i|i') ≈ 1
@test float(j'|j) ≈ 1

v = (i|C|j)
@test float(v) ≈ 1/√3
