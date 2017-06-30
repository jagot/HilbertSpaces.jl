module HilbertSpaces

import Base: |, identity

abstract HilbertSpace

abstract Bra
abstract Ket

abstract Operator

abstract Scalar

|(b::Bra, k::Ket) = scalar(k)(b, identity(op(k)), k)

|(b::Bra, O::Operator) = (b,O)
|(O::Operator, k::Ket) = (O,k)

|(b::Bra, Ok::Tuple{Operator, Ket}) = scalar(b)(b, Ok[1], Ok[2])
|(bO::Tuple{Bra, Operator}, k::Ket) = scalar(k)(bO[1], bO[2], k)

include("spherical.jl")

export ket, bra, op, space

end # module
