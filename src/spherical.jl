# The Hilbert space spanned by the spherical harmonics obey the
# angular momentum algebra, and contractions can be expressed as
# matrix elements of spherical tensors, which can be calculated using
# Wigner 3-j symbols.

using WIGXJPF
using UnicodeFun

import Base: transpose, |, identity, float, ==, show

abstract AngularSpace <: HilbertSpace

type AngBra <: Bra
    j::Rational
    m::Rational
end
bra(::Type{AngularSpace}) = AngBra
==(b1::Bra, b2::Bra) = (b1.j == b2.j) && (b1.m == b2.m)

type AngKet <: Ket
    j::Rational
    m::Rational
end
ket(::Type{AngularSpace}) = AngKet
==(k1::Ket, k2::Ket) = (k1.j == k2.j) && (k1.m == k2.m)

transpose(b::AngBra) = AngKet(b.j, b.m)
transpose(k::AngKet) = AngBra(k.j, k.m)

type SphericalTensor <: Operator
    k::Rational
    q::Rational
end
==(C1::SphericalTensor, C2::SphericalTensor) =
    (C1.k == C2.k) && (C1.q == C2.q)

op(::Union{Type{AngularSpace},AngBra,AngKet}) = SphericalTensor
identity(::Union{Type{AngularSpace},SphericalTensor,Type{SphericalTensor}}) =
    SphericalTensor(0, 0)
space(::Union{AngBra,AngKet,SphericalTensor}) = AngularSpace

type AngScalar <: Scalar
    b::AngBra
    C::SphericalTensor
    k::AngKet
end

scalar(::Union{Type{AngularSpace},AngBra,AngKet}) = AngScalar

function float(s::AngScalar)
    full = wig3j(s.b.j, s.C.k, s.k.j,
                 -s.b.m, s.C.q, s.k.m)
    red = wig3j(s.b.j, s.C.k, s.k.j,
                0, 0, 0)
    (-1)^(2s.b.j-s.b.m)*√((2s.b.j+1)*(2s.k.j+1))*full*red
end


function show(stream::IO, s::AngScalar)
    rat(v) = den(v) == 1 || num(v) == 0 ?
        string(num(v)) : to_fraction(num(v), den(v))
    write(stream::IO, @sprintf("〈%s;%s|%s;%s|%s;%s〉 ≈ %0.5g",
                               rat(s.b.j), rat(s.b.m),
                               rat(s.C.k), rat(s.C.q),
                               rat(s.k.j), rat(s.k.m),
                               float(s)))
end

export AngularSpace
