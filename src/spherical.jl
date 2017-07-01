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
    AngBra(j,m) = abs(m) > j ? error("Invalid bra, |$(m)| > $(j)") : new(j, m)
end
bra(::Type{AngularSpace}) = AngBra
==(b1::Bra, b2::Bra) = (b1.j == b2.j) && (b1.m == b2.m)

type AngKet <: Ket
    j::Rational
    m::Rational
    AngKet(j,m) = abs(m) > j ? error("Invalid ket, |$(m)| > $(j)") : new(j, m)
end
ket(::Type{AngularSpace}) = AngKet
==(k1::Ket, k2::Ket) = (k1.j == k2.j) && (k1.m == k2.m)

transpose(b::AngBra) = AngKet(b.j, b.m)
transpose(k::AngKet) = AngBra(k.j, k.m)

type SphericalTensor <: Operator
    k::Integer # Rank
    q::Integer
    SphericalTensor(k, q) = abs(q) > k ? error("Invalid tensor, |$(q)| > $(k) ") : new(k, q)
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

const ells = "spdfghiklmnopqrtuvwxyz"

function partial_wave_str(ell, m)
    rat(v) = den(v) == 1 || num(v) == 0 ?
        string(num(v)) : to_fraction(num(v), den(v))

    if ell == 0
        "s"
    elseif ell == 1
        if m == 0
            "p"
        elseif m == 1
            "p⁺"
        else
            "p⁻"
        end
    elseif den(ell) == 1 && den(m) == 1
        @sprintf("%s;%i", ells[num(ell)+1], m)
    else
        @sprintf("%s;%s", rat(ell), rat(m))
    end
end

show(stream::IO, b::AngBra) =
    write(stream::IO, @sprintf("〈%s|", partial_wave_str(b.j, b.m)))
show(stream::IO, k::AngKet) =
    write(stream::IO, @sprintf("|%s〉", partial_wave_str(k.j, k.m)))
show(stream::IO, C::SphericalTensor) =
    write(stream::IO, @sprintf("C%s%s",
                               to_subscript(C.q),
                               to_superscript("($(C.k))")))

function show(stream::IO, Ck::Tuple{SphericalTensor, AngKet})
    show(stream, Ck[1])
    show(stream, Ck[2])
end

function show(stream::IO, bC::Tuple{AngBra, SphericalTensor})
    show(stream, bC[1])
    show(stream, bC[2])
end

function show(stream::IO, s::AngScalar)
    show(stream, s.b)
    if s.C.k == 0 && s.C.q == 0
        write(stream, partial_wave_str(s.k.j, s.k.m))
        write(stream, "〉")
    else
        show(stream, s.C)
        show(stream, s.k)
    end
    write(stream::IO, @sprintf(" ≈ %0.5g", float(s)))
end

export AngularSpace
