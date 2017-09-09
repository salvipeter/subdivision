module Subdivision

# Part II. Dramatis Personae

# Section 11: An introduction to some regularly-appearing characters

"""
A subdivision scheme is defined by its mask and its arity.
The mask can be given in integers, as in `Scheme(2, [1,3,3,1])`.
A suitable divisor for normalizing is computed by `divisor`.
"""
struct Scheme
    arity :: Int
    mask  :: Array{Float64,1}
end

cubic_bspline = Scheme(2, [1,4,6,4,1])
quadratic_bspline = Scheme(2, [1,3,3,1])
ternary_quadratic_bspline = Scheme(3, [1,3,6,7,6,3,1])
ternary_neither = Scheme(3, [1,3,5,5,3,1])
fourpoint = Scheme(2, [-1,0,9,16,9,0,-1])

# A list of the above schemes, for convenience
schemes = [cubic_bspline, quadratic_bspline, ternary_quadratic_bspline, ternary_neither, fourpoint]

"Returns the number with which the mask is divided."
divisor(s) = sum(s.mask) / s.arity

"This is the mask to be used in computations."
dividedmask(s) = s.mask / divisor(s)

"The result of `stencils(s)` is a list of `s.arity` lists."
function stencils(s)
    n = length(s.mask)
    [[s.mask[j] for j in n-i:-s.arity:1] for i in 0:s.arity-1]
end

# Part III. Analyses

# Section 12: Support

"""
The call `subdivide(s,d)` subdivides `d` using the scheme `s`.
The returned value is a list of `s.arity * length(d)` elements.
When the third `dividep` argument is `false`, the mask is not normalized
(this is for its use in the `square` function below).
"""
function subdivide(s, d; dividep=true)
    n = length(s.mask)
    result = zeros(n + (length(d) - 1) * s.arity)
    for i in eachindex(d)
        result[(i-1)*s.arity+1:(i-1)*s.arity+n] += s.mask * d[i]
    end
    dividep ? result / divisor(s) : result
end

"""
Uses the scheme on itself (interpreted as 1-dimensional control points).
This results in a more dense scheme with the same support.
Iterated use of this function gives a good approximation of the basis function.
With an additional argument `n` it squares the scheme `n` times.
"""
square(s, n=1) = n < 1 ? s : square(Scheme(s.arity^2, subdivide(s, s.mask, dividep=false)), n-1)

"Returns the support of the scheme. Note that it may not be integral."
support(s) = (length(s.mask) - 1) / (s.arity - 1)

"""
`practical_supports(s, n)` approximates the support of scheme `s`,
where the basis function has considerable impact.
The result is a list of three values corresponding to the tolerances 1%, 2%, and 5%.
It uses `n` squaring operations to generate a good approximation of the basis function,
and assumes that the scheme is symmetric.
"""
function practical_supports(s, n)
    tolerances = [0.01, 0.02, 0.05]
    s2 = square(s, n)
    mask = dividedmask(s2)
    small = [findfirst(x -> abs(x) > tol, mask) for tol in tolerances] - 1
    (length(mask) - small * 2 - 1) / s2.arity
end

# Section 13: Enclosure

"""
`norm(s, n)` approximates the l_\infty norm of the scheme `s`.
It uses `n` squaring operations to generate a good approximation of the basis function.
"""
function norm(s, n)
    s2 = square(s, n)
    d = divisor(s2)
    maximum(map(row -> sum(abs(row / d)), stencils(s2)))
end

# Section 14: Continuity 1 - at Support Ends

"""
Returns the Holder-continuity at the ends of the basis function.
This is only an upper bound on the continuity of the limit curve.
"""
function endcontinuity(s)
    y0 = abs(s.mask[1]) / divisor(s)
    k = -log(s.arity, y0)
    d = ceil(k) - 1
    (d, k - d)
end

# Section 15: Continuity 2 - Eigenanalysis

"""
Returns the core parts of the subdivision matrix around mark points.
Then `eig` can be used to compute the eigenvalues and eigenvectors.
"""
matrices(s) = error("TODO")

# Section 16: Continuity 3 - Difference Schemes

"Returns a lower bound for the number of continuous derivatives."
continuity_lowerbound(s) = error("TODO")

# Section 17: Continuity 4 - Difference Eigenanalysis

const Kernel = Array{Float64,1}

kernel(s) = error("TODO")

"Returns an upper bound for the Holder-continuity of the limit curve, based on the kernel."
continuity_upperbound(s) = error("TODO")

end
