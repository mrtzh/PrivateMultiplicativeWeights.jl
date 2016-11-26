immutable GosperIterator
    n::Int64
    k::Int64
end

"""
    gosper(n, k)

Returns iterator over numbers whose binary representation correspond to the
indicator vector of a subset of {1,...,n} of size k.

Example:
   collect(gosper(3,2)) == [3, 5, 6]
"""
function gosper(n::Int, k::Int)
    @assert 0 < k <= n <= 62
    return GosperIterator(n, k)
end

typealias GosperState Tuple{Int64, Int64}

eltype(it::GosperIterator) = Int64
length(it::GosperIterator) = binomial(it.n, it.k)

function start(it::GosperIterator)
    (2^(it.k) - 1, 2^(it.n) -1)
end

function next(it::GosperIterator, state::GosperState)
    # Gosper's hack: Finds the next smallest number with exactly k bits
    # set to 1 in its binary representation.
    o = state[1]
    u = state[1] & -state[1]
    v = u + state[1]
    y = v + (div(v $ state[1], u) >> 2)
    return o, (y, state[2])
end

done(it::GosperIterator, state::GosperState) = state[1] > state[2]
