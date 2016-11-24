typealias GosperState Tuple{Int64, Int64}

immutable GosperIterator
    n::Int64
    k::Int64
end

eltype(it::GosperIterator) = Int64
length(it::GosperIterator) = binomial(it.n, it.k)

function gosper(n::Int, k::Int)
    @assert 0 < k <= n <= 62
    return GosperIterator(n, k)
end

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

# Iterate over all subsets of a collection with a given size
typealias BinomialState Tuple{Array{Int64, 1}, Bool}

immutable Binomial{T}
    xs::Array{T, 1}
    n::Int64
    k::Int64
end

eltype(it::Binomial) = Array{eltype(it.xs), 1}
length(it::Binomial) = binomial(it.n, it.k)
subsets(xs, k) = Binomial(xs, length(xs), k)
start(it::Binomial) = (Int64[1:it.k], false)

function next(it::Binomial, state::BinomialState)
    idx = state[1]
    set = it.xs[idx]
    i = it.k
    while(i>0)
        if idx[i] < it.n - it.k + i
            idx[i] += 1
            idx[i+1:it.k] = [idx[i]+1:idx[i]+it.k-i]
            break
        else
            i -= 1
        end
    end
    if i==0
        return set, (idx, true)
    else
        return set, (idx, false)
    end
end

done(it::Binomial, state::BinomialState) = state[2]
