struct GosperIterator
    n::Int
    k::Int
end

"""
    gosper(n, k)

Returns iterator over numbers whose binary representation correspond to the
indicator vector of a subset of {1,...,n} of size k.

Example:
   collect(gosper(3,2)) == [3, 5, 6]
"""
function gosper(n::Integer, k::Integer)
    0 < k <= n <= 62 || throw(ArgumentError("n and k must satisfy 0 < k <= n <= 62"))
    return GosperIterator(Int(n), Int(k))
end

const GosperState = Tuple{Int, Int}

eltype(::Type{GosperIterator}) = Int
length(it::GosperIterator) = binomial(it.n, it.k)

function start(it::GosperIterator)
    (2^(it.k) - 1, 2^(it.n) -1)
end

function iterate(it::GosperIterator, state::GosperState = start(it))

    if state[1] > state[2]
        return nothing
    end

        # Gosper's hack: Finds the next smallest number with exactly k bits
    # set to 1 in its binary representation.
    o = state[1]
    u = state[1] & -state[1]
    v = u + state[1]
    y = v + (div(xor(v, state[1]), u) >> 2)
    return o, (y, state[2])
end


struct BinaryIndexIterator
    d::Int # dimension
    indices::Vector{Int} # indices of 1's
    init_value::Int # value that the iterator adds the subnumbers to
end

"""
For some `a` in `{0,1}^d` with ones at `indices`, iterate over every subset of
`indices` and return its binary representation.
"""
function BinaryItr(d::Integer, indices::AbstractVector{<:Integer}, init_value::Integer=0)
    0 < d <= 62 || throw(ArgumentError("dimension must be between 1 and 62"))
    all(i -> 1 <= i <= d, indices) || throw(ArgumentError("indices must be within the dimension"))
    allunique(indices) || throw(ArgumentError("indices must be unique"))
    0 <= init_value < (1 << d) || throw(ArgumentError("init_value must be within the dimension"))
    all(i -> init_value & (1 << (i - 1)) == 0, indices) ||
        throw(ArgumentError("init_value must not overlap with indices"))

    BinaryIndexIterator(Int(d), Int.(indices), Int(init_value))
end

"""
Now the indices of the 1's are the number itself with those indices
"""
function BinaryItr(d::Integer, indices::Integer, init_value::Integer=0)
    0 < d <= 62 || throw(ArgumentError("dimension must be between 1 and 62"))
    0 <= indices < (1 << d) || throw(ArgumentError("bit mask is outside the dimension"))

    _, arr_idx = norm_1_with_indices(indices)

    BinaryItr(d, arr_idx, init_value)
end

const BinaryIndexState = Int

eltype(::Type{BinaryIndexIterator}) = Int
length(it::BinaryIndexIterator) = 2^(length(it.indices))

function start(it::BinaryIndexIterator)
   0
end

function iterate(it::BinaryIndexIterator, state::BinaryIndexState = start(it))

    if state == length(it)
        return nothing
    end

    sub_binary_number = it.init_value

    temp_state = state

    for i=1:length(it.indices)
        if temp_state == 0
            break
        end
        if temp_state%2 == 1
            sub_binary_number += 1<<(it.indices[i]-1)
        end
        temp_state >>= 1
    end

    state += 1

    return sub_binary_number, state

end

"""
Calculate the norm1 of binary number (the number of 1's in the number).
return: -norm1(alpha)
        -array with the indices of the 1's
"""
function norm_1_with_indices(alpha::Integer)
    alpha >= 0 || throw(ArgumentError("alpha must be nonnegative"))

    count = 0
    index = 1
    idx = Int[]
    while alpha > 0
        if alpha % 2 == 1
            count += 1
            push!(idx, index)
        end
        alpha >>= 1
        index += 1
    end
    count, idx
end

"""
Calculate the norm1 of binary number (the number of 1's in the number).
return: -norm1(alpha)
        -array with the indices of the 1's
"""
function norm_1(alpha::Integer)
    alpha >= 0 || throw(ArgumentError("alpha must be nonnegative"))

    count = 0
    while alpha > 0
        if alpha % 2 == 1
            count += 1
        end
        alpha >>= 1
    end
    count
end
