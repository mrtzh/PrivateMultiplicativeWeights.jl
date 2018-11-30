struct GosperIterator
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

const GosperState = Tuple{Int64, Int64}

eltype(it::GosperIterator) = Int64
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
    d::Int64#  dimension
    indices::Array{Int64}# indices of 1's
    init_value::Int64# value that the iterator add the subnumbers to it
end

"""
For some a in {0,1}^d, that have 1's in the places "indices", iterate over all substes 
of indices in "indices" and return the number of it's binary representation
"""
function BinaryItr(d::Int64, indices::Array{Int64}, init_value::Int64 = 0)
    @assert 0 < d <= 62

    BinaryIndexIterator(d, indices, init_value)
end

"""
Now the indices of the 1's are the number itself with those indices
"""
function BinaryItr(d::Int64, indices::Int64, init_value::Int64 = 0)
    @assert 0 < d <= 62

    _, arr_idx = norm_1_with_indices(indices)

    BinaryIndexIterator(d, Array{Int64}(arr_idx), init_value)
end

const BinaryIndexState = Int64

eltype(it::BinaryIndexIterator) = Int64
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
function norm_1_with_indices(alpha::Int64)

    @assert alpha >= 0

    count = 0
    index = 1
    idx = []
    while alpha > 0
        if alpha % 2 == 1
            count += 1
            push!(idx, index)
        end
        alpha >>= 1
        index += 1
    end
    count, Array{Int64}(idx)
end

"""
Calculate the norm1 of binary number (the number of 1's in the number).
return: -norm1(alpha)
        -array with the indices of the 1's
"""
function norm_1(alpha::Int64)

    @assert alpha >= 0

    count = 0
    while alpha > 0
        if alpha % 2 == 1
            count += 1
        end
        alpha >>= 1
    end
    count
end