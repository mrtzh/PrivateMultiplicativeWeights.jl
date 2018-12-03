using DataFrames

"""
translate the element (row in the dataset) to index in the 
histogram of the domain, according to the columns dimension sizes.
params:
    element:    iterable that contain the values of the its columns.
    domains:    the domain sizes of the columns.
"""
function get_flatten_index(element, domains)
    counter = 0
    pow = 1
    for (i,s) in zip(element,domains)
        counter += i*pow
        pow *= s
    end
    counter
end

"""
translate the index in the histogram of the domain to the corresponding 
element in the dataset (row in the dataset), according to the columns dimension sizes.
params:
    index:    index in the histogram (integer).
    domains:    the domain sizes of the columns.
"""
function reverse_flatten_index(index, dim)

    new_index = zeros(Int64, length(dim))

    for (i,s) in enumerate(dim)
        new_index[i] += index % s
        index -= index % s
        index /= s
    end

    new_index

end

"""
get a list of values, and return how many times the each value is found.
return: dictionary{value : counter}
"""
function countMembers(itr)
    d = Dict{eltype(itr), Int64}()
    for val in itr
        if isa(val, Number) && isnan(val)
            continue
        end
        d[val] = get!(d, val, 0) + 1
    end
    return d
end

"""
create a mapping for the elements of a domain, into integers
    params:
        domain_elements:    the elements of the domain
    return:
        A tuple of two mappings: (index => elemnt, elemnt => index)
"""
function getMapping(domain_elements::Array)

    direct_mapping = Dict{Int64, Any}()
    reverse_mapping2 = Dict{Any, Int64}()

    for (i,element) in enumerate(domain_elements)
        direct_mapping[i-1] = element
        reverse_mapping2[element] = i-1
    end

    direct_mapping, reverse_mapping2

end

"""
Get a dataset and make it into a histogram of its domain.
    params:
        data:               A dataframe containing the dataset
        domain_elements:    list of the domain elements of every column (a list of lists)
    return:
        histogram:  the histogram of the dataset over the domain.
        mapping:    a list of tuples of mappings (index => elemnt, elemnt => index) for every column domain.
"""
function createHistogram(data::DataFrames.DataFrame, domains_elements::Array)

    n, d = size(data)

    domain_sizes = [length(domain)  for domain in domains_elements]
    histogram_size = prod(domain_sizes)

    histogram = zeros(Float64, histogram_size)

    mappings = []

    for domain_elements in domains_elements
        push!(mappings, getMapping(domain_elements))
    end
    
    for i = 1:n
        x = [data[i,j] for j=1:d]
        index = [mappings[j][2][x[j]] for j=1:d]
        histogram_index = get_flatten_index(index, domain_sizes) + 1
        histogram[histogram_index] += 1
    end

    histogram, mappings

end

""""
create histogram from a matrix of data.
"""
function createHistogram(data::AbstractArray, domains_elements::Array)

    createHistogram(DataFrame(samples), domains_elements)

end

"""
Recreate a database from histogram according to the mapping of the original elements of the dataset domain.
(the mapping from the creating of the histogram. see "createHistogram" function)
params:
    histogram:  a histogram of a dataset.
    mapping:    a list of tuples of mappings (index => elemnt, elemnt => index) for every column domain of the dataset.
    domain_elements:    list of the domain elements of every column (a list of lists)

return:
    df: A dataframe containing the recreated dataset
"""
function recreateDataset(histogram::Array{Float64, 1}, mappings::Array, domains_elements::Array)

    n = sum(histogram)
    d = length(domains_elements)

    domain_sizes = [length(domain)  for domain in domains_elements]

    df = DataFrame()
    for i=1:d
        df[i] = []
    end

    for index = 1:length(histogram)
        domain_element_encoded = reverse_flatten_index(index-1, domain_sizes)
        domain_element = [mappings[j][1][domain_element_encoded[j]] for j=1:d]
        for i=1:histogram[index]
            push!(df, copy(domain_element))
        end
    end

    df

end

function readContingencyTable(data::DataFrames.DataFrame, dimension::Int64)

    @assert length(size(data)) == 2

    n,temp = size(data)

    @assert temp == 2
    @assert n > 0

    hist = zeros(Float64, 2^dimension)

    for i in 1:n
        dig = digits(data[[:1]][1][i], base = 10)
        index = sum([dig[k]*2^(k-1) for k=1:length(dig)])
        hist[index] += data[[:2]][1][i]
    end

    hist
end