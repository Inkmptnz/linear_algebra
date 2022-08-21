module GenerateMatrices

using LinearAlgebra 
using AbstractAlgebra
using Polynomials
using DataStructures

# public functions

# gives Matrices with whole, non imaginary EVs

matrix_whole_ev(dim, possible_values, is_whole_ev)

matrix_good_ev!(dim, start, stop) = matrix_good_ev!(dim, start:stop)

function matrix_good_ev!(dim, possible_values, is_good_ev::Function)
    mat = rand(dim, dim)
    matrix_good_ev!(mat, possible_values)
end

matrix_good_ev!(mat::AbstractArray, start, stop) = matrix_good_ev!(mat, start:stop)

function matrix_good_ev!(mat::AbstractArray, possible_values)
    evs = zeros(size(mat, 1))
    evcache = similar(mat)
    evs[1] = 0.1
    while !is_whole_ev(evs)
        evs = random_try!(mat, evcache, possible_values)
    end
    convert(Matrix{Int64}, mat)
end

is_whole_ev(::Vector{ComplexF64}) = false

function is_whole_ev(evs)
    for ev in evs
        if ev ≠ floor(ev)
            return false
        end
    end
    return true
end

# gives matrices with only imaginary EVs, but whole
matrix_whole_imag_ev!(dim,start,stop) = matrix_whole_imag_ev!(dim, start:stop)

function matrix_whole_imag_ev!(dim, possible_values)
    mat = rand(dim, dim)
    matrix_whole_imag_ev!(mat, possible_values)
end

matrix_whole_imag_ev!(mat::AbstractArray,start,stop) = matrix_whole_imag_ev!(mat, start:stop)

function matrix_whole_imag_ev!(mat::AbstractArray, possible_values)
    evs = zeros(size(mat, 1))
    evcache = similar(mat)
    evs[1] = 0.1 # ensures that while runs at least ones
    while !is_whole_imag_ev(evs)
        evs = random_try!(mat, evcache, possible_values)
    end
    convert(Matrix{Int64}, mat)
end

is_whole_imag_ev(::Vector{Float64}) = false

function is_whole_imag_ev(evs)
    for ev in evs
        if imag(ev) ≠ floor(imag(ev))
            return false
        end
    end
    true
end

matrix_multiple_whole_ev!(dim,start,stop,equal_count) = matrix_whole_imag_ev!(dim, start:stop, equal_count)

function matrix_multiple_whole_ev!(dim, possible_values, equal_count)
    mat = rand(dim, dim)
    matrix_whole_imag_ev!(mat, possible_values, equal_count)
end

matrix_multiple_whole_ev!(mat::AbstractArray,start,stop,equal_count) = matrix_whole_imag_ev!(mat, start:stop, equal_count)

function matrix_multiple_whole_ev!(mat::AbstractArray, possible_values, equal_count)
    evs = zeros(size(mat, 1))
    evcache = similar(mat)
    evs[1] = 0.1 # ensures that while runs at least ones
    while !whole_equal_ev(evs, equal_count)
        evs = random_try!(mat, evcache, possible_values)
    end
    convert(Matrix{Int64}, mat)
end

whole_equal_ev(::Vector{ComplexF64}, _) = false

function whole_equal_ev(evs, equal_count)
    c = collect(counter(evs))
    if whole_ev(evs)
        if equal_count == size(evs)[1]
            return size(c)[1] == 1
        end
        return equal_count ∈ [c[i][2] for i in 1:size(c)[1]]
    end
    false
end

# private functions

function random_try!(mat, evcache, possible_values)
    myrand!(mat, possible_values)
    copyto!(evcache, mat)
    eigvals!(evcache)
end

function myrand!(mat, possible_values)
    # overrides every entry of mat
    for i in eachindex(mat)
        mat[i] = rand(possible_values)
    end
end

