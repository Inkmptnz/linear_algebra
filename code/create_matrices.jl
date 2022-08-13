using LinearAlgebra
using AbstractAlgebra
using Polynomials
using DataStructures

# naive attempt - takes more space, is similar fast - eigvals takes longest
matrix_whole_ev(dim, start, stop) = matrix_whole_ev!(dim, start:stop, false)
matrix_whole_ev(dim, possible_values) = matrix_whole_ev!(dim, possible_values, false)
matrix_whole_ev(dim, start, stop, with_i) = matrix_whole_ev!(dim, start:stop, with_i)

function matrix_whole_ev(dim, possible_values, with_i)
    evs = zeros(dim)
    A = zeros(dim,dim)

    evs[1] = 0.1 #ensures that while runs at least ones
    while !whole_ev(evs, with_i)
        A = rand(possible_values,dim,dim)
        evs = eigvals(A)
    end
    return A
end

#optimized version
matrix_whole_ev!(dim, start, stop) = matrix_whole_ev!(dim, start:stop, false)
matrix_whole_ev!(dim, possible_values) = matrix_whole_ev!(dim, possible_values, false)
matrix_whole_ev!(dim, start, stop, with_i) = matrix_whole_ev!(dim, start:stop, with_i)

function matrix_whole_ev!(dim, possible_values, with_i)
    mat = rand(dim, dim)
    matrix_whole_ev!(mat, possible_values, with_i)
end

function matrix_whole_ev!(mat::AbstractArray, possible_values, with_i)
    # mat allocates memory
    evs = zeros(size(mat, 1))
    evcache = similar(mat) # chache mat
    evs[1] = 0.1 # ensures that while runs at least ones
    while !whole_ev(evs, with_i)
        evs = random_try!(mat, evcache, possible_values)
    end
    convert(Matrix{Int64}, mat)
end

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

whole_ev(evs) = whole_ev(evs, false)

function whole_ev(evs::Vector{ComplexF64}, with_i)
    if !with_i
        return false
    end
    for ev in evs
        if imag(ev) ≠ floor(imag(ev))
            return false
        else
            if real(ev) ≠ floor(real(ev)) 
                return false
            end
        end
    end
    return true
end

function whole_ev(evs, _)
    for ev in evs
        if ev ≠ floor(ev)
            return false
        end
    end
    return true
end

matrix_whole_imag_ev!(dim,start,stop) = matrix_whole_imag_ev!(dim, start:stop)
matrix_whole_imag_ev!(mat::AbstractArray,start,stop) = matrix_whole_imag_ev!(mat, start:stop)

function matrix_whole_imag_ev!(dim, possible_values)
    mat = rand(dim, dim)
    matrix_whole_imag_ev!(mat, possible_values)
end

function matrix_whole_imag_ev!(mat::AbstractArray, possible_values)
    evs = zeros(size(mat, 1))
    evcache = similar(mat)
    evs[1] = 0.1 # ensures that while runs at least ones
    while !whole_imag_ev(evs)
        evs = random_try!(mat, evcache, possible_values)
    end
    convert(Matrix{Int64}, mat)
end

whole_imag_ev(::Vector{Float64}) = false

function whole_imag_ev(evs)
    for ev in evs
        if imag(ev) ≠ floor(imag(ev))
            return false
        end
    end
    true
end


matrix_multiple_whole_ev!(dim,start,stop,equal_count) = matrix_whole_imag_ev!(dim, start:stop, equal_count)
matrix_multiple_whole_ev!(mat::AbstractArray,start,stop,equal_count) = matrix_whole_imag_ev!(mat, start:stop, equal_count)

function matrix_multiple_whole_ev!(dim, possible_values, equal_count)
    mat = rand(dim, dim)
    matrix_whole_imag_ev!(mat, possible_values, equal_count)
end

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
