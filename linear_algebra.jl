using LinearAlgebra
using AbstractAlgebra

# naive attempt - takes more space, is similar fast - eigvals takes longest
# matrix_good_ev(dim, start, stop) = matrix_good_ev(dim, start:stop)
# function matrix_good_ev(dim, possible_values)
#     evs = zeros(dim)
#     A = zeros(dim,dim)

#     evs[1] = 0.1 #ensures that while runs at least ones
#     while !good_ev(evs)
#         A = rand(possible_values,dim,dim)
#         evs = eigvals(A)
#     end
#     return A
# end

#optimized version

matrix_good_ev!(dim, start, stop) = matrix_good_ev!(dim, start:stop)

function matrix_good_ev!(dim, possible_values)
    mat = rand(dim, dim)
    matrix_good_ev!(mat, possible_values)
end

function matrix_good_ev!(mat::AbstractArray, possible_values)
    # mat allocates memory
    evs = zeros(size(mat, 1))
    evcache = similar(mat) # chache mat
    evs[1] = 0.1 # ensures that while runs at least ones
    while !good_ev(evs)
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

good_ev(::Vector{ComplexF64}) = false
function good_ev(evs)
    for ev in evs
        if ev ≠ floor(ev)
            return false
        end
    end
    return true
end

# uncomplete
function find_multiple_good_matrices(count, dim, range)
    good_matrices = fill(zeros(dim, dim), count)
    index = 1
    for i in 1:count
        good_matrix = matrix_good_ev!(dim, range)
        if good_matrix ∉ good_matrices
            good_matrices[index] = good_matrix
            index += 1
        end

    end
    good_matrices
end

function number_good_matrices(matrices) 
    count = 0
    for matrix in matrices
        if matrix != last(matrices)
            # display(matrix)
            count += 1
        end
    end
    return count
end

function char_poly(A) 
    B = matrix(ZZ, A);
    Zx, x = ZZ["x"]
    return charpoly(Zx,B)
end
# count = 100000
# dim = 4
# range = (0,1)
# #A = matrix_good_ev(dim, (-3,3))
# matrices = find_many_good_matrices(count,dim,range)
# #println(matrices)
# println(number_good_matrices(matrices))


