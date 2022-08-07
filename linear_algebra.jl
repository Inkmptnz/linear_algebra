using LinearAlgebra
using AbstractAlgebra
using Polynomials
# naive attempt - takes more space, is similar fast - eigvals takes longest
matrix_good_ev(dim, start, stop) = matrix_good_ev(dim, start:stop)
function matrix_good_ev(dim, possible_values)
    evs = zeros(dim)
    A = zeros(dim,dim)

    evs[1] = 0.1 #ensures that while runs at least ones
    while !good_ev(evs)
        A = rand(possible_values,dim,dim)
        evs = eigvals(A)
    end
    return A
end

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

# function format_char_pol(evs)
#     if good_ev(evs)
#         return [Polynomial([-ev,1], :λ) for ev ∈ evs]
#     else
#         good_evs = get_good_ev(evs)
#         return vcat([Polynomial([-ev,1], :λ) for ev ∈ good_evs], [Polynomial([-ev,1], :λ) for ev in evs if ev ∉ good_evs])
#     end
# end

# function get_good_ev(evs)
#     if good_ev(evs)
#         return evs
#     else
#         return [ev for ev in evs if good_ev(ev)]
#     end
# end

orthonormal_base(V::Matrix{Int64}) = orthonormal_base(convert(Matrix{Float64}, V))
function orthonormal_base(V::Matrix{Float64})
    (n, k) = size(V)
    U = zeros(n,k)
    U[:,1] = V[:,1]/norm(V[:,1])
    for i in 2:k
        U[:,i] = V[:,i]
        for j in 1:(i-1)
            U[:,i] = U[:,i] - (transpose(U[:,j])*U[:,i]) * U[:,j];
        end
        U[:,i] = U[:,i] / norm(U[:,i])
    end
    U
end


# count = 100000
# dim = 4
# range = (0,1)
# #A = matrix_good_ev(dim, (-3,3))
# matrices = find_many_good_matrices(count,dim,range)
# #println(matrices)
# println(number_good_matrices(matrices))


