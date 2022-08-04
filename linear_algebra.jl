
using LinearAlgebra
using AbstractAlgebra

function matrix_good_ev(dim, possible_values)
    evs = zeros(dim)
    A = zeros(dim,dim)
    while !good_ev(evs)
        A = rand(possible_values,dim,dim)
        evs = eigvals(A)
    end
    return A
end

function matrix_good_ev(dim, possible_values)
    mat = rand(dim, dim)
    matrix_good_ev!(mat, possible_values)
    mat
end

function matrix_good_ev!(mat::AbstractArray, possible_values)
    evs = zeros(size(mat, 1))
    evcache = similar(mat)
    evs[1] = 0.1 #ensure to run while at least 1
    while !good_ev(evs)
        evs = random_try!(mat, evcache, possible_values)
    end
    mat
end

function random_try!(mat, evcache, possible_values)
    myrand!(mat, possible_values)
    copyto!(evcache, mat)
    eigvals!(evcache)
end

function myrand!(mat, possible_values)
    for i in eachindex(mat)
        mat[i] = rand(possible_values)
    end
end



matrix_good_ev(dim, start, stop) = matrix_good_ev(dim, start:stop)

good_ev(::Vector{ComplexF64}) = false

function good_ev(evs)
    for ev in evs
        if ev ≠ floor(ev)
            return false
        end
    end
    return true
end

matrix_good_ev(2, 1, 5)


function char_poly(A) 
    B = matrix(ZZ, A);
    Zx, x = ZZ["x"]
    return charpoly(Zx,B)
end

function find_many_good_matrices(count,dim,range)
    good_matrices = fill(zeros(dim,dim),count)
    index = 1
    for i in 1:count
        good_matrix = matrix_good_ev(dim,range)
        if good_matrix ∉ good_matrices
            good_matrices[index] = good_matrix
            index += 1
        end

    end
    return good_matrices
end

function number_good_matrices(matrices) 
    count = 0
    for matrix in matrices
        if matrix != last(matrices)
            #display(matrix)
            count += 1
        end
    end
    return count
end
count = 100000
dim = 4
range = (0,1)
#A = matrix_good_ev(dim, (-3,3))
matrices = find_many_good_matrices(count,dim,range)
#println(matrices)
println(number_good_matrices(matrices))


