
using LinearAlgebra
using AbstractAlgebra

function matrix_good_ew(dim, range)
    evs = zeros(n)
    evs[1] = 0.1 #ensure to run while at least 1
    A = zeros(dim,dim)
    while !good_ew(evs)
        A = rand(range[1]:range[2],dim,dim)
        evs = eigvals(A)
    end
    return A
end

function good_ew(evs)
    if isa(evs,Vector{ComplexF64})
        return false
    end
    for ev in evs
        if ev ≠ floor(ev)
            return false
        end
    end
    return true
end


function char_poly(A) 
    B = matrix(ZZ, A);
    Zx, x = ZZ["x"]
    return charpoly(Zx,B)
end

function find_many_good_matrices(count,dim,range)
    good_matrices = fill(zeros(dim,dim),count)
    index = 1
    for i in 1:count
        good_matrix = matrix_good_ew(dim,range)
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
count = 10000
dim = 4
range = (0,1)
#A = matrix_good_ew(dim, (-3,3))
matrices = find_many_good_matrices(count,dim,range)
#println(matrices)
println(number_good_matrices(matrices))


