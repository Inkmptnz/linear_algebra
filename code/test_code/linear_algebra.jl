using LinearAlgebra, AbstractAlgebra, Polynomials

include("create_matrices.jl")
# uncomplete

function daily_math()
    dims = Dict(2 => -6:6, 3 => -4:4, 4 => -2:2, 5 => 0:2)
    println("First Matrices without Complex Numbers.")
    good_matrices = sort([matrix_whole_ev!(dim[1],dim[2]) for dim in dims ],
    lt = (x,y) -> size(x) < size(y))
    print_math(good_matrices, false)

    println("Good job, now you get only matrices with imaginary numbers.")
    good_matrices = sort([matrix_whole_imag_ev!(dim[1],dim[2]) for dim in dims],
    lt = (x,y) -> size(x) < size(y))
    print_math(good_matrices, true)
end

function print_math(matrices, with_i)#
    if  with_i
        for mat in matrices
            display(mat)
            while true
                print("Want to see the results? (y): ")
                input = readline()
                if input == "y"
                    print_char_ev(mat)
                    println("Would you like to continue? (press anything)")
                    input = readline()
                    break
                end
            end
        end
    end
    for mat in matrices
        display(mat)
        while true
            println("What are your Eigenvalues?")
            input = readline()
            if input == "cheat"
                print_char_ev(mat)
                break
            end
            evs = [parse(Int16, i) for i in split(input, ",")]
            if sort(evs) == sort(eigvals(mat))
                println("Nice, now the next!")
                break
            else
                println("Try again!")
            end
        end
    end
end  

function find_multiple_good_matrices(count, dim, range)
    good_matrices = fill(zeros(dim, dim), count)
    index = 1
    for i in 1:count
        good_matrix = matrix_whole_ev!(dim, range)
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

function print_char_ev(A)
    wolframAlpha_copy(A)
    println(char_poly(A))
    println(eigvals(A))
end

function wolframAlpha_copy(A)
    string = "{"
    for i in 1:(size(A)[1])
        string *= "{"
        row = A[i,:]
        for j in row[1:end-1]
            string *= "$j,"
        end
        string *= "$(row[end])},"
    end
    println(string[1:(end-1)] * "}")
end
# function format_char_pol(evs)
#     if whole_ev(evs)
#         return [Polynomial([-ev,1], :λ) for ev ∈ evs]
#     else
#         whole_evs = get_whole_ev(evs)
#         return vcat([Polynomial([-ev,1], :λ) for ev ∈ whole_evs], [Polynomial([-ev,1], :λ) for ev in evs if ev ∉ whole_evs])
#     end
# end

# function get_whole_ev(evs)
#     if whole_ev(evs)
#         return evs
#     else
#         return [ev for ev in evs if whole_ev(ev)]
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
# #A = matrix_whole_ev(dim, (-3,3))
# matrices = find_many_good_matrices(count,dim,range)
# #println(matrices)
# println(number_good_matrices(matrices))


