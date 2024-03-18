using LinearAlgebra
using SparseArrays
using RowEchelon
using Plots

D1 = [ -1 1 0 0; 0 -1 1 0; 0 0 -1 1 ]


using SparseArrays
N = 3
D1 = sparse(I,N-1,N) - spdiagm(N-1,N, 1=>ones(N-1))
D = [ kron(D1, sparse(I,N,N)); kron(sparse(I,N,N), D1) ]
i, j = 1 , 9
b = zeros(N^2); b[i], b[j] = 1, -1
v = (D' * D) \ b
v[i] - v[j]

D'*D


A = Float64[-1 1 0 0;
-1 0 1 0;
0 -1 1 0;
-1 0 0 1;
0 -1 0 1;
0 0 -1 1 ]

h = lattice_initialise(10,1)

f = find_clusters_trial_test_allocs_fix(h)

# heatmap_plot(f)
my_plot(f)
# savefig("resis_net_10x10.png")


HUI = lattice_to_matrix(h)


A = HUI


A = Float64[ 1 -1 0 0 0 0;
            0 1 -1 0 0 0;
            0 1 0 0 -1 0;
            0 0 0 1 -1 0;
            0 0 0 0 1 -1
            ]



M = A' * A

# M_n = M[1:8,1:8]

M*ones(6)

b = zeros(6)
for i=1:25:625
    # println(i)
    b[i] = 1
    b[i+24] = -1
end
b[1] = 1
b[end] = -2
b[3] = 1
b[4] = -1

b = zeros(100)
b[12] = 1
b[77] = -1


b



x = sparse(M) \ b #[ 1 ;0 ;0 ]

x[6] - x[1]

A*x


A* push!(x,0)550


# eigenvector()



function create_random_resistor_matrix(n,p) #size is the input
    rd_lat = zeros(2*n^2 - 2*n , n^2)
    for i=1:r
        rd_lat[]
    end    
end


zeros(4,3)

h = lattice_initialise(100,0.5)

f = find_clusters_trial_test_allocs_fix(h)

heatmap_plot(f)


t= deepcopy(h)

push!(t[1].neighbours, 2)

t

h


function lattice_to_matrix(lat)
    n = trunc(Int,sqrt(size(lat)[1])) # sqrt(size(lat)[1])

    lat_copy = deepcopy(lat)

    mat = zeros(2*n^2 - 2*n , n^2)
    j = 1
    for i=1:n^2
        site_neighbours = deepcopy(lat_copy[i].neighbours)
        for k=1:size(site_neighbours)[1]
            mat[j,i] = -1; mat[j,site_neighbours[k]] = 1
            deleteat!(lat_copy[i].neighbours, findfirst(==(site_neighbours[k]), lat_copy[i].neighbours))
            deleteat!(lat_copy[site_neighbours[k]].neighbours, findfirst(==(i), lat_copy[site_neighbours[k]].neighbours))
            # push!(lat[i].neighbours,i-n)
            # push!(lat[i-n].neighbours, i)
            j+=1
            println(i)
        end
    end
    return mat
end

HUI = lattice_to_matrix(h)


# B = [1 3 3 2;2 6 9 7;-1 -3 3 4]

# C = [4 3; 6 3]

# rref(A)[:,:]

# L,U,p = lu(A)

# L,U,p = lu(B,check=false)

# factorize(A)


# hehe = lu(A)

# lu(A,check =false)


# sparse(I, 3, 3)

# (2 => ones(4))

# sparse(I,10,10)

# kron([1 2 ; 3 4 ], [5 6 ; 7 8])