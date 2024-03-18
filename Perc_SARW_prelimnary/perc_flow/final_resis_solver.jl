using LinearAlgebra
using SparseArrays
include("initialisation_module.jl")
using .initialise_module

function lattice_initialise_bondnum(n, p)
    len = n*n
    lat = [lat_site(i,[],0) for i=1:len]

    bond_num = 0

    for i=n+1:len # Giving top-bottom neighbours with some probability p
        if rand() < p
            push!(lat[i].neighbours,i-n)
            push!(lat[i-n].neighbours, i)
            bond_num += 1
        end
    end   

    for i=1:len # Giving left-right neighbours with some probability p
        if i%n==0
            continue
        end
        if rand() < p
            push!(lat[i].neighbours,i+1)
            push!(lat[i+1].neighbours, i)
            bond_num += 1
        end
    end   

    return lat, bond_num
end

function return_site_pos(lat, site)
    for i=1:size(lat)[1]
        if lat[i].site_id == site
            return i
        end
    end
end



function lattice_to_matrix(lat,n,bond_num)
    # n = trunc(Int,sqrt(size(lat)[1])) # sqrt(size(lat)[1])

    lat_copy = deepcopy(lat)

    # mat = zeros(2*n^2 - 2*n , n^2)
    mat = zeros(bond_num, size(lat)[1])
    j = 1
    t = 1
    for i=1:size(lat)[1]
        site_neighbours = deepcopy(lat_copy[i].neighbours)

        if size(site_neighbours)[1] == 0
            continue
        end

        # println(i)
        
        for k=1:size(site_neighbours)[1]
            mat[j,t] = -1
            
            # r = return_site_pos(lat_copy, site_neighbours[k])

            # mat[j,site_neighbours[k]] = 1
            r = return_site_pos(lat_copy, site_neighbours[k])
            
            mat[j,r] = 1
            
            

            deleteat!(lat_copy[i].neighbours, findfirst(==(site_neighbours[k]), lat_copy[i].neighbours))
            
            # r = return_site_pos(lat_copy, site_neighbours[k])
            # println(r)
            deleteat!(lat_copy[r].neighbours, findfirst(==(lat_copy[i].site_id), lat_copy[r].neighbours))
            # push!(lat[i].neighbours,i-n)
            # push!(lat[i-n].neighbours, i)
            j+=1
            # println(i)
        end
        t += 1
    end
    return mat
end

function perc_prob(lat)
    left_side = Int[]
    right_side = Int[]
    top_side = Int[]
    bottom_side = Int[]

    n = trunc(Int, sqrt(size(lat)[1]))

    for i=1:n
        push!(left_side, lat[1 + (i-1)*n].cluster_id)
        push!(right_side, lat[i*n].cluster_id)
        push!(top_side, lat[n^2 - n + i].cluster_id)
        push!(bottom_side, lat[i].cluster_id)        
    end

    deleteat!(left_side, left_side.==0)
    deleteat!(right_side, right_side.==0)
    deleteat!(top_side, top_side.==0)
    deleteat!(bottom_side, bottom_side.==0)

    if size(intersect(left_side, right_side))[1] != 0 || size(intersect(top_side, bottom_side))[1] != 0
        return 1
    end

    return 0
end


function perc_cluster_elements(lat)
    left_side = Int[]
    right_side = Int[]
    top_side = Int[]
    bottom_side = Int[]

    n = trunc(Int, sqrt(size(lat)[1]))

    for i=1:n
        push!(left_side, lat[1 + (i-1)*n].cluster_id)
        push!(right_side, lat[i*n].cluster_id)
        push!(top_side, lat[n^2 - n + i].cluster_id)
        push!(bottom_side, lat[i].cluster_id)        
    end

    deleteat!(left_side, left_side.==0)
    deleteat!(right_side, right_side.==0)
    deleteat!(top_side, top_side.==0)
    deleteat!(bottom_side, bottom_side.==0)

    # if size(intersect(left_side, right_side))[1] != 0 && size(intersect(top_side, bottom_side))[1] != 0
    #     return ...
    # elseif size(intersect(left_side, right_side))[1] != 0
    #     ...
    # elseif size(intersect(top_side, bottom_side))[1] != 0
    #     ...
    # end

    # print(intersect(left_side, right_side))

    

    if size(intersect(left_side, right_side))[1] == 1
        # print("yess",intersect(left_side, right_side))
        temp_arr = []
        for i=1:size(lat)[1]
            if lat[i].cluster_id == intersect(left_side, right_side)[1]
                push!(temp_arr,lat[i])
            end
        end
        return temp_arr
    end

    return 0    
end

# m = [2]



for i=1:9
    if i%2==0
        continue
    end
    println(i)
end



h = lattice_initialise_bondnum(15,0.51)

f = find_clusters_trial_test_allocs_fix(h[1])

g = perc_cluster_elements(f)

q = lattice_to_matrix(g, 50, h[2])


# g = []

# for i in f
#     if i.cluster_id == 2
#         push!(g, i)
#     end
# end

# g

my_plot(f)
# heatmap_plot(f)



# f

# deleteat!(f, [100])

# f


A = q#lattice_to_matrix(g,5)

M = A' * A

shee = size(M)[1]

b = zeros(shee)
for i=1:50:2500
    # println(i)
    b[i] = 1
    b[i+49] = -1
end

# b = zeros(2500)
b[1:5] .= 1

b[end-4:end] .= -1

sum(b)

sum(M[1:end-1,1:end-1]*ones(size(M)[1]-1))

x = sparse(M.+1) \ b 

maximum(x)

x[1] - x[end]

x = sparse(M)\b

x[1] - x[end]



g = [(i,j) for i=1:5, j=5:8]

