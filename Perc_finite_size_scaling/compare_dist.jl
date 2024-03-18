include("initialisation_module.jl")
using .initialise_module
using StaticArrays
using Plots
using StatsBase
using Random




function p_value_perc(n::Int, rd_len) #inputs are size and prob of having those objs
    len = n*n
    lat = [lat_site(i,[],0) for i=1:len]

    num_bonds = 0

    num_rdws = 0

    # rd_len = 90 #trunc(Int,len*p)

    cl_id = 0

    lat_copy = find_cluster_latest(deepcopy(lat))
    while perc_prob(lat_copy) != 1

        ## Getting the list of sites who are not already used for cluster
        req_arr = Int64[]

        for j=1:len
            if lat[j].cluster_id == 0 
                push!(req_arr, j)
            end
        end

        rd_site = rand(req_arr)

        rd_walk_obj = gen_SARW_fix_length_new(rd_len, rd_site, n)        

        for i=2:rd_len
            if !(rd_walk_obj[i] in lat[rd_walk_obj[i-1]].neighbours) 
                push!(lat[rd_walk_obj[i-1]].neighbours, rd_walk_obj[i]) ; 
                push!(lat[rd_walk_obj[i]].neighbours, rd_walk_obj[i-1])        
                num_bonds += 1
                # len_rdw += 1
            end
        end
    
        num_rdws += 1
        # println(num_rdws)
        # sleep(0.01)
        lat_copy = find_cluster_latest(deepcopy(lat))
    end

    for i=1:len
        lat[i].cluster_id = 0
    end

   return lat , num_bonds, num_rdws
end

function find_cluster_latest_modified(lat)

    arr_sites = Array{Int64}([ i for i=1:size(lat)[1] ])

    for i=1:size(lat)[1]
        lat[i].cluster_id = 0
    end

    cl_id = 0

    # for i=1:size(lat)[1]
    #     if lat[i].cluster_id > cl_id
    #         cl_id = lat[i].cluster_id
    #     end
    # end


    while length(arr_sites) != 0 
        cl_id += 1
                
        start_site = arr_sites[1]

        deleteat!(arr_sites, 1)

        lat[start_site].cluster_id = cl_id #set cl_id of first element in cluster list
        
        track = size(lat[start_site].neighbours)[1] #tracking the length of neighbours set using this variable 
        
        neighbour_arr = SVector{track,Int64}(lat[start_site].neighbours)

        while track!=0      #terminate loop when can't find any viable next site for cluster
            neighbours_set = zeros(Int64, 2*size(neighbour_arr)[1]+4)
            iter_track::Int64 = 0
            for obj in neighbour_arr
                lat[obj].cluster_id = cl_id

                
                for nnn in lat[obj].neighbours
                    
                    if lat[nnn].cluster_id == 0 && !(nnn in neighbours_set)
                            iter_track += 1
                            neighbours_set[iter_track] = nnn
                    end
                end

                deleteat!(arr_sites, findfirst(==(obj), arr_sites))                
            end      
            deleteat!(neighbours_set, neighbours_set.==0)
            track = size(neighbours_set)[1]
            neighbour_arr = SVector{track,Int64}(neighbours_set)
        end
    end
    return lat
end


function find_p_value_dist(ens)
    n = 50
    p_value_arr = zeros(ens)
    num_rdws_arr = zeros(ens)
    
    Threads.@threads for i=1:ens
        h = p_value_perc(n,35)
        # push!(p_value_arr,h[])
        p_value_arr[i] = h[2]/(2*n^2-2*n)
        num_rdws_arr[i] = h[3]

        # println(i)
    end

    return p_value_arr, num_rdws_arr
end




dat = @time find_p_value_dist(10^4)


dat_1 = @time find_p_value_dist(5*10^3)


dat_2 = @time find_p_value_dist(10^4)

dat_3 = @time find_p_value_dist(10^4)

open("geek.txt", "w") do file
    write(file, dat_1)
end

open("datfile_150x150_105l_10^4.txt","a") do io
    println(io,dat)
end

open("datfile_100x100_75l_10^4.txt","a") do io
    println(io,dat_2)
end

open("datfile_50x50_35l_10^4.txt","a") do io
    println(io,dat_3)
end


p_vl_1 = vcat(dat_1[1], dat_1[1])

Plots.scalefontsizes(1/1.2)

histogram(dat[1], normalize=true, label="150x150; l=105", alpha=0.5, size=(2200,1500))
histogram!(dat_2[1], normalize=true, label="100x100; l=75", alpha=0.5, size=(2200,1500))
histogram!(dat_3[1], normalize=true, label="50x50; l=50", alpha=0.5, size=(2200,1500))


xlabel!("p_value")
ylabel!("distribution")
title!("distribution of p values when first percolates")



p_vl_2 = vcat(dat_1[2], dat_1[2])


histogram(p_vl_2,normalize=false)


histogram(dat[2], normalize=true, label="150x150; l=105", alpha=0.5, size=(2200,1500))
histogram!(dat_2[2], normalize=true, label="100x100; l=75", alpha=0.5, size=(2200,1500))
histogram!(dat_3[2], normalize=true, label="50x50; l=50", alpha=0.5, size=(2200,1500))




xlabel!("number of SARWS")
ylabel!("distribution")
title!("distribution of #SARWS when first percolates")



savefig("p_value_dist_compare.png")
savefig("num_rdws_hist_compare.png")


for i=1:10
    h = lattice_initialise(300,0.5)

    g = find_cluster_latest(h)
    
    l = find_cluster_latest_modified(g)
end



h = lattice_initialise(300,0.5)

g = find_cluster_latest(h)

h[1].neighbours = [2]

g

h

l = find_cluster_latest_modified(g)


my_plot(g)

g[2]

push!(g[12].neighbours, 22)
push!(g[22].neighbours, 12)

g[3]

l = find_cluster_latest_modified(g)

my_plot(l)


