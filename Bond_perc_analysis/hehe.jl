include("initialisation_module.jl")
using .initialise_module
using StaticArrays
using Plots
using StatsBase

function find_clusters_trial_test_allocs_fix(lat)  # Hmmmm...Idk this solved the major problem of allocs (at high lattice number). Can I do anymore optimisation?!

    arr_sites = Array{Int64}([i for i=1:(size(lat)[1])])
    cl_id = 0

    while length(arr_sites) != 0 
        cl_id += 1
                
        start_site = arr_sites[1]

        deleteat!(arr_sites, 1)#findfirst(==(start_site), arr_sites))

        lat[start_site].cluster_id = cl_id #set cl_id of first element in cluster list
        
        track = size(lat[start_site].neighbours)[1] #tracking the length of neighbours set using this variable 
        
        neighbour_arr = SVector{track,Int64}(lat[start_site].neighbours)

        while track!=0      #terminate loop when can't find any viable next site for cluster

            # neighbours_set = Int64[] #next neighbour

            neighbours_set = zeros(Int64, 2*size(neighbour_arr)[1]+4)
            # println((neighbour_arr))
            iter_track::Int64 = 0

            for obj in neighbour_arr
                lat[obj].cluster_id = cl_id

                
                for nnn in lat[obj].neighbours
                    
                    if lat[nnn].cluster_id == 0 && !(nnn in neighbours_set)
                        # lat[obj].cluster_id = cl_id
                        # if !(nnn in neighbours_set)
                            #push!(neighbours_set, nnn)
                            
                            # println(iter_track)
                            iter_track += 1

                            neighbours_set[iter_track] = nnn
                            
                            # println(neighbours_set)
                        # end
                    end
                end

                
                # println(neighbours_set)

                #arr_sites = deleteat!(arr_sites, findfirst(==(obj), arr_sites))
                # if obj in arr_sites
                # println(arr_sites, " ", obj)
                deleteat!(arr_sites, findfirst(==(obj), arr_sites))
                # end

                # sleep(0.01)
                # println(length(arr_sites))
                
            end            
            

            deleteat!(neighbours_set, neighbours_set.==0)
            track = size(neighbours_set)[1]
            # println(track)
            neighbour_arr = SVector{track,Int64}(neighbours_set)
   
        end

    end
    return lat
end
#############################################################################################


n = 10

f = @time lattice_initialise(n,0.1);


h = @time find_clusters_trial_test_allocs_fix(f);

my_plot(h)


heatmap_plot(h)


cluster_size_freq(h)


#S is the cluster size, hence the following is plot of avg cluster size over all range of p values
function S_plotting()

    N_p = 100 # 50 values of p
    S_values = zeros(Float64, N_p)

    Threads.@threads for i=1:N_p
    # for i=1:N_p
        S = 0
        for j=1:5
            n = 200
            lat = find_clusters_trial_test_allocs_fix(lattice_initialise(n,i/N_p))
            
            my_data = countmap([obj.cluster_id for obj in lat])#distribution_cluster_size(lat)
            #n_s = collect(values(my_data))  
            # s = collect(keys(my_data)) 

            temp = countmap(collect(values(my_data)))

            n_s = collect(values(temp))/(n^2)

            s = collect(keys(temp))/(n^2)


            S += sum(s.^2 .* n_s) / ( sum(s .* n_s) )
        end
        S_values[i] =  S/5
        #println(100*i/N_p,"%")
    end

    #g = distribution_cluster_size(find_clusters_new_algo_not_condi(lattice_initialise(n,0.6)))    

    #S =  sum(s.^2 .* n_s) / ( sum(s .* n_s) )

    
    plot(LinRange(0,1,N_p),S_values )
    savefig("s_value_perc.png")
end

@time S_plotting()


#following is the plot of cluster sizes freq for a given value of p
function cluster_size_freq(h)
    cluster_size = [obj.cluster_id for obj in h]
    
    my_data = countmap(cluster_size)
    s = collect(values(my_data))
    
    temp = sort(countmap(s))
        
    n_s = collect(values(temp))
    
    s_size = collect(keys(temp)) 
        
    scatter(log.(s_size), log.(n_s))
    savefig("freq_vs_clustersize.png")        
end  

cluster_size_freq(h)
