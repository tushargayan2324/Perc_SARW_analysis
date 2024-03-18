include("initialisation_module.jl")
using .initialise_module
using StaticArrays
using Plots


function my_my_plot(f)
	plt = scatter(0,0);
    n = trunc(Int64,sqrt(size(f)[1]))

	for obj in f
    	x1 = x_coord(obj.site_id,n)
    	y1 = y_coord(obj.site_id,n)
    		plt = scatter!((x1,y1),c=:prism,marker_z=obj.cluster_id, markersize=10,legend=false, cbar=false)#,axis=([], false))
    	for nei in obj.neighbours
        	    xs = [x1,x_coord(nei,n)]
        	    ys = [y1,y_coord(nei,n)]
        	    plt = plot!(xs,ys,color=:black)
    	end
	end
	# return plt
    savefig("perc_10_wbond.png")
end

function lattice_initialise_new(n, p)
    len = n*n
    lat = [lat_site(i,[],0) for i=1:len]

    for i=1:len
        if rand() < p 
            
            next_site = rand([i+n,i-n,i+1,i-1])
            prev_site = i


            for j=1:rand([4,5,6,7])
                if (1<x_coord(next_site,n)<n) && (1<y_coord(next_site,n)<n)
                
                # next_site = rand([i+n,i-n,i+1,i-1])
                # push!(lat[i].neighbours, next_site) ; push!(lat[next_site].neighbours, i) 
                # for j=1:5#rand([4,5,6,7])
                    #next_site = rand([i+n,i-n,i+1,i-1])
                    if !(next_site in lat[prev_site].neighbours) 
                        push!(lat[prev_site].neighbours, next_site) ; push!(lat[next_site].neighbours, prev_site)
                        
                    end
                    # 
                # end
                prev_site = next_site
                next_site = rand([next_site+n,next_site-n,next_site+1,next_site-1])    
                else
                    continue
                end
                
            end

        end
    end



    # for i=n+1:len # Giving top-bottom neighbours with some probability p
    #     if rand() < p
    #         push!(lat[i].neighbours,i-n)
    #         push!(lat[i-n].neighbours, i)
    #     end
    # end   

    # for i=1:len # Giving left-right neighbours with some probability p
    #     if i%n==0
    #         continue
    #     end
    #     if rand() < p
    #         push!(lat[i].neighbours,i+1)
    #         push!(lat[i+1].neighbours, i)
    #     end
    # end   

    return lat
end


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

n = 500

f = @time lattice_initialise(n,0.5);


h = @time find_clusters_trial_test_allocs_fix(f);


heatmap_plot(h)


my_my_plot(h)





function multithreaded_solve_half_horiz(lat)
    half_way = size(lat)[1]/2
    lat1 = lat[1:half_way]
    lat2 = lat[half_way:end]
    return 0
end