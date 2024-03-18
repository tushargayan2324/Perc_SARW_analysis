include("/home/tushar/julia_scripts/perc_question/initialisation_module.jl")
using .initialise_module
using StaticArrays
using Plots

mutable struct lat_site_astar
    site::lat_site
    g_score::Int
    h_score::Int
end


function find_clusters_trial_test_allocs_fix(lat)

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
        temp_arr = lat_site[]
        for i=1:size(lat)[1]
            if lat[i].cluster_id == intersect(left_side, right_side)[1]
                push!(temp_arr,lat[i])
            end
        end
        return temp_arr
    end

    return 0    
end


f = lattice_initialise(500,0.5)
g = find_clusters_trial_test_allocs_fix(f)
h = perc_cluster_elements(g)

cl_heatmap_plot(h,500)


function cl_heatmap_plot(f,n)
    # n = trunc(Int64,sqrt(size(f)[1]))
    # f = reshape(f, (n,n))
    Mat = zeros(Float64, (n,n))
	# for i=1:n
    #     for j=1:n
    #         Mat[i,j] = sin(11*f[j,i].cluster_id)^2
    #     end
    # end
    for obj in f
        s_id = obj.site_id
        Mat[x_coord(s_id,n),y_coord(s_id,n)] = obj.cluster_id #sin(11*f[j,i].cluster_id)^2
    end
    
    heatmap(Mat,size=(1500,1500), legend=false) #c=:prism tab20 glasbey_hv_n256

    Plots.savefig("hehe.png")
end


function manhattan_dist(x1,y1,x2,y2)
    return abs(x1-x2) + abs(y1-y2)
end

function h_score(lat_id1, lat_id2, n)
    # return abs(x1-x2) + abs(y1-y2)
    return manhattan_dist(x_coord(lat_id1,n), y_coord(lat_id1,n), 
                x_coord(lat_id2,n), y_coord(lat_id2,n))
end



manhattan_dist(1,2,3,4)

h_score(4,14,4)

function A_star_algo_percolating_cluster(lat,n, st_id, end_id) #input is the percolating cluster
    id_list = [obj.site_id for obj in lat ]

    st_index = findfirst(obj -> obj==st_id, id_list )
    
    r = size(lat)[1]

    g = 0

    open_list = lat_site_astar[lat_site_astar(lat[st_index],g,h_score(st_id,end_id,n))]
    
    closed_list = lat_site_astar[]

    crr_id = st_id

    while crr_id != end_id 
        crr_node = x_coord()
    
    end

    return 0
    
end

# Getting list of left starting points
# for i=1:r
#     if lat[i].cluster_id % n == 1
#         push!(open_list,[r,g])
#     end
# end
    


lat_site_astar(lat_site(1,[2,3],0), 4,5)

hehe = lattice_initialise(10,0.5)
gege = find_clusters_trial_test_allocs_fix(hehe)

map(x -> x.cluster_id, gege)


open_list = [[gege[1],10]]

push!(open_list, [gege[2], 10])

shee = rand(gege, 10)

findfirst(obj -> obj==88, [obj.site_id for obj in shee ])



# i=1
# while true
    

#     if i==10
#         break
#     end

#     i+=1

#     println(i)
#     sleep(0.3)
# end

