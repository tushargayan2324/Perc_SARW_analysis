module initialise_module
    
    using ExportAll
    using Plots
    using StaticArrays


mutable struct lat_site
    const site_id::Int64
    neighbours::Array{Int64,1}
    cluster_id::Int64
end

function lattice_initialise(n, p)
    len = n*n
    lat = [lat_site(i,[],0) for i=1:len]

    for i=n+1:len # Giving top-bottom neighbours with some probability p
        if rand() < p
            push!(lat[i].neighbours,i-n)
            push!(lat[i-n].neighbours, i)
        end
    end   

    for i=1:len # Giving left-right neighbours with some probability p
        if i%n==0
            continue
        end
        if rand() < p
            push!(lat[i].neighbours,i+1)
            push!(lat[i+1].neighbours, i)
        end
    end   

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


function difference_array(X,Y) #return array with ele in X but not in Y
    for i in X
        if !(i in Y)
            continue
        end
        deleteat!(X, X .== i);
    end
    return X
end

function union_array_assigned(X,Y) #NOT A UNION LIKE SET!!!
    n = size(X)[1]+size(Y)[1]
    M = zeros(Int64, n)
    for i=1:size(X)[1]
        M[i] = X[i]        
    end
    for j=size(X)[1]+1:size(X)[1]+size(Y)[1]
        M[j] = Y[j-size(X)[1]]
    end
    return M
end



function union_array_push(X,Y)
    for i in Y
        if !(i in X)
            push!(X,i)
        end
    end
    return X
end

function del_element_first(arr, x)
    for i in arr
        if i!=x
            continue
        end
        deleteat!(arr, findfirst(==(i), arr))
    end

    return arr    
end

function mod_num(x,n)
    if x%n == 0
        return n
    else
        return x%n
    end
end



# function y_coord(y,n)
#     if y%n == 0
#         return trunc(Int,y/n) - 1
#     else
#         return trunc(Int,y/n)
#     end
# end

function x_coord(x,n)   # These elegant form suggested to me by @budbak
    return (x-1)%n + 1    
end


function y_coord(y,n)
    return trunc(Int64, (y-1)/n) + 1    
end


function heatmap_plot(f)
    n = trunc(Int64,sqrt(size(f)[1]))
    f = reshape(f, (n,n))
    Mat = zeros(Float64, (n,n))
	for i=1:n
        for j=1:n
            Mat[i,j] = sin(11*f[j,i].cluster_id)^2
        end
    end
    
    heatmap(Mat,c=:prism,size=(1500,1500), legend=false) #c=:prism tab20 glasbey_hv_n256

    Plots.savefig("hehe.png")
end




function my_plot(f)
	plt = scatter(0,0);
    n = trunc(Int64,sqrt(size(f)[1]))

	for obj in f
    	x1 = x_coord(obj.site_id,n)
    	y1 = y_coord(obj.site_id,n)
    		plt = scatter!((x1,y1), c=:glasbey_hv_n256,marker_z=sin(11*obj.cluster_id)^2, markersize=10,legend=false, cbar=false, size=(1000,600))#,axis=([], false))
    	for nei in obj.neighbours
        	    xs = [x1,x_coord(nei,n)]
        	    ys = [y1,y_coord(nei,n)]
        	    plt = plot!(xs,ys,color=:black, linewidth=5,size=(1000,600))
    	end
	end
	return plt
end



function my_plot_no_line(f)
	plt = scatter(0,0);
    n = trunc(Int64,sqrt(size(f)[1]))

	for obj in f
    	x1 = x_coord(obj.site_id,n)
    	y1 = y_coord(obj.site_id,n)
    		plt = scatter!((x1,y1), marker_z=obj.cluster_id, markersize=10,legend=false)#,axis=([], false))
#    	for nei in obj.neighbours
#        	    xs = [x1,mod_num(nei,n)]
#        	    ys = [y1,1+y_coord(nei,n)]
#        	    plt = plot!(xs,ys,color=:black)
#    	end
	end
	#return plt
	Plots.savefig("hehe.png")
end

    @exportAll

end
