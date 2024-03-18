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

function find_cluster_latest(lat)

    arr_sites = Array{Int64}([i for i=1:(size(lat)[1])])
    cl_id = 0

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


function find_cluster_latest_modified(lat)

    # arr_sites = Array{Int64}([i for i=1:(size(lat)[1])])
    arr_sites = Int64[]
    for i=size(lat)[1]
        if lat[i].cluster_id == 0
            push!(arr_sites, lat[i].site_id)
        end
    end

    cl_id = 0

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
	end
	#return plt
	Plots.savefig("hehe.png")
end


function next_site_NPBC(site::Int, n::Int)
    if (1<x_coord(site,n)<n) && (1<y_coord(site,n)<n)
        # return rand([site+n, site-n, site+1, site-1])
        # println("1")
        return [site+n, site-n, site+1, site-1]


    elseif (1<y_coord(site,n)<n) && (x_coord(site, n) == 1 || x_coord(site, n) == n) 
        # return rand([ site+n, site-n, site - sign( x_coord(site,n) - trunc(Int, (n+1)/2) ) ])
        # println("2")
        return [ site+n, site-n, site - sign( x_coord(site,n) - trunc(Int, (n+1)/2) ) ]


    elseif (1<x_coord(site,n)<n) && ( y_coord(site, n) == 1 || y_coord(site, n) == n )
        # return rand([ site+1, site-1, site - n*sign( y_coord(site,n) - trunc(Int, (n+1)/2) ) ])
        # println("3")
        return [ site+1, site-1, site - n*sign( y_coord(site,n) - trunc(Int, (n+1)/2) ) ]


    else
        # return rand([site - sign( x_coord(site,n) - trunc(Int, (n+1)/2) ) , site - n*sign( y_coord(site,n) - trunc(Int, (n+1)/2) ) ])
        # println("4")
        return [site - sign( x_coord(site,n) - trunc(Int, (n+1)/2) ) , site - n*sign( y_coord(site,n) - trunc(Int, (n+1)/2) ) ]
    end
end

function cl_heatmap_plot(f,n)
    # n = trunc(Int64,sqrt(size(f)[1]))
    Mat = zeros(Float64, (n,n))
    for obj in f
        s_id = obj.site_id
        Mat[x_coord(s_id,n),y_coord(s_id,n)] = obj.cluster_id #sin(11*f[j,i].cluster_id)^2
    end
    
    heatmap(Mat,size=(1500,1500), legend=false) #c=:prism tab20 glasbey_hv_n256

    Plots.savefig("hehe.png")
end



function gen_SARW_fix_length_new(max_itr, start_site, n) #Generating and getting back the length of generated SARW
    # str_site = [0,0]
    
    next_site = copy(start_site)

    sites_arr = zeros(Int, max_itr)#[0 for i=1:max_itr]

    sites_arr[1] = start_site

    k = true
    while k==true
        for i=2:max_itr
            prev_site = next_site        
            next_site_arr = next_site_NPBC(prev_site, n)

            choose_site_arr = []

            for s in next_site_arr
                if !(s in sites_arr)
                    push!(choose_site_arr, s)
                end
            end

            if size(choose_site_arr)[1] == 0
                next_site = copy(start_site)
                sites_arr = zeros(Int, max_itr)            
                sites_arr[1] = start_site
                break
            elseif i==max_itr
                k=false
            end

            next_site = rand(choose_site_arr)
            sites_arr[i] = next_site
        end
        
    end


    return sites_arr
end


function initialise_lattice_sarw_fix_len_while(n::Int, p, rd_len) #inputs are size and prob of having those objs
    len = n*n
    lat = [lat_site(i,[],0) for i=1:len]

    num_bonds = 0

    num_rdws = p#trunc(Int,len*p)

    # rd_len = 90 #trunc(Int,len*p)

    cl_id = 0
    for i=1:num_rdws

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
                push!(lat[rd_walk_obj[i-1]].neighbours, rd_walk_obj[i]) ; push!(lat[rd_walk_obj[i]].neighbours, rd_walk_obj[i-1])        
                num_bonds += 1
                # len_rdw += 1
            end
        end
    end

    for i=1:len
        lat[i].cluster_id = 0
    end

   return lat , num_bonds
end

# function initialise_lattice_sarw_fix_len(n::Int, p, rd_len) #inputs are size and prob of having those objs
#     len = n*n
#     lat = [lat_site(i,[],0) for i=1:len]

#     num_bonds = 0

#     num_rdws = p#trunc(Int,len*p)

#     # rd_len = 90 #trunc(Int,len*p)


#     cl_id = 0
#     for i=1:num_rdws

#         ## Getting the list of sites who are not already used for cluster
#         req_arr = Int64[]
#         for j=1:len
#             if lat[j].cluster_id == 0 
#                 push!(req_arr, j)
#             end
#         end

#         rd_site = rand(req_arr)

#         len_rdw = 0

#         prev_site = rd_site
#         next_site = rand(next_site_NPBC(rd_site, n)) #write a stupid fn for this
        
#         cl_id += 1
#         lat[next_site].cluster_id = cl_id
#         lat[prev_site].cluster_id = cl_id

#         for k=1:rd_len

#             if !(next_site in lat[prev_site].neighbours) 
#                 push!(lat[prev_site].neighbours, next_site) ; push!(lat[next_site].neighbours, prev_site)        
#                 num_bonds += 1
#                 len_rdw += 1
#             end
#                 # 
#             # end
#             prev_site = next_site
#             next_site = next_site_NPBC(next_site,n) #rand([next_site+n,next_site-n,next_site+1,next_site-1])    
#             # deleteat!(next_site, neighbours_set.==0)
#             pick_arr = Int[]

            
#             for obj in next_site
#                 if lat[obj].cluster_id != cl_id
#                     #deleteat!(next_site, findfirst(==(obj), next_site)) 
#                     push!(pick_arr, obj )
#                 end
#             end

#             if size(pick_arr)[1] == 0
#                 break
#             end

        
            
#             next_site = rand(pick_arr)

#             lat[next_site].cluster_id = cl_id

#         end



#     end

#     for i=1:len
#         lat[i].cluster_id = 0
#     end

#    return lat , num_bonds
# end

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

function perc_prob_ens(n,num_ens)
    n_p = 30 # 0-1 in how many steps
    prob_avg_ens = zeros(n_p)
    p = LinRange(0,1,n_p)

    Threads.@threads for j=1:n_p #LinRange(0,1,n_p)
        for i=1:num_ens
            f = find_clusters_trial_test_allocs_fix(lattice_initialise(n,p[j]))
            prob_avg_ens[j] += perc_prob(f)
        end
        println(j)
    end

    return prob_avg_ens/num_ens
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





function S_plotting()

    N_p = 100 # 50 values of p
    S_values = zeros(Float64, N_p)

    for i=1:N_p
    # for i=1:N_p
        S = 0
        for j=1:5
            n = 200
            lat = find_cluster_latest(initialise_lattice_sarw_fix_len_while(200,90+i,150)[1])
            
            my_data = countmap([obj.cluster_id for obj in lat])#distribution_cluster_size(lat)
            #n_s = collect(values(my_data))  
            # s = collect(keys(my_data)) 

            temp = countmap(collect(values(my_data)))

            n_s = collect(values(temp))/(n^2)

            s = collect(keys(temp))/(n^2)


            S += sum(s.^2 .* n_s) / ( sum(s .* n_s) )
        end
        S_values[i] =  S/5
        println(100*i/N_p,"%")
    end
    
    plot(LinRange(0,1,N_p),S_values )
    savefig("s_value_perc.png")
end

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




    @exportAll

end
