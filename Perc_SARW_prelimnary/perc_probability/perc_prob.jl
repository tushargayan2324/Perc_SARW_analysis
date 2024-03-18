include("initialisation_module.jl")
using .initialise_module
using StaticArrays
using Plots
using StatsBase
using Random

function my_my_plot(f)
	plt = scatter(0,0);
    n = trunc(Int64,sqrt(size(f)[1]))

	for obj in f
    	x1 = x_coord(obj.site_id,n)
    	y1 = y_coord(obj.site_id,n)
    		plt = scatter!((x1,y1),size=(800,800),c=:prism,marker_z=obj.cluster_id, markersize=12,legend=false, cbar=false)#,axis=([], false))
    	for nei in obj.neighbours
        	    xs = [x1,x_coord(nei,n)]
        	    ys = [y1,y_coord(nei,n)]
        	    plt = plot!(xs,ys,size=(800,800),color=:black)
    	end
	end
	# return plt
    savefig("perc_10_wbond.png")
end


# This is not self avoiding random walks. I have let all the objects be independent of 
# each other. Can do the former one too, but its a pain.

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


# function lattice_initialise_new_sarw(n, p)
#     len = n*n
#     lat = [lat_site(i,[],0) for i=1:len]

#     cl_id = 0
#     for i=1:len
#         if rand() < p 
            
#             next_site = rand([i+n,i-n,i+1,i-1])
#             prev_site = i

#             cl_id += 1
#             lat[prev_site].cluster_id = cl_id 


#             for j=1:rand([4,5,6,7])
                
#                 if (1<x_coord(next_site,n)<n) && (1<y_coord(next_site,n)<n) 
                
#                 # next_site = rand([i+n,i-n,i+1,i-1])
#                 # push!(lat[i].neighbours, next_site) ; push!(lat[next_site].neighbours, i) 
#                 # for j=1:5#rand([4,5,6,7])
#                     #next_site = rand([i+n,i-n,i+1,i-1])
#                     if !(next_site in lat[prev_site].neighbours) 
#                         push!(lat[prev_site].neighbours, next_site) ; push!(lat[next_site].neighbours, prev_site)
                        
#                     end
#                     # 
#                 # end
#                 prev_site = next_site
#                 next_site = rand([next_site+n,next_site-n,next_site+1,next_site-1])    
#                 else
#                     continue
#                 end
                
#             end

#         end
#     end

#     return lat
# end


# m = Int[]
# push!(m,9)

# trunc(Int, 144*0.1 )


# Int == Int64

# m = 8^20

# typeof(m) == Int

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

# for i=1:16
#     println(next_site_NPBC(i, 4))
# end


# sign( y_coord(site,n) - trunc(Int, n/2) )

# y_coord(7,3)



function initialise_lattice_sarw(n::Int, p, rd_len) #inputs are size and prob of having those objs
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

        prev_site = rd_site
        next_site = rand(next_site_NPBC(rd_site, n)) #write a stupid fn for this
        
        cl_id += 1
        lat[next_site].cluster_id = cl_id
        lat[prev_site].cluster_id = cl_id

        for k=1:rd_len

            if !(next_site in lat[prev_site].neighbours) 
                push!(lat[prev_site].neighbours, next_site) ; push!(lat[next_site].neighbours, prev_site)        
                num_bonds += 1
            end
                # 
            # end
            prev_site = next_site
            next_site = next_site_NPBC(next_site,n) #rand([next_site+n,next_site-n,next_site+1,next_site-1])    
            # deleteat!(next_site, neighbours_set.==0)
            pick_arr = Int[]

            
            for obj in next_site
                if lat[obj].cluster_id != cl_id
                    #deleteat!(next_site, findfirst(==(obj), next_site)) 
                    push!(pick_arr, obj )
                end
            end

            if size(pick_arr)[1] == 0
                break
            end
        
            
            next_site = rand(pick_arr)

            lat[next_site].cluster_id = cl_id

        end



    end

    for i=1:len
        lat[i].cluster_id = 0
    end

   return lat , num_bonds
end

function initialise_lattice_rw(n::Int, p, rd_len) #inputs are size and prob of having those objs
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

        prev_site = rd_site
        next_site = rand(next_site_NPBC(rd_site, n)) #write a stupid fn for this
        
        cl_id += 1
        lat[next_site].cluster_id = cl_id
        lat[prev_site].cluster_id = cl_id

        for k=1:rd_len

            if !(next_site in lat[prev_site].neighbours) 
                push!(lat[prev_site].neighbours, next_site) ; push!(lat[next_site].neighbours, prev_site)        
                num_bonds += 1
            end
                # 
            # end
            prev_site = next_site
            next_site = rand(next_site_NPBC(next_site,n)) #rand([next_site+n,next_site-n,next_site+1,next_site-1])    
            # deleteat!(next_site, neighbours_set.==0)
            # pick_arr = Int[]

            
            # for obj in next_site
            #     if lat[obj].cluster_id != cl_id
            #         #deleteat!(next_site, findfirst(==(obj), next_site)) 
            #         push!(pick_arr, obj )
            #     end
            # end

            # if size(pick_arr)[1] == 0
            #     break
            # end
        
            
            # next_site = rand(pick_arr)

            lat[next_site].cluster_id = cl_id

        end



    end

    for i=1:len
        lat[i].cluster_id = 0
    end

   return lat , num_bonds
end


# initialise_lattice_sarw(3, 0.5)

# for i=1:9
#     if i%2==0
#         continue
#     end
#     println(i)
# end

# m = [ i%2 == 0 ? "YES" : i%3 == 0 ? "NOO" : 0 for i=1:10 ]

# println(m)

# hehe_arr = map( x -> x, [1,2,3] )


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


n = 100

# Random.seed!(4444444449999999)

# f = @time lattice_initialise_new(n,0.1);
f = @time initialise_lattice_sarw(n,40,100)
h = @time find_clusters_trial_test_allocs_fix(f[1]);
# my_heatmap_plot(h)
my_heatmap_plot_nocolor(h)

f[1].cluster_id


for i=1:20
    f = @time initialise_lattice_sarw(n,1,300);
    h = @time find_clusters_trial_test_allocs_fix(f);
    my_heatmap_plot(h)
    sleep(1)    
end


my_plot_new(h)
savefig("RW_attempt_1.png")



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

n = 100

f = @time lattice_initialise(n,0.5);
@time find_clusters_trial_test_allocs_fix(f);

heatmap_plot(f)

perc_prob(f)


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

dat_perc_prob = @time perc_prob_ens(100, 10^2)

plot(LinRange(0,1,size(dat_perc_prob)[1]), dat_perc_prob)

savefig("prob_perc_normal.png")


#############################################################################
#############################################################################


function perc_prob_ens_new(n,num_ens)
    n_p = 150 # 0-1 in how many steps
    prob_avg_ens = zeros(n_p)
    p = zeros(n_p)#LinRange(0,1,n_p)

    Threads.@threads for j=80:n_p #LinRange(0,1,n_p)
        for i=1:num_ens
            h = initialise_lattice_sarw(n,j,100)
            f = find_clusters_trial_test_allocs_fix(h[1])
            prob_avg_ens[j] += perc_prob(f)
            p[j] += h[2]
        end
        println(j)
    end

    return prob_avg_ens/num_ens, p/(n^2 * num_ens)
end

dat_perc_prob = @time perc_prob_ens_new(100, 5*10)

plot(dat_perc_prob[2], dat_perc_prob[1])


savefig("perc_500_prob_sarw_10ens_new.png")




############################################################################

