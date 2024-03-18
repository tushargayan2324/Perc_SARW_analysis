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

   return lat 
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


n = 150

# Random.seed!(4444444449999999)

# f = @time lattice_initialise_new(n,0.1);
f = @time initialise_lattice_sarw(n,60,500);
h = @time find_clusters_trial_test_allocs_fix(f);
my_heatmap_plot(h)
my_heatmap_plot_nocolor(h)


for i=1:20
    f = @time initialise_lattice_sarw(n,1,300);
    h = @time find_clusters_trial_test_allocs_fix(f);
    my_heatmap_plot(h)
    sleep(1)    
end


my_plot_new(h)
savefig("RW_attempt_1.png")

scatter(0,0)
# plot!([1,2], [3,4])

plot!([1,2, 2], [3,3, 4])
# plot!([1,1], [3,4])



heatmap_plot(h)

my_heatmap_plot_nocolor(h)
my_heatmap_plot(h)


function my_heatmap_plot(f)
    n = trunc(Int64,sqrt(size(f)[1]))
    plt = scatter(0,0);

    g = reshape(f, (n,n))
    Mat = zeros(Float64, (n,n))
	for i=1:n
        for j=1:n
            if size(g[j,i].neighbours)[1] != 0
                Mat[i,j] = sin(11*g[j,i].cluster_id)
            else
                Mat[i,j] = (1 + (-1)^(i+j))/2
            end
        end
    end
    
    plt = heatmap(Mat,c=:devon,size=(1500,1500), legend=false) #c=:prism tab20 glasbey_hv_n256
    
    x_arr = []
    y_arr = []
    # c_arr = []
    for obj in f
        if size(obj.neighbours)[1] != 0    
            push!(x_arr, x_coord(obj.site_id,n))
            push!(y_arr, y_coord(obj.site_id,n))
            # push!(c_arr, obj.cluster_id)
        end
    end

    plt = scatter!(x_arr, y_arr,size=(1500,1500))


    for obj in f
        for nei in obj.neighbours
            x1 = x_coord(obj.site_id,n)
            y1 = y_coord(obj.site_id,n)    
            xs = [x1,x_coord(nei,n)]
            ys = [y1,y_coord(nei,n)]
            plt = plot!(xs,ys,color=:black, linewidth=5,size=(1500,1500), legend=false)
        end    
    end
    


    Plots.savefig("hehe.png")
end


function my_heatmap_plot_nocolor(f)
    n = trunc(Int64,sqrt(size(f)[1]))
    plt = scatter(0,0);

    g = reshape(f, (n,n))
    Mat = zeros(Float64, (n,n))
	# for i=1:n
    #     for j=1:n
    #         if size(g[j,i].neighbours)[1] != 0
    #             Mat[i,j] = sin(11*g[j,i].cluster_id)
    #         else
    #             Mat[i,j] = (1 + (-1)^(i+j))/2
    #         end
    #     end
    # end
    
    # plt = heatmap(Mat,c=:devon,size=(1500,1500), legend=false) #c=:prism tab20 glasbey_hv_n256
    
    for obj in f
        for nei in obj.neighbours
            x1 = x_coord(obj.site_id,n)
            y1 = y_coord(obj.site_id,n)    
            xs = [x1,x_coord(nei,n)]
            ys = [y1,y_coord(nei,n)]
            plt = plot!(xs,ys,color=:black, linewidth=5,size=(1500,1500), legend=false)
        end    
    end
    


    Plots.savefig("hehe.png")
end


function my_plot_new(f)
	plt = scatter(0,0);
    n = trunc(Int64,sqrt(size(f)[1]))

    x_arr = []
    y_arr = []
    # c_arr = []
    for obj in f
        if size(obj.neighbours)[1] != 0    
            push!(x_arr, x_coord(obj.site_id,n))
            push!(y_arr, y_coord(obj.site_id,n))
            # push!(c_arr, obj.cluster_id)
        end
    end

    plt = scatter!(x_arr, y_arr,size=(1500,1500))

    for obj in f
        for nei in obj.neighbours
            x1 = x_coord(obj.site_id,n)
            y1 = y_coord(obj.site_id,n)    
            xs = [x1,x_coord(nei,n)]
            ys = [y1,y_coord(nei,n)]
            plt = plot!(xs,ys,color=:black, linewidth=3,size=(1500,1500), legend=false)
        end    
    end

    return plt
end















my_plot(h)

cluster_size_freq(h)



my_my_plot(h)





# println(f)

# push!(f[2].neighbours,3)

# My_set = Set{Int64}([1,2,3])

# push!(My_set, 1)

# for obj in f
#     if rand() < 0.1
#         # println(1)
#         start_site_nei = lat[i].neighbours,rand([i+n,i-n,i+1,i-1])
#         #push!(lat[i].neighbours,rand([i-n,i+1,i-1]))

#         for j=1:5#rand([3,4,5])
#             println(j)
#             next_site = rand([i+n,i-n,i+1,i-1])
#             if !(next_site in obj.neighbours)
#                 push!(lat[i].neighbours,next_site)
#             end
#         end
#     end
# end


# for i=1:rand([3,4,5])
#     println(i)
# end


#S is the cluster size, hence the following is plot of avg cluster size over all range of p values
function S_plotting()

    N_p = 100 # 50 values of p
    S_values = zeros(Float64, N_p)

    for i=1:N_p
    # for i=1:N_p
        S = 0
        for j=1:5
            n = 200
            lat = find_clusters_trial_test_allocs_fix(lattice_initialise_new(n,i/N_p))
            
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


## Create function for plotting the probability of finding percolating cluster in some
## particular direction. Either just check cluster_ids of sites or create a new function 
## to just find if its percolating or not