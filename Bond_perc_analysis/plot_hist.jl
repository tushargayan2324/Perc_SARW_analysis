function distribution_cluster_size(lat)
    cluster_size = [obj.cluster_id for obj in lat]
    #histogram(cluster_size)
    return countmap(cluster_size)
end


#my_lattice_unique = @time find_clusters_new_algo_not_condi(lattice_initialise(n,0.5))

my_data = distribution_cluster_size(l)

#g = distribution_cluster_size(find_clusters_new_algo_not_condi(lattice_initialise(n,0.6)))

n_s = collect(values(my_data))  

s = collect(keys(my_data)) 


S =  sum(s.^2 .* n_s) / ( sum(s .* n_s) )

N_p = 50 # 50 values of p

S_values = zeros(Float64, N_p)




function S_plotting()

    N_p = 60 # 50 values of p
    S_values = zeros(Float64, N_p)

    for i=1:N_p
        S = 0
        for j=1:5
            n = 200
            lat = find_clusters_trial_test(lattice_initialise(n,i/N_p))
            my_data = distribution_cluster_size(lat)
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

S_plotting()


n = 100

f = @time lattice_initialise(n,0.5);

g = @time lattice_initialise(n,0.5);

h = @time find_clusters_trial(f);

l = @time find_clusters_trial_test(g);


lat = @time find_clusters_trial_test(lattice_initialise(400,0.50))


# my_my_plot(lat)

cluster_size = [obj.cluster_id for obj in lat]

# println(cluster_size)

# countmap(cluster_size)

my_data = distribution_cluster_size(lat)
s = collect(values(my_data))

temp = sort(countmap(s))

# sort(temp)

n_s = collect(values(temp))

s_size = collect(keys(temp)) 

# histogram(temp)

# somedata = zip([1,2,3], [0.5,0.7,0.1])

scatter(log.(s_size), log.(n_s))
# savefig("freq_vs_clustersize.png")


s_xi = -0.3

scatter( (s_size.^2 .* n_s ), (s_size .* (s_xi) .* ones(size(s_size)[1])) )

# plot(somedata)

# histogram(s_size,n_s)
#s = collect(keys(my_data)) 

# S_values = zeros(Float64, N_p)


# histogram( log.(s) , collect(values(h)) , xlims=(0,200), ylims=(0,400))

# scatter( log.(collect(keys(h))/(n*n)) , collect(values(h)) )#, xlims=(0,200), ylims=(0,400))

# scatter((collect(values(h))), log.(collect(keys(h))/(n*n)), xlims=(0,200))


# histogram((s, log.(n_s)), )


