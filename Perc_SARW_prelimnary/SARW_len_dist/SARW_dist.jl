using Plots

function length_SARW(max_itr) #Generating and getting back the length of generated SARW
    str_site = [0,0]
    
    next_site = copy(str_site)

    sites_arr = [[0,0] for i=1:max_itr]

    neigh_add = [ [0,1], [0,-1], [1,0], [-1,0]  ]

    sites_arr[1] = str_site

    for i=2:max_itr
        prev_site = next_site
        next_site_arr = [prev_site + obj for obj in neigh_add] #rand([ [0,1], [0,-1], [1,0], [-1,0]  ])
        # println(next_site_arr)
        choose_site_arr = []

        for s in next_site_arr
            if !(s in sites_arr)
                push!(choose_site_arr, s)
            end
        end

        if size(choose_site_arr)[1] == 0
            return i
            break
        end

        next_site = rand(choose_site_arr)
        sites_arr[i] = next_site
    end


    return max_itr
end





function dist_len_SARW(Num_ens, max_itr)
    m = zeros(Int, Num_ens)
    
    Threads.@threads for i=1:Num_ens
        m[i] = length_SARW(max_itr)
    end

    return m
end


function rg_length_SARW(max_itr) # Finding Radius of Gyration for a given length of SARW
    str_site = [0,0]
    
    next_site = copy(str_site)

    sites_arr = [[0,0] for i=1:max_itr]

    neigh_add = [ [0,1], [0,-1], [1,0], [-1,0]  ]

    sites_arr[1] = str_site

    for i=2:max_itr
        prev_site = next_site
        next_site_arr = [prev_site + obj for obj in neigh_add] #rand([ [0,1], [0,-1], [1,0], [-1,0]  ])
        # println(next_site_arr)
        choose_site_arr = []

        for s in next_site_arr
            if !(s in sites_arr)
                push!(choose_site_arr, s)
            end
        end

        if size(choose_site_arr)[1] == 0
            # return i, 1/(i).*sum(sites_arr .- [sum(sites_arr)/i for i=1:max_itr])
            return i, ( 1/(i)*sum(sum(sites_arr - [sum(sites_arr)/i for i=1:max_itr]).^2) )^0.5
            break
        end

        next_site = rand(choose_site_arr)
        sites_arr[i] = next_site
    end


    return max_itr, 1/(max_itr).*sum(sum(sites_arr - [sum(sites_arr)/max_itr for i=1:max_itr]).^2)
end

rg_length_SARW(500)

function dist_len_rad_gyr_SARW(Num_ens, max_itr, len)
    m = []#zeros(Int, Num_ens)
    
    Threads.@threads for i=1:Num_ens
        temp = rg_length_SARW(max_itr)
        if temp[1] == len
            push!(m, temp[2])
        end
    end

    return m
end


hist_dat_20 = dist_len_rad_gyr_SARW(10^4, 500, 20)

hist_dat_40 = dist_len_rad_gyr_SARW(10^4, 500, 40)

hist_dat_50 = dist_len_rad_gyr_SARW(10^4, 500, 50)


hist_dat = dist_len_SARW(3*10^4, 10^4) 




histogram(hist_dat_50)


# [[0,1] , [1,0]] .- [[1,0], [0,2]]


# str_site = [0,0]

# str_site += [0,1]


# sites_arr = [[0,0] for i=1:10]#Vector{Vector{Int64}}[]

# sites_arr[1] = (0,1)#str_site

# sites_arr[1]
# sites_arr[7] = [2,9]

# m = [2,9]

# findfirst((2,9),sites_arr)
# sites_arr .== (2,9)

# findfirst(==(m), sites_arr)

# m in sites_arr


# hehehe = ([0,1])