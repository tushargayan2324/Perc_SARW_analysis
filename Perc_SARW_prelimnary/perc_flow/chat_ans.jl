using JuMP
using GLPK

function equivalent_resistance(network::Vector{Tuple{Int, Int, Float64}}, nodes::Tuple{Int, Int})
    num_nodes = maximum(maximum.(network[:, 1:2]))
    source, target = nodes
    
    # Create a JuMP model
    model = Model(GLPK.Optimizer)
    
    # Define variables: current through each resistor
    @variable(model, I[1:size(network, 1)] >= 0)
    
    # Define variable: total resistance
    @variable(model, Req >= 0)
    
    # Kirchhoff's Current Law (KCL) for each node except the target and source
    for i in 1:num_nodes
        if i != source && i != target
            @constraint(model, sum(I[j] for j in findall(x->x==i, network[:, 1])) ==
                                sum(I[j] for j in findall(x->x==i, network[:, 2])))
        end
    end
    
    # KCL for the source node
    @constraint(model, sum(I[j] for j in findall(x->x==source, network[:, 1])) ==
                        sum(I[j] for j in findall(x->x==source, network[:, 2])) - 1)
    
    # KCL for the target node
    @constraint(model, sum(I[j] for j in findall(x->x==target, network[:, 1])) ==
                        sum(I[j] for j in findall(x->x==target, network[:, 2])) + 1)
    
    # Ohm's Law: Req = V/I
    @constraint(model, Req == sum(I[j] * network[j, 3] for j in 1:size(network, 1)))
    
    # Objective: minimize total resistance
    @objective(model, Min, Req)
    
    # Solve the optimization problem
    optimize!(model)
    
    # Return the optimal equivalent resistance
    return value(Req)
end

# Example usage:
# Define the network as a list of tuples: (node1, node2, resistance)
network = [(1, 2, 5.0), (2, 3, 3.0), (3, 4, 4.0), (4, 1, 6.0)]
# Define the nodes between which you want to find the equivalent resistance
nodes = (1, 3)

# Call the function to find the equivalent resistance
eq_resistance = equivalent_resistance(network, nodes)
println("Equivalent Resistance between nodes $nodes: $eq_resistance")
