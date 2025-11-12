using Random, BenchmarkTools, SparseArrays

# Define Union-Find structure
struct UnionFind
    parent::Vector{Int}
    rank::Vector{Int}

    function UnionFind(n)
        new(collect(1:n), zeros(Int, n))
    end
end

function find(uf::UnionFind, x::Int)
    if uf.parent[x] != x
        uf.parent[x] = find(uf, uf.parent[x])  # Path compression
    end
    return uf.parent[x]
end

function union!(uf::UnionFind, x::Int, y::Int)
    rootX = find(uf, x)
    rootY = find(uf, y)

    if rootX != rootY
        if uf.rank[rootX] > uf.rank[rootY]
            uf.parent[rootY] = rootX
        elseif uf.rank[rootX] < uf.rank[rootY]
            uf.parent[rootX] = rootY
        else
            uf.parent[rootY] = rootX
            uf.rank[rootX] += 1
        end
    end
end

# Function to find islands using Union-Find
function find_islands(buses, lines)
    bus_list = collect(buses)
    bus_index = Dict(b => i for (i, b) in enumerate(bus_list))  # Map bus numbers to indices
    n = length(buses)
    
    uf = UnionFind(n)

    for (u, v) in lines
        union!(uf, bus_index[u], bus_index[v])
    end

    # Group buses by their root parent
    island_map = Dict{Int, Vector{Int}}()
    for bus in bus_list
        root = find(uf, bus_index[bus])
        if root ∉ keys(island_map)
            island_map[root] = []
        end
        push!(island_map[root], bus)
    end

    return collect(values(island_map))
end

# # Generate a random synthetic Polish 2736-bus system
# function generate_random_grid(num_buses, num_lines)
#     buses = Set(1:num_buses)
#     lines = Set{Tuple{Int,Int}}()

#     rng = MersenneTwister(42)  # Fixed seed for reproducibility
#     while length(lines) < num_lines
#         u, v = rand(rng, 1:num_buses), rand(rng, 1:num_buses)
#         if u != v
#             push!(lines, (min(u, v), max(u, v)))  # Avoid duplicate/self-loops
#         end
#     end

#     return buses, collect(lines)
# end


# Efficient function to extract load & line data for each island
function extract_island_data(islands::Vector{Vector{Int}}, Load::Matrix, Line::Matrix)
    island_loads = Dict{Int, Matrix}()
    island_lines = Dict{Int, Matrix}()

    # Create a fast lookup table for bus-to-island mapping
    bus_to_island = Dict{Int, Int}()
    for (i, island) in enumerate(islands)
        for bus in island
            bus_to_island[bus] = i
        end
    end

    # Convert Load and Line columns to sets for efficient filtering
    load_buses = Set(Load[:, 1])   # Set of all buses with loads
    line_buses = Set(vcat(Line[:, 1], Line[:, 2]))  # Set of all buses in lines

    # Extract loads for each island (using filtering)
    for (i, island) in enumerate(islands)
        island_buses = Set(island)  # Convert to Set for O(1) lookups
        island_loads[i] = Load[in.(Load[:, 1], Ref(island_buses)), :]
    end

    # Extract lines for each island (checking if both buses exist in the island)
    for (i, island) in enumerate(islands)
        island_buses = Set(island)
        island_lines[i] = Line[in.(Line[:, 1], Ref(island_buses)) .& in.(Line[:, 2], Ref(island_buses)), :]
    end

    return island_loads, island_lines
end





# Example Usage:
# islands = find_islands(buses, lines)  # Assume islands already computed
# Load and Line matrices are from MATPOWER
# Test on a large system (Polish-like: 2736 buses, ~3500 lines)

# num_buses = 2736
# num_lines = 3500
# buses, lines = generate_random_grid(num_buses, num_lines)

buses =  Set(B[:, 1])
lines =  Tuple{Int, Int}.(eachrow(L[:, 1:2])) 
trans =  Tuple{Int, Int}.(eachrow(T[:, 1:2])) 
lines =  vcat(lines, trans)

# Run and time the function
@time islands = find_islands(buses, lines)
println("Number of islands detected: ", length(islands))

# More precise benchmark
@benchmark find_islands($buses, $(collect(lines)))

@time island_loads, island_lines = extract_island_data(islands, Ld, L)

# Print results
for i in keys(island_loads)
    println("Island $i: Load size = ", size(island_loads[i]), ", Line size = ", size(island_lines[i]))
end



# function find_equivalent_rows(matrix)
#     row_tuples = [Tuple(row) for row in eachrow(matrix)]  # Convert rows to tuples
#     row_groups = countmap(row_tuples)  # Count occurrences of each unique row
    
#     # Find indices of duplicate rows
#     equivalent_rows = Dict()
#     for (i, row) in enumerate(row_tuples)
#         if row_groups[row] > 1  # Only consider duplicates
#             get!(equivalent_rows, row, []).append!(i)
#         end
#     end
    
#     return values(equivalent_rows) |> collect
# end