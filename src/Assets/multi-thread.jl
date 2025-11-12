using Base.Threads

# Number of available threads
num_threads = nthreads()
println("Running on $num_threads threads")

# Example: Parallel computation of an array
n = 10
result = zeros(n)  # Preallocate array

Threads.@threads for i in 1:n
    delay
    result[i] = i^2  # Example computation
end

println("Result: ", result)