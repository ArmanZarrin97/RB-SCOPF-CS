import Pkg; Pkg.activate("ST-SCOPF");
Pkg.instantiate()
# Pkg.update()
# Pkg.resolve()
using KLU,JuMP,BSON, SparseArrays, Ipopt,OdsIO, DataFrames, Plots,LinearAlgebra,JLD2,PlotThemes,Random,Statistics,CPLEX,HiGHS,HSL_jll, ProgressLogging,BenchmarkTools,Profile,Base.Threads,SplitApplyCombine,CairoMakie #, ProfileCanvas, ProfileView,Cthulhu
#********************************************************************************************************************
# JuMP: Modeling Optimization Problem
#Ipopt: Optimization solver paccage
# OdsIO: Input & output for .ods files
# DataFrames: Data structure
# Plots: Plotting package
# LinearAlgebra: Linear Algebra package
# JLD2: Data saving package
# PlotThemes: Plotting themes
# Random: Random number generator
# HSL: ma97 solver package
# Statistics: Statistical package
#********************************************************************************************************************


# using ProgressLogging
# ProgressLogging.@progress name="My loop" for i in 1:100
# 	sleep(0.1)
# end

# function profile_test(n)
#     for i = 1:n
#         A = randn(100,100,20)
#         m = maximum(A)
#         Am = mapslices(sum, A; dims=2)
#         B = A[:,:,5]
#         Bsort = mapslices(sort, B; dims=1)
#         b = rand(100)
#         C = B.*b
#     end
# end

# @profview profile_test(1)
# ProfileCanvas.@profview profile_test(1)

            # @descend Flow_Expr(LinMethod,Lininit,nTP,raw,Vars,'L';SenCoef=SensModel)
            # ProfileView.@profview Flow_Expr(LinMethod,Lininit,nTP,raw,Vars,'L';SenCoef=SensModel)


#             function LIBHSL_isfunctional()
#             @ccall libhsl.LIBHSL_isfunctional()::Bool
#             end
#             bool = LIBHSL_isfunctional()


#             using HSL
# using LinearAlgebra
# using SparseArrays
# A = sprand(100, 100, 0.1)
# b = rand(100)
# x = ma97_solve(A, b)


# using JuMP
# import Ipopt, HSL_jll
# # An optimization problem
# model = Model(Ipopt.Optimizer)
# @variable(model, x)
# @objective(model, Min, (x - 2) ^ 2)

# # Load the HSL solvers
# set_attribute(model, "hsllib", HSL_jll.libhsl_path)
# set_attribute(model, "linear_solver", "ma97")
# optimize!(model)