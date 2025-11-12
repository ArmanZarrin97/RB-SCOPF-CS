### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ e6027319-14fd-41a5-9852-4af244215506
import Pkg; Pkg.add("ProfileCanvas")

# ╔═╡ 0bf3ecd3-c5db-4b78-af95-c563f2e19a3d
using ProfileCanvas

# ╔═╡ 76cd0aa8-a9f9-4b0e-af60-ab6c9e3b9704
function a(b,c)
return b+c
end

# ╔═╡ 4cccab1c-64a1-462d-a1a4-c1e1c5586f55
@profview a(2,3)

# ╔═╡ Cell order:
# ╠═e6027319-14fd-41a5-9852-4af244215506
# ╠═0bf3ecd3-c5db-4b78-af95-c563f2e19a3d
# ╠═76cd0aa8-a9f9-4b0e-af60-ab6c9e3b9704
# ╠═4cccab1c-64a1-462d-a1a4-c1e1c5586f55
