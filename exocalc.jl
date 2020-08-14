### A Pluto.jl notebook ###
# v0.11.5

using Markdown
using InteractiveUtils

# ╔═╡ c9ac27ee-dac0-11ea-2a8c-2d144b034a82
md"""
# Exoplanet Calculator 🪐
"""

# ╔═╡ b2286b26-dac2-11ea-1ce0-c7da562aa641
md"Given exoplanet and host star parameters from the literature, calculate derived values relevant for detection of the planet's atmosphere."

# ╔═╡ 19b35ef4-dac3-11ea-2d25-97e5482ff6a0
md"### Literature values"

# ╔═╡ 07db65d6-dd99-11ea-103b-33d0317af127
md"[Sources](https://exoplanetarchive.ipac.caltech.edu/overview/HAT-P-23%20b#legend)" 

# ╔═╡ 7cff6dfc-dd9f-11ea-1fdf-7b6aaa9435b4
md"### Results"

# ╔═╡ 0b6821a4-dac3-11ea-27d7-911521f0d3c0
md"### How parameters were calculated"

# ╔═╡ f8281da6-dd9f-11ea-1b6c-d32702215397
md"This first checks if there are missing input parameters and then calls the appropriate function to calculate them for each study. The resulting derived parameters are then calculated from them."

# ╔═╡ 49f75dea-dda0-11ea-1a85-bbdd4750b878
md"""These are the functions used to calculate each parameter based on the combination of inputs given. No `if` statements or default "`None`" keyword arguments needed thanks to Julia's multiple dispatch!"""

# ╔═╡ c5c5ea28-dd9e-11ea-1f89-5b1371831177
md"### Function for displaying the results"

# ╔═╡ 8e5811ae-dd9e-11ea-127e-b9812511492b
md"### Structure used to hold the possible input parameters used by a study"

# ╔═╡ a8df7ad0-dd9e-11ea-2a6a-f16683371016
md"### Structure used to hold a summary of all input and derived params"

# ╔═╡ 7db94ad6-dda1-11ea-2f33-1da144f1b7ad
md"Libraries for using things like physical constants and units."

# ╔═╡ 02bfa078-d62b-11ea-15df-d701431829b9
begin
	using Measurements, Unitful, UnitfulAstro, Parameters, Markdown
	using PhysicalConstants.CODATA2018: G, k_B, m_u, σ
	const amu, k = m_u, k_B
end;

# ╔═╡ 3f79c516-da77-11ea-1f6b-d3e7191a95d8
begin	
	# Star radius
	get_Rₛ(; Lₛ, Tₛ) = (Lₛ / (4.0π*σ*Tₛ^4))^(1//2)
	
	# Star-to-planet radius ratio
	get_RₚRₛ(; Rₚ, Rₛ) = Rₚ / Rₛ
	
	# Planet radius
	get_Rₚ(; RₚRₛ, Rₛ) =  RₚRₛ * Rₛ
	
	# Semi-major axis / Star density
	get_aRₛ(ρₛ::Unitful.Density, P::Unitful.Time) = ((G * P^2 * ρₛ)/(3.0π))^(1//3)
	get_aRₛ(a::Unitful.Length, Rₛ::Unitful.Length) = a / Rₛ
	
	# Star density
	get_ρₛ(P::Unitful.Time, aRₛ::Measurement) = (3.0π / (G * P^2)) * aRₛ^3
	get_ρₛ(Mₛ::Unitful.Mass, Rₛ::Unitful.Length) = Mₛ / ((4.0/3.0)π * Rₛ^3)

	# Star mass
	get_Mₛ(; ρₛ, Rₛ) = ρₛ * (4.0/3.0) * π * Rₛ^3.0
	get_Mₛ(Lₛ) = ((Lₛ / (u"Lsun"))^0.25)u"Msun" # MS-relation
	
	# Star luminosity
	get_Lₛ(; Tₛ, Rₛ) = 4.0π * Rₛ^2 * σ * Tₛ^4
	get_Tₛ(; Lₛ, Rₛ) = (L / (4.0π * Rₛ^2 * σ))^(1//4)
	
	#Planet mass
	get_Mₚ(; K, i, P, Mₛ) = (K/sin(i)) * (P / (2.0π*G))^(1//3) * Mₛ^(2//3)

	# Star surface gravity
	get_gₛ(; Mₛ, Rₛ) = G * Mₛ / Rₛ^2

	# Planet surface gravity
	get_gₚ(; Mₚ, RₚRₛ, Rₛ) = G * Mₚ / (RₚRₛ^2 * Rₛ^2)

	# Planet equilibrium temperature
	get_Tₚ(; Tₛ, aRₛ, α) = Tₛ * (1.0 - α)^(1//4) * (0.5/aRₛ)^(1//2)

	# Planet scale height
	get_H(; μ, Tₚ, gₚ) = k * Tₚ / (μ * gₚ)

	# Estimated signal from planet atmosphere
	get_ΔD(; H, RₚRₛ, Rₛ) = 2.0 * H * RₚRₛ/Rₛ
end;

# ╔═╡ 1110fdee-dde8-11ea-268c-8728e704c9c3
begin
	aRs = (4.26 ± 0.14)
	i = (87.4 ± 2.0)u"°"
	b = aRs * cos(i)
	b = uconvert(NoUnits, b)
	Measurements.value(b)
end

# ╔═╡ db28dbd2-db12-11ea-28e4-2b6cf30bd102
@with_kw_noshow struct Study @deftype Union{Nothing, Quantity}
	# Reference name (i.e. Ciceri et al. 2015)
	name::String = "Custom"
	
	# Orbital params
	RₚRₛ::Union{Measurement, Nothing} = nothing # Planet to star radius ratio
	aRₛ::Union{Measurement, Nothing} = nothing  # Semi-major axis to star radius ratio
	a::Union{Measurement, Nothing} = nothing    # Semi-major axis
	P = nothing                                 # Period
	K = nothing                                 # RV semi-amplitude
	i = nothing                                 # Inclination
	
	# Planet params
	μ = nothing                                 # Mean molecular weight
	α::Union{Measurement, Nothing} = nothing    # Albedo
	Tₚ = nothing                                # Equilibrium temperature
	ρₚ = nothing                                # Density
	Mₚ = nothing                                # Mass
	Rₚ = nothing                                # Radius
	gₚ = nothing                                # Surface gravity
	
	# Star params
	Tₛ = nothing                                # Effective temperature
	ρₛ = nothing                                # Density
	Mₛ = nothing                                # Mass
	Rₛ = nothing                                # Radius
	gₛ = nothing                                # Surface gravity
	Lₛ = nothing                                # Luminosity
end;

# ╔═╡ 17302d74-d63b-11ea-3de3-49f0df0554ca
# Input params from studies to explore
studies = [
	Study(
		name = "WASP-43/b: Weaver et al. (2020)",
		μ    = 2.0*amu,
		α    = 0.0 ± 0.0,
		i    = (1.433 ± 0.1)u"rad",
		P    = (0.813473978 ± 3.5e-8)u"d",
		Tₛ   = (4520 ± 120)u"K",
		Rₛ   = (0.667 ± 0.010)u"Rsun",
		aRₛ  = 4.872 ± 0.14,
		Mₛ   = (0.717 ± 0.025)u"Msun",
		Tₚ   = (1440.0 ± 40.0)u"K",
		Mₚ   = (2.052 ± 0.053)u"Mjup",
		K = (551.7 ± 4.7)u"m/s",
		Rₚ   =  (1.036 ± 0.012)u"Rjup",
	),
	Study(
		name = "HAT-P-23/b: Ciceri et al. (2015)",
		μ    = 2.0*amu,
		α    = 0.0 ± 0.0,
		K    = (368.5 ± 17.6)u"m/s",
		i    = (85.74 ± 0.95)u"°",
		P    = (1.21288287 ± 0.00000017)u"d",
		RₚRₛ = 0.11616 ± 0.00081,
		Tₛ   = (5885.0 ± 72.0)u"K",
		Rₛ   = (1.089 ± 0.028)u"Rsun",
		aRₛ  = (4.5459 ± 0.0919),
	),
	Study(
		name = "HAT-P-23/b: Stassun et al. (2017, GAIA DR1)",
		μ    = 2.0*amu,
		α    = 0.0 ± 0.0,
		K    = (346.0 ± 21)u"m/s", # latest RV data, from B17
		i    = (85.1 ± 1.5)u"°",  # latest RV data, from B17
		P    = (1.212880 ± 0.000002)u"d", # latest transit data: (S&R16)
		RₚRₛ = 0.1113 ± 0.0010, # latest transit data: (S&R16)
		Tₛ   = (5905.0 ± 80.0)u"K",
		ρₛ   = (0.92 ± 0.18)u"g/cm^3",
		Rₛ   = (0.960±0.200)u"Rsun",
	),
	Study(
		name = "HAT-P-23/b: GAIA DR2",
		μ    = 2.0*amu,
		α    = 0.0 ± 0.0,
		K    = (346.0 ± 21)u"m/s", # latest RV data, from B17
		i    = (85.1 ± 1.5)u"°",  # latest RV data, from B17
		P    = (1.2128867 ± 0.0000002)u"d", # latest transit data: (S&R16)
		RₚRₛ = 0.1113 ± 0.0010, # latest transit data: (S&R16)
		Tₛ   = (5734.0 ± 100.0)u"K",
		Rₛ   = (1.1858169 ±	0.0424133)u"Rsun",
		Lₛ   = (10.0^(0.13656067 ± 0.00864667))u"Lsun",
		#Mₚ   = (1.34 ± 0.59)u"Mjup", # DR1 mass, gives inconsistent ΔD=550ppm result
	),
	Study(
		name = "HAT-P-23/b: TICv8",
		μ    = 2.0*amu,
		α    = 0.0 ± 0.0,
		K    = (368.0 ± 21)u"m/s", # latest RV data, from B17
		i    = (85.1 ± 1.5)u"°",  # latest RV data, from B17
		P    = (1.2128867 ± 0.0000002)u"d", # latest transit data: (S&R16)
		RₚRₛ = 0.1113 ± 0.0010, # latest transit data: (S&R16)
		Tₛ   = (5918.230 ± 136.811)u"K",
		Rₛ   = (1.1517600 ± 0.0596583)u"Rsun",
		ρₛ = (0.99471000 ± 0.23240140)u"g/cm^3",
		Mₛ = (1.078000 ± 0.136618)u"Msun",
		#Lₛ = (10.0^(0.1661873 ± 0.0191600))u"Lsun",
		gₛ = (10.0^(4.3479600 ± 0.0819789))u"cm/s^2",
	),
];

# ╔═╡ bd752a9e-dd80-11ea-141c-779c5135d4d8
@with_kw_noshow struct Derived @deftype Union{Nothing, Quantity}
	name::String = "Custom"
	#Orbital params
	RₚRₛ::Union{Measurement, Nothing} = nothing
	inputs_RₚRₛ::Union{String, Nothing} = nothing
	P   = nothing
	inputs_P::Union{String, Nothing} = nothing
	aRₛ::Union{Measurement, Nothing} = nothing
	inputs_aRₛ::Union{String, Nothing} = nothing
	K   = nothing
	inputs_K::Union{String, Nothing} = nothing
	i   = nothing
	inputs_i::Union{String, Nothing} = nothing
	
	# Planet params
	μ   = nothing
	inputs_μ::Union{String, Nothing} = nothing
	α::Union{Measurement, Nothing} = nothing
	inputs_α::Union{String, Nothing} = nothing
	gₚ  = nothing
	inputs_gₚ::Union{String, Nothing} = nothing
	Mₚ  = nothing
	inputs_Mₚ::Union{String, Nothing} = nothing
	Rₚ  = nothing
	inputs_Rₚ::Union{String, Nothing} = nothing
	Tₚ  = nothing
	inputs_Tₚ::Union{String, Nothing} = nothing
	
	# Star Params
	ρₛ  = nothing
	inputs_ρₛ::Union{String, Nothing} = nothing
	gₛ  = nothing
	inputs_gₛ::Union{String, Nothing} = nothing
	Mₛ  = nothing
	inputs_Mₛ::Union{String, Nothing} = nothing
	Rₛ  = nothing
	inputs_Rₛ::Union{String, Nothing} = nothing
	Tₛ  = nothing
	inputs_Tₛ::Union{String, Nothing} = nothing
	Lₛ  = nothing
	inputs_Lₛ::Union{String, Nothing} = nothing
	
	# Signal
	H   = nothing
	ΔD  = nothing
end;

# ╔═╡ 3833772c-d63f-11ea-09b5-f36d68e512ea
begin
	results = []
	for st in studies
		# Rₛ, Tₛ, Lₛ
		if !isnothing(st.Rₛ)
			Rₛ = st.Rₛ
			inputs_Rₛ = "given"
			if any((!isnothing).([st.Tₛ, st.Lₛ]))
				if isnothing(st.Tₛ)
					Lₛ = st.Lₛ
					inputs_Lₛ = "given"
					Tₛ = get_Tₛ(Lₛ=Lₛ, Rₛ=Rₛ)
					inputs_Tₛ = "st.Lₛ, st.Rₛ"
				else
					Tₛ = st.Tₛ
					inputs_Tₛ = "given"
					Lₛ = get_Lₛ(Tₛ=Tₛ, Rₛ=Rₛ)
					inputs_Lₛ = "st.Tₛ, st.Rₛ"
				end
			else
				error("Rₛ was given. Lₛ or Tₛ must also be given.")
			end
		else
			if all((!isnothing).([st.Lₛ, st.Tₛ]))
				Lₛ, Tₛ = st.Lₛ, st.Tₛ
				inputs_Lₛ, inputs_Tₛ = "given", "given"
				Rₛ = get_Rₛ(Lₛ=Lₛ, Tₛ=Tₛ)
				inputs_Rₛ = "(st.Lₛ, st.Tₛ)"	
			else
				error("Rₛ was not given. Lₛ mad Tₛ must be given then.")
			end
		end
		
		# Mₛ
		if !isnothing(st.Mₛ)
			Mₛ = st.Mₛ
			inputs_Mₛ = "given"
		else
			Mₛ = get_Mₛ(Lₛ)
			inputs_Mₛ = "Lₛ $(inputs_Lₛ)"
		end
		
		# RₚRₛ and Rₚ
		if !isnothing(st.Rₚ)
			Rₚ = st.Rₚ
			inputs_Rₚ = "given"
			RₚRₛ = get_RₚRₛ(Rₚ=Rₚ, Rₛ=Rₛ)
			inputs_RₚRₛ = "st.Rₚ, $(inputs_Rₛ)"
		elseif !isnothing(st.RₚRₛ)
			RₚRₛ = st.RₚRₛ
			inputs_RₚRₛ = "given"
			Rₚ = get_Rₚ(RₚRₛ=RₚRₛ, Rₛ=Rₛ)
			inputs_Rₚ = "st.RₚRₛ, Rₛ $(inputs_Rₛ)"
		else
			error("RₚRₛ or Rₚ must be given.")
		end
		
		# P
		if !isnothing(st.P)
			P = st.P
			inputs_P = "given"
		else
			error("Please specify a period.")
		end
		
		# aRₛ, ρₛ
		if all((!isnothing).([st.aRₛ, st.ρₛ]))
			error("Inconsistent inputs. Only aRₛ or ρₛ can be given.")
		end			
		if !isnothing(st.aRₛ)
			aRₛ = st.aRₛ
			inputs_aRₛ = "given"
			ρₛ = get_ρₛ(P, aRₛ)
			inputs_ρₛ = "P, aRₛ"
		elseif !isnothing(st.a)
			a = st.a
			inputs_a = "given"
			aRₛ = get_aRₛ(a, Rₛ)
			inputs_aRₛ = "a, Rₛ"
			ρₛ = get_ρₛ(P, aRₛ)
			inputs_ρₛ = "P, aRₛ"
		elseif !isnothing(st.ρₛ)
			ρₛ = st.ρₛ
			inputs_ρₛ = "given"
			aRₛ = get_aRₛ(ρₛ, P)
			inputs_aRₛ = "ρₛ, P"
		elseif !isnothing(st.Mₛ)
			Mₛ = st.Mₛ
			inputs_Mₛ = "given"
			ρₛ = get_ρₛ(Mₛ, Rₛ)
			inputs_ρₛ = "Mₛ, Rₛ"
			aRₛ = get_aRₛ(ρₛ, P)
			inputs_aRₛ = "ρₛ, P"
		elseif !isnothing(st.Lₛ)
			Lₛ = st.Lₛ
			inputs_Lₛ = "given"
			Mₛ = get_Mₛ(Lₛ)
			inputs_Mₛ = "Lₛ"
			ρₛ = get_ρₛ(Mₛ, Rₛ)
			inputs_ρₛ = "Mₛ, Rₛ"
			aRₛ = get_aRₛ(ρₛ, P)
			inputs_aRₛ = "ρₛ, P"
		else
			error("Params not defined for ρₛ, P, a!")
		end
		
		# Calculate remaining params if not given
		if !isnothing(st.Mₚ)
			Mₚ = st.Mₚ
			inputs_Mₚ = "given"
		else
			Mₚ = get_Mₚ(K=st.K, i=st.i, P=P, Mₛ=Mₛ)
			inputs_Mₚ = "st.K, st.i, P=P, Mₛ"
		end
			
		Tₚ = (isnothing(st.Tₚ)) ? get_Tₚ(Tₛ=st.Tₛ, aRₛ=aRₛ, α=st.α) : st.Tₚ
		gₛ = (isnothing(st.gₛ)) ? get_gₛ(Mₛ=Mₛ, Rₛ=Rₛ) : st.gₛ
		gₚ = (isnothing(st.gₚ)) ? get_gₚ(Mₚ=Mₚ, RₚRₛ=RₚRₛ, Rₛ=Rₛ) : st.gₚ
		
		# Calculate depth
		H  = get_H(μ=st.μ, Tₚ=Tₚ, gₚ=gₚ)
		ΔD = get_ΔD(H=H, RₚRₛ=RₚRₛ, Rₛ=Rₛ)
		
		# Store summary
		summary = Derived(
			name = st.name,
			
			#Orbital params
			RₚRₛ = RₚRₛ,
			inputs_RₚRₛ = inputs_RₚRₛ, 
			P   = P,
			aRₛ = aRₛ,
			K   = st.K,
			i   = st.i,

			# Planet params
			μ   = st.μ,
			α   = st.α,
			gₚ  = gₚ,
			Mₚ  = Mₚ,
			Rₚ  = Rₚ,
			inputs_Rₚ = inputs_Rₚ,
			Tₚ  = Tₚ,
			H   = H,

			# Star Params
			ρₛ  = ρₛ,
			gₛ  = gₛ,
			Mₛ  = Mₛ,
			Rₛ  = Rₛ,
			Tₛ  = st.Tₛ,

			# Signal
			ΔD  = ΔD,
		)
		push!(results, summary)
	end
end;

# ╔═╡ 33fc58d0-dbd9-11ea-3c45-83f4b5a2a818
function display_summary(d::Derived)
	md"""
	###### **$(d.name):**
	**Orbital params** \
	RₚRₛ( $(d.inputs_RₚRₛ) ) = $(uconvert(NoUnits, d.RₚRₛ)) \
	P   = $(uconvert(u"d", d.P)) \
	aRₛ = $(uconvert(NoUnits, d.aRₛ)) \
	K   = $(uconvert(u"m/s", d.K)) \
	i   = $(uconvert(u"°", d.i)) 
	
	**Planet params** \
	μ   = $(uconvert(u"u", d.μ)) \
	α   = $(uconvert(NoUnits, d.α)) \
	gₚ  = $(uconvert(u"m/s^2", d.gₚ)) \
	Mₚ  = $(uconvert(u"Mjup", d.Mₚ)) \
	Rₚ( $(d.inputs_Rₚ) ) = $(uconvert(u"Rjup", d.Rₚ)) \
	Tₚ  = $(uconvert(u"K", d.Tₚ)) \
	H   = $(uconvert(u"km", d.H))

	**Star Params** \
	ρₛ  = $(uconvert(u"g/cm^3", d.ρₛ)) \
	log gₛ (cm/s²) = $(log10(ustrip(uconvert(u"cm/s^2", d.gₛ)))) \
	Mₛ  = $(uconvert(u"Msun", d.Mₛ)) \
	Rₛ  = $(uconvert(u"Rsun", d.Rₛ)) \
	Tₛ  = $(d.Tₛ)

	**Signal** \
	ΔD  = $(5 * uconvert(NoUnits, d.ΔD) * 1e6) ppm
	"""
end;

# ╔═╡ 4bfaf322-dbd9-11ea-0449-87d9aa07311f
display_summary.(results)

# ╔═╡ Cell order:
# ╟─c9ac27ee-dac0-11ea-2a8c-2d144b034a82
# ╟─b2286b26-dac2-11ea-1ce0-c7da562aa641
# ╟─19b35ef4-dac3-11ea-2d25-97e5482ff6a0
# ╟─07db65d6-dd99-11ea-103b-33d0317af127
# ╠═17302d74-d63b-11ea-3de3-49f0df0554ca
# ╟─7cff6dfc-dd9f-11ea-1fdf-7b6aaa9435b4
# ╠═4bfaf322-dbd9-11ea-0449-87d9aa07311f
# ╟─0b6821a4-dac3-11ea-27d7-911521f0d3c0
# ╟─f8281da6-dd9f-11ea-1b6c-d32702215397
# ╠═3833772c-d63f-11ea-09b5-f36d68e512ea
# ╟─49f75dea-dda0-11ea-1a85-bbdd4750b878
# ╠═3f79c516-da77-11ea-1f6b-d3e7191a95d8
# ╠═1110fdee-dde8-11ea-268c-8728e704c9c3
# ╟─c5c5ea28-dd9e-11ea-1f89-5b1371831177
# ╠═33fc58d0-dbd9-11ea-3c45-83f4b5a2a818
# ╟─8e5811ae-dd9e-11ea-127e-b9812511492b
# ╠═db28dbd2-db12-11ea-28e4-2b6cf30bd102
# ╟─a8df7ad0-dd9e-11ea-2a6a-f16683371016
# ╠═bd752a9e-dd80-11ea-141c-779c5135d4d8
# ╟─7db94ad6-dda1-11ea-2f33-1da144f1b7ad
# ╠═02bfa078-d62b-11ea-15df-d701431829b9
