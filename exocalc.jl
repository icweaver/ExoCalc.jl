### A Pluto.jl notebook ###
# v0.11.8

using Markdown
using InteractiveUtils

# ╔═╡ 02bfa078-d62b-11ea-15df-d701431829b9
begin
	using Measurements, Unitful, UnitfulAstro, Parameters, Markdown
	using PhysicalConstants.CODATA2018: G, k_B, m_u, σ
	const amu, k = m_u, k_B
end;

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

# ╔═╡ 38a61304-e0fe-11ea-14b2-17d9b9e13c7b
md"### Calculate parameters"

# ╔═╡ 0b6821a4-dac3-11ea-27d7-911521f0d3c0
md"### How parameters were calculated"

# ╔═╡ f8281da6-dd9f-11ea-1b6c-d32702215397
md"""
The calculator first checks if there are missing input parameters and then calls the appropriate function to calculate them for each study. The resulting derived parameters are then calculated from them.

Everything was done with a "star first" approach, meaning that all stellar parameters were determined first, and then the planet parameters were determined self-consistently from that. If conflicting parameters are given, the calculator will try to give priority to using direct observables to perform calculations and error otherwise. 

For transparency, the inputs used for each calculation are show in parenthesis next to each parameter.
"""

# ╔═╡ 49f75dea-dda0-11ea-1a85-bbdd4750b878
md"""These are the functions used to calculate each parameter based on the combination of inputs given. No `if` statements or default "`None`" keyword arguments needed thanks to Julia's multiple dispatch!"""

# ╔═╡ 3f79c516-da77-11ea-1f6b-d3e7191a95d8
begin
	# Star radius
	get_Rₛ(; Lₛ, Tₛ) = (Lₛ / (4.0π*σ*Tₛ^4))^(1//2)

	# Star-to-planet radius ratio
	get_RₚRₛ(; Rₚ, Rₛ) = Rₚ / Rₛ

	# Semi-major axis / Star density
	get_aRₛ(ρₛ::Unitful.Density, P::Unitful.Time) = ((G * P^2 * ρₛ)/(3.0π))^(1//3)
	get_aRₛ(a::Unitful.Length, Rₛ::Unitful.Length) = a / Rₛ
	
	# Semi-major axis
	get_a(; aRₛ, Rₛ) = aRₛ * Rₛ
	
	# Impact parameter
	get_b(; i, aRₛ) = aRₛ * cos(i)

	# Star density
	get_ρₛ(P::Unitful.Time, aRₛ::Measurement) = (3.0π / (G * P^2)) * aRₛ^3
	get_ρₛ(Mₛ::Unitful.Mass, Rₛ::Unitful.Length) = Mₛ / ((4.0/3.0)π * Rₛ^3)

	# Star mass
	get_Mₛ(; ρₛ, Rₛ) = ρₛ * (4.0/3.0) * π * Rₛ^3.0

	# Star luminosity
	get_Lₛ(; Tₛ, Rₛ) = 4.0π * Rₛ^2 * σ * Tₛ^4
	
	# Star temperature 
	get_Tₛ(; Lₛ, Rₛ) = (L / (4.0π * Rₛ^2 * σ))^(1//4)

	# Planet mass
	get_Mₚ(; K, i, P, Mₛ) = (K/sin(i)) * (P / (2.0π*G))^(1//3) * Mₛ^(2//3)
	
	# Planet radius
	get_Rₚ(; RₚRₛ, Rₛ) =  RₚRₛ * Rₛ
	
	# Planet density
	get_ρₚ(; Mₚ, Rₚ) = Mₚ / ((4.0/3.0)π * Rₚ^3)

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

# ╔═╡ c5c5ea28-dd9e-11ea-1f89-5b1371831177
md"### Function to display the results"

# ╔═╡ 8e5811ae-dd9e-11ea-127e-b9812511492b
md"### Structure used to hold the possible input parameters used by a study"

# ╔═╡ db28dbd2-db12-11ea-28e4-2b6cf30bd102
@with_kw_noshow struct Study @deftype Union{Nothing, Quantity, Measurement}
	# Reference name (e.g. Ciceri et al. 2015)
	name::String = "Custom"
	
	# Star params
	Tₛ = nothing
	ρₛ = nothing
	Mₛ = nothing
	Rₛ = nothing
	gₛ = nothing
	Lₛ = nothing

	# Orbital params
	RₚRₛ = nothing
	aRₛ = nothing
	a = nothing
	b = nothing
	P = nothing
	K = nothing
	i = nothing

	# Planet params
	μ = nothing
	α = nothing
	Tₚ = nothing
	ρₚ = nothing
	Mₚ = nothing
	Rₚ = nothing
	gₚ = nothing
	N_scales::Float64 = 5.0 # Number of scale heights (default 5)
end;

# ╔═╡ 17302d74-d63b-11ea-3de3-49f0df0554ca
# Input params from studies to explore
studies = [
	Study(
		name = "WASP-43/b: Weaver et al. (2020)",
		μ	 = 2.0*amu,
		α	 = 0.0 ± 0.0,
		i	 = (1.433 ± 0.1)u"rad",
		P	 = (0.813473978 ± 3.5e-8)u"d",
		Tₛ	 = (4520 ± 120)u"K",
		Rₛ	 = (0.667 ± 0.010)u"Rsun",
		aRₛ  = 4.872 ± 0.14,
		Mₛ	 = (0.717 ± 0.025)u"Msun",
		Tₚ	 = (1440.0 ± 40.0)u"K",
		Mₚ	 = (2.052 ± 0.053)u"Mjup",
		K = (551.7 ± 4.7)u"m/s",
		Rₚ	 =	(1.036 ± 0.012)u"Rjup",
	),
	Study(
		name = "HAT-P-23/b: Ciceri et al. (2015)",
		μ	 = 2.0*amu,
		α	 = 0.0 ± 0.0,
		K	 = (368.5 ± 17.6)u"m/s",
		i	 = (85.74 ± 0.95)u"°",
		P	 = (1.21288287 ± 0.00000017)u"d",
		RₚRₛ = 0.11616 ± 0.00081,
		Tₛ	 = (5885.0 ± 72.0)u"K",
		Rₛ	 = (1.089 ± 0.028)u"Rsun",
		aRₛ  = (4.5459 ± 0.0919),
	),
	Study(
		name = "HAT-P-23/b: Sada & Ramón-Fox (2016)",
		μ	 = 2.0*amu,
		α	 = 0.0 ± 0.0,
		K	 = (346.0 ± 21.0)u"m/s", # latest RV data, from B17
		i	 = (85.1 ± 1.5)u"°",  # latest RV data, from B17
		P	 = (1.212880 ± 0.000002)u"d", # latest transit data: (S&R16)
		RₚRₛ = 0.1113 ± 0.0010, # latest transit data: (S&R16)
		Tₛ	 = (5905.0 ± 80.0)u"K",
		Rₛ	 = (0.960±0.200)u"Rsun",
		aRₛ  = 4.26 ± 0.14,
	),
	Study(
		name = "HAT-P-23/b: Stassun et al. (2017, GAIA DR1)",
		μ	 = 2.0*amu,
		α	 = 0.0 ± 0.0,
		K	 = (368.5 ± 17.6)u"m/s",
		i	 = (85.1 ± 1.5)u"°",  # latest RV data, from B17
		P	 = (1.212880 ± 0.000002)u"d",
		RₚRₛ = 0.1113 ± 0.0010,
		Tₛ	 = (5905.0 ± 80.0)u"K",
		ρₛ	 = (0.92 ± 0.18)u"g/cm^3",
		Rₛ	 = (0.960±0.200)u"Rsun",
	),
	Study(
		name = "HAT-P-23/b: GAIA DR2",
		μ	 = 2.0*amu,
		α	 = 0.0 ± 0.0,
		K	 = (346.0 ± 21)u"m/s", # latest RV data, from B17
		i	 = (85.1 ± 1.5)u"°",  # latest RV data, from B17
		P	 = (1.2128867 ± 0.0000002)u"d", # latest transit data: (S&R16)
		RₚRₛ = 0.1113 ± 0.0010, # latest transit data: (S&R16)
		Tₛ	 = (5734.0 ± 100.0)u"K",
		Rₛ	 = (1.1858169 ± 0.0424133)u"Rsun",
		Lₛ	 = (10.0^(0.13656067 ± 0.00864667))u"Lsun",
		#Mₚ	 = (1.34 ± 0.59)u"Mjup", # DR1 mass, gives inconsistent ΔD=550ppm result
		aRₛ  = (4.26 ± 0.14), # latest transit data: (S&R16)
	),
	Study(
		name = "HAT-P-23/b: TICv8",
		μ	 = 2.0*amu,
		α	 = 0.0 ± 0.0,
		K	 = (346.0 ± 21.0)u"m/s", # latest RV data, from B17
		i	 = (85.1 ± 1.5)u"°",  # latest RV data, from B17
		P	 = (1.2128867 ± 0.0000002)u"d", # latest transit data: (S&R16)
		RₚRₛ = 0.1113 ± 0.0010, # latest transit data: (S&R16)
		Tₛ	 = (5918.230 ± 136.811)u"K",
		Rₛ	 = (1.1517600 ± 0.0596583)u"Rsun",
		ρₛ = (0.99471000 ± 0.23240140)u"g/cm^3",
		Mₛ = (1.078000 ± 0.136618)u"Msun",
		Lₛ = (10.0^(0.1661873 ± 0.0191600))u"Lsun",
		gₛ = (10.0^(4.3479600 ± 0.0819789))u"cm/s^2",
	),
];

# ╔═╡ c01eb856-e0f9-11ea-01d5-07593189ce46
function calculate_params(st::Study)
	# Rₛ, Tₛ, Lₛ
	if !isnothing(st.Rₛ)
		Rₛ = st.Rₛ
		inputs_Rₛ = (Rₛ=Rₛ,)
		if any((!isnothing).([st.Tₛ, st.Lₛ]))
			if all((!isnothing).([st.Tₛ, st.Lₛ]))
				Lₛ = st.Lₛ
				inputs_Lₛ = (Lₛ=Lₛ,)
				Tₛ = st.Tₛ
				inputs_Tₛ = (Tₛ=Tₛ,)
			elseif isnothing(st.Tₛ)
				Lₛ = st.Lₛ
				inputs_Lₛ = (Lₛ=Lₛ,)
				inputs_Tₛ = (Lₛ=Lₛ, Rₛ=Rₛ)
				Tₛ = get_Tₛ(; inputs_Tₛ...)
			else
				Tₛ = st.Tₛ
				inputs_Tₛ = (Tₛ=Tₛ,)
				inputs_Lₛ = (Tₛ=Tₛ, Rₛ=Rₛ)
				Lₛ = get_Lₛ(; inputs_Lₛ...)
			end
		else
			error("Rₛ was given. Lₛ or Tₛ must also be given.")
		end
	else
		if all((!isnothing).([st.Lₛ, st.Tₛ]))
			Lₛ, Tₛ = st.Lₛ, st.Tₛ
			inputs_Lₛ, inputs_Tₛ = (Lₛ=Lₛ,), (Tₛ=Tₛ,)
			inputs_Rₛ = (Lₛ=Lₛ, Tₛ=Tₛ)
			Rₛ = get_Rₛ(; inputs_Rₛ...)
		else
			error("Rₛ was not given. Lₛ and Tₛ must be given then.")
		end
	end

	# RₚRₛ and Rₚ
	if !isnothing(st.Rₚ)
		Rₚ = st.Rₚ
		inputs_Rₚ = (Rₚ=Rₚ,)
		inputs_RₚRₛ = (Rₚ=Rₚ, Rₛ=Rₛ)
		RₚRₛ = get_RₚRₛ(; inputs_RₚRₛ...)
	elseif !isnothing(st.RₚRₛ)
		RₚRₛ = st.RₚRₛ
		inputs_RₚRₛ = (RₚRₛ=RₚRₛ,)
		inputs_Rₚ = (RₚRₛ=RₚRₛ, Rₛ=Rₛ)
		Rₚ = get_Rₚ(; inputs_Rₚ...)
	else
		error("Please specify either RₚRₛ or Rₚ.")
	end

	# P
	if !isnothing(st.P)
		P = st.P
		inputs_P = (P=P,)
	else
		error("Please specify a period.")
	end

	# ρₛ, aRₛ, a
	if all((!isnothing).([st.aRₛ, st.ρₛ]))
		error("Conflicting inputs. Only aRₛ or ρₛ can be given.")
	end
	if all((!isnothing).([st.aRₛ, st.b]))
		error("Conflicting inputs. Only aRₛ or b can be given.")
	end
	if !isnothing(st.ρₛ)
		ρₛ = st.ρₛ
		inputs_ρₛ = (ρₛ=ρₛ,)
		inputs_aRₛ = (ρₛ=ρₛ, P=P)
		aRₛ = get_aRₛ(inputs_aRₛ...)
		inputs_a = (aRₛ=aRₛ, Rₛ=Rₛ)
		a = get_a(; inputs_a...)
	elseif !isnothing(st.aRₛ)
		aRₛ = st.aRₛ
		inputs_aRₛ = (aRₛ=aRₛ,)
		inputs_a = (aRₛ=aRₛ, Rₛ=Rₛ)
		a = get_a(; inputs_a...)
		inputs_ρₛ = (P=P, aRₛ=aRₛ)
		ρₛ = get_ρₛ(inputs_ρₛ...)
	elseif !isnothing(st.a)
		a = st.a
		inputs_a = (a=a,)
		inputs_aRₛ = (a=a, Rₛ=Rₛ)
		aRₛ = get_aRₛ(inputs_aRₛ...)
		inputs_ρₛ = (P=P, aRₛ=aRₛ)
		ρₛ = get_ρₛ(inputs_ρₛ...)
	else
		error("ρₛ or (aRₛ or a) must be given for $(st.name)")
	end

	# Mₛ
	if !isnothing(st.Mₛ)
		Mₛ = st.Mₛ
		inputs_Mₛ = (Mₛ=Mₛ,)
	else
		inputs_Mₛ = (ρₛ=ρₛ, Rₛ=Rₛ)
		Mₛ = get_Mₛ(; inputs_Mₛ...)
	end

	# Calculate remaining params if not given/calculated
	if isnothing(st.i)
		error("Must provide inclination (i).")
	else
		i = st.i
		inputs_i = (i=i,)
	end
	if isnothing(st.K)
		error("Must provide RV semi-amplitude (K).")
	else
		K = st.K
		inputs_K = (K=K,)
	end
	if isnothing(st.α)
		error("Must provide albedo (α).")
	else
		α = st.α
		inputs_α = (α=α,)
	end
	if !isnothing(st.b)
		b = st.b
		inputs_b = (b=b,)
	else
		inputs_b = (aRₛ=aRₛ, i=i)
		b = get_b(; inputs_b...)
	end
	if !isnothing(st.Mₚ)
		Mₚ = st.Mₚ
		inputs_Mₚ = (Mₚ=Mₚ,)
	else
		inputs_Mₚ = (K=K, i=i, P=P, Mₛ=Mₛ)
		Mₚ = get_Mₚ(; inputs_Mₚ...)
	end
	if !isnothing(st.Tₚ)
		Tₚ = st.Tₚ
		inputs_Tₚ = (Tₚ=Tₚ,)
	else
		inputs_Tₚ = (Tₛ=Tₛ, aRₛ=aRₛ, α=α)
		Tₚ = get_Tₚ(; inputs_Tₚ...)
	end
	if !isnothing(st.gₛ)
		gₛ = st.gₛ
		inputs_gₛ = (gₛ=gₛ,)
	else
		inputs_gₛ = (Mₛ=Mₛ, Rₛ=Rₛ)
		gₛ = get_gₛ(; inputs_gₛ...)
	end
	if !isnothing(st.gₚ)
		gₚ = st.gₚ
		inputs_gₚ = (gₚ=gₚ,)
	else
		inputs_gₚ = (Mₚ=Mₚ, RₚRₛ=RₚRₛ, Rₛ=Rₛ)
		gₚ = get_gₚ(; inputs_gₚ...)
	end
	if !isnothing(st.ρₚ)
		ρₚ = st.ρₚ
	else
		inputs_ρₚ = (Mₚ=Mₚ, Rₚ=Rₚ)
		ρₚ = get_ρₚ(; inputs_ρₚ...)
	end

	# Calculate signal
	if isnothing(st.μ)
		error("Must provide mean molecula weight (μ).")
	else
		μ = st.μ
		inputs_μ = (μ=μ,)
	end
	inputs_H = (μ=μ, Tₚ=Tₚ, gₚ=gₚ)
	H  = get_H(; inputs_H...)
	ΔD = get_ΔD(H=H, RₚRₛ=RₚRₛ, Rₛ=Rₛ)
	
	# Store results
	params = (
		# Star Params
		ρₛ	 = ρₛ,
		gₛ	 = gₛ,
		Mₛ	 = Mₛ,
		Rₛ	 = Rₛ,
		Tₛ	 = Tₛ,
		Lₛ   = Lₛ,

		#Orbital params
		RₚRₛ = RₚRₛ,
		P	 = P,
		aRₛ  = aRₛ,
		a    = a,
		K	 = K,
		i	 = i,
		b    = b,

		# Planet params
		μ	 = μ,
		α	 = α,
		gₚ	 = gₚ,
		Mₚ	 = Mₚ,
		Rₚ	 = Rₚ,
		ρₚ   = ρₚ,
		Tₚ	 = Tₚ,
		H	 = H,

		# Signal
		N_scales = st.N_scales,
		ΔD	= ΔD,
	)
	
	params_inputs = (
		# Star Params
		inputs_ρₛ   = inputs_ρₛ,
		inputs_gₛ   = inputs_gₛ,
		inputs_Mₛ   = inputs_Mₛ,
		inputs_Rₛ   = inputs_Rₛ,
		inputs_Tₛ   = inputs_Tₛ,
		inputs_Lₛ   = inputs_Lₛ,

		#Orbital params
		inputs_RₚRₛ = inputs_RₚRₛ,
		inputs_P    = inputs_P,
		inputs_aRₛ  = inputs_aRₛ,
		inputs_a    = inputs_a,
		inputs_K    = inputs_K,
		inputs_i    = inputs_i,
		inputs_b    = inputs_b,

		# Planet params
		inputs_μ    = inputs_μ,
		inputs_α    = inputs_α,
		inputs_gₚ   = inputs_gₚ,
		inputs_Mₚ   = inputs_Mₚ,
		inputs_Rₚ   = inputs_Rₚ,
		inputs_ρₚ   = inputs_ρₚ,
		inputs_Tₚ   = inputs_Tₚ,
		inputs_H    = inputs_H,
	)
	
	return params, params_inputs;
end

# ╔═╡ a8df7ad0-dd9e-11ea-2a6a-f16683371016
md"### Structure to hold a summary of all parameters"

# ╔═╡ bd752a9e-dd80-11ea-141c-779c5135d4d8
@with_kw_noshow struct Derived @deftype Quantity
	# Reference name (e.g. Ciceri et al. 2015)
	name::String = "Custom"
	
	# Star Params
	ρₛ
	gₛ
	Mₛ
	Rₛ
	Tₛ
	Lₛ

	#Orbital params
	RₚRₛ::Measurement
	P
	aRₛ::Measurement
	a
	b::Measurement
	K
	i

	# Planet params
	μ
	α::Measurement
	gₚ
	Mₚ
	Rₚ
	ρₚ
	Tₚ
	H

	# Signal
	N_scales::Float64
	ΔD
end;

# ╔═╡ 855e7c4c-e0fe-11ea-1bbb-1b9db42a984d
md"### Structure to hold the inputs used to calculate each parameter"

# ╔═╡ 410f5804-e0ef-11ea-0576-e1692cd42b1b
@with_kw_noshow struct Derived_inputs @deftype NamedTuple
	# Reference name (e.g. Ciceri et al. 2015)
	name::String = "Custom"
	
	# Star Params
	inputs_ρₛ
	inputs_gₛ
	inputs_Mₛ
	inputs_Rₛ
	inputs_Tₛ
	inputs_Lₛ
	
	#Orbital params
	inputs_RₚRₛ
	inputs_P
	inputs_aRₛ
	inputs_a
	inputs_b
	inputs_K
	inputs_i
	
	# Planet params
	inputs_μ
	inputs_α
	inputs_gₚ
	inputs_Mₚ
	inputs_Rₚ
	inputs_ρₚ
	inputs_Tₚ
	inputs_H
end;

# ╔═╡ 3833772c-d63f-11ea-09b5-f36d68e512ea
begin
	results, results_inputs = [], []
	for st in studies
		# Calculate parameters
		params, params_inputs = calculate_params(st)
		
		# Store summary
		summary = Derived(; name=st.name, params...)
		push!(results, summary)
		
		# Store summary inputs
		summary_inputs = Derived_inputs(; name=st.name, params_inputs...)
		push!(results_inputs, summary_inputs)
	end
end;

# ╔═╡ 33fc58d0-dbd9-11ea-3c45-83f4b5a2a818
function display_summary(d::Derived, d_i::Derived_inputs)
	md"""
	###### **$(d.name):**
	**Star Params** \
	Rₛ $(keys(d_i.inputs_Rₛ)) = $(uconvert(u"Rsun", d.Rₛ)) \
	Mₛ $(keys(d_i.inputs_Mₛ)) = $(uconvert(u"Msun", d.Mₛ)) \
	Tₛ $(keys(d_i.inputs_Tₛ)) = $(uconvert(u"K", d.Tₛ)) \
	Lₛ $(keys(d_i.inputs_Lₛ)) = $(uconvert(u"Lsun", d.Lₛ)) \
	ρₛ $(keys(d_i.inputs_ρₛ)) = $(uconvert(u"g/cm^3", d.ρₛ)) \
	log gₛ (cm/s²) $(keys(d_i.inputs_gₛ)) = 
	$(log10(ustrip(uconvert(u"cm/s^2", d.gₛ))))
	
	**Orbital params** \
	K $(keys(d_i.inputs_K)) = $(uconvert(u"m/s", d.K)) \
	i $(keys(d_i.inputs_i)) = $(uconvert(u"°", d.i)) \
	RₚRₛ $(keys(d_i.inputs_RₚRₛ)) = $(uconvert(NoUnits, d.RₚRₛ)) \
	aRₛ $(keys(d_i.inputs_aRₛ)) = $(uconvert(NoUnits, d.aRₛ)) \
	P $(keys(d_i.inputs_P)) = $(uconvert(u"d", d.P)) \
	b $(keys(d_i.inputs_b)) = $(d.b)

	**Planet params** \
	μ $(keys(d_i.inputs_μ)) = $(uconvert(u"u", d.μ)) \
	α $(keys(d_i.inputs_α)) = $(uconvert(NoUnits, d.α)) \
	Rₚ $(keys(d_i.inputs_Rₚ)) = $(uconvert(u"Rjup", d.Rₚ)) \
	Mₚ $(keys(d_i.inputs_Mₚ)) = $(uconvert(u"Mjup", d.Mₚ)) \
	ρₚ $(keys(d_i.inputs_ρₚ)) = $(uconvert(u"g/cm^3", d.ρₚ)) \
	Tₚ $(keys(d_i.inputs_Tₚ)) = $(uconvert(u"K", d.Tₚ)) \
	gₚ $(keys(d_i.inputs_gₚ)) = $(uconvert(u"m/s^2", d.gₚ)) \
	H $(keys(d_i.inputs_H)) = $(uconvert(u"km", d.H))

	**Signal at $(d.N_scales) scale heights** \
	ΔD = $(d.N_scales * uconvert(NoUnits, d.ΔD) * 1e6) ppm
	"""
end;

# ╔═╡ 4bfaf322-dbd9-11ea-0449-87d9aa07311f
display_summary.(results, results_inputs)

# ╔═╡ 7db94ad6-dda1-11ea-2f33-1da144f1b7ad
md"Libraries for using things like physical constants and units."

# ╔═╡ Cell order:
# ╟─c9ac27ee-dac0-11ea-2a8c-2d144b034a82
# ╟─b2286b26-dac2-11ea-1ce0-c7da562aa641
# ╟─19b35ef4-dac3-11ea-2d25-97e5482ff6a0
# ╟─07db65d6-dd99-11ea-103b-33d0317af127
# ╠═17302d74-d63b-11ea-3de3-49f0df0554ca
# ╟─7cff6dfc-dd9f-11ea-1fdf-7b6aaa9435b4
# ╠═4bfaf322-dbd9-11ea-0449-87d9aa07311f
# ╟─38a61304-e0fe-11ea-14b2-17d9b9e13c7b
# ╠═3833772c-d63f-11ea-09b5-f36d68e512ea
# ╟─0b6821a4-dac3-11ea-27d7-911521f0d3c0
# ╟─f8281da6-dd9f-11ea-1b6c-d32702215397
# ╠═c01eb856-e0f9-11ea-01d5-07593189ce46
# ╟─49f75dea-dda0-11ea-1a85-bbdd4750b878
# ╠═3f79c516-da77-11ea-1f6b-d3e7191a95d8
# ╟─c5c5ea28-dd9e-11ea-1f89-5b1371831177
# ╠═33fc58d0-dbd9-11ea-3c45-83f4b5a2a818
# ╟─8e5811ae-dd9e-11ea-127e-b9812511492b
# ╠═db28dbd2-db12-11ea-28e4-2b6cf30bd102
# ╟─a8df7ad0-dd9e-11ea-2a6a-f16683371016
# ╠═bd752a9e-dd80-11ea-141c-779c5135d4d8
# ╟─855e7c4c-e0fe-11ea-1bbb-1b9db42a984d
# ╠═410f5804-e0ef-11ea-0576-e1692cd42b1b
# ╟─7db94ad6-dda1-11ea-2f33-1da144f1b7ad
# ╠═02bfa078-d62b-11ea-15df-d701431829b9
