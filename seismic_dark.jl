using GLMakie, GeometryBasics
using LinearAlgebra, Colors, FileIO, Downloads
using CSV, DataFrames


#################################################################
########################## EARTHQUAKES ##########################
earthquakes = DataFrame(CSV.File("query.csv"));

## depth unit, km
function toCartesian(lon, lat; r = 1.02, cxyz = (0, 0, 0))
    x = cxyz[1] + (r + 1500_000) * cosd(lat) * cosd(lon)
    y = cxyz[2] + (r + 1500_000) * cosd(lat) * sind(lon)
    z = cxyz[3] + (r + 1500_000) * sind(lat)
    return (x, y, z) ./ 1500_000
end


# Earthquakes data for plot
lons, lats = earthquakes.longitude, earthquakes.latitude
depth = earthquakes.depth
mag = earthquakes.mag
toPoints3D = [Point3f([toCartesian(lons[i], lats[i];
    r = -depth[i] * 1000)...]) for i in 1:length(lons)]
ms = (exp.(mag) .- minimum(exp.(mag))) ./ maximum(exp.(mag) .- minimum(exp.(mag)));



#################################################################
############################# EARTH #############################

# Earth for Plot
n = 1024 ÷ 4 # 2048
θ = LinRange(0, pi, n)
φ = LinRange(-pi, pi, 2 * n)
xe = [cos(φ) * sin(θ) for θ in θ, φ in φ]
ye = [sin(φ) * sin(θ) for θ in θ, φ in φ]
ze = [cos(θ) for θ in θ, φ in φ];


# Earth image (make sure it is an RGBA image)
earth_img = load("2k_custom_earth_daymap_alpha.png")


#################################################################
############################# PLOT ##############################

# activate and theme
GLMakie.activate!()
set_theme!(backgroundcolor = :black)
fig=Figure(; resolution=(1000, 1000))

# scene
ax = LScene(fig[1, 1]; 
    show_axis=false,
    )

# title
Label(fig[1, 1, Top()], 
    "Seismic Events on Earth", 
    padding = (0, 0, 20, 0), 
    color= :white, 
    fontsize = 36)

# scatter earthquakes
sc = meshscatter!(ax, toPoints3D; 
    markersize = ms / 20 .+ 0.001, 
    color = mag,
    # colormap = :nuuk
    # colormap = :linear_bmy_10_95_c78_n256
    colormap = :plasma
    )

# earth globe
surface!(ax, xe, ye, ze; 
    color = earth_img*0.6,
    transparency=true
    )

# Colorbar
Colorbar(fig[1, 2], 
    sc, 
    label="Magnitudes", 
    labelpadding = 14,
    labelcolor= :white, 
    labelsize = 24,
    tickcolor = :white,
    ticklabelcolor =:white,
    height=Relative(0.5)
    )


save("seismic_earth.png", fig;)