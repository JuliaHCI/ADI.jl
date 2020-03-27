using FITSIO
using PyPlot
using ADI

# plotly()

filepath = "/Users/miles/dev/research/WFIRST-CGI"

cube = read(FITS(joinpath(filepath, "cube_cent_crop_bpfix_skysub.fits"))[1])
angles = read(FITS(joinpath(filepath, "parallactic_angles.fits"))[1])

dc = DataCube(cube, angles)

@time Threads.@threads for ncomps in 1:size(cube, 3)
    reduced = reduce(PCA, dc, maxoutdim = ncomps)

    imshow(reduced, cmap = "inferno", title = "ncomps=$ncomps")
    savefig(p, joinpath(@__DIR__, "figures/HR8799_n$ncomps.png"))
end
