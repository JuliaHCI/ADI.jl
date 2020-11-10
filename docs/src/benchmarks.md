# Benchmarks

The large scale image-processing required for ADI algorithms can lead to concerns about runtime efficiency. To this end, *ADI.jl* (and the associated JuliaHCI packages) are developed with performance in mind. These packages do not aim to be as fast as possible; rather they focus on being as fast as *is convenient* (for the users and the devs).

The [Vortex Imaging Pipeline](https://github.com/vortex-exoplanet/vip) (VIP) is the inspiration for *ADI.jl*. It is one of the major Python HCI packages and it offers many more features than *ADI.jl*. Some of the common uses for both packages include full-frame ADI processing, S/N maps, and contrast curves.

### System/Setup Information

The benchmarks here can be found in the [`bench/`](https://github.com/JuliaHCI/ADI.jl/blob/master/bench/) folder organized into Julia files. The benchmarks utilize BenchmarkTools.jl, PyCall.jl with `virtualenv`, and CSV.jl for accuracy, reproducibility, and organization.

```
Julia Version 1.5.0
Commit 96786e22cc (2020-08-01 23:44 UTC)
Platform Info:
  OS: macOS (x86_64-apple-darwin18.7.0)
  CPU: Intel(R) Core(TM) i5-8259U CPU @ 2.30GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-9.0.1 (ORCJIT, skylake)
Environment:
  JULIA_NUM_THREADS = 4
```
For reproducibility, there is a `Manifest.toml` file in `bench/`. To reproduce this environment, first activate it, then instantiate it
```
$ julia --project=bench -e 'using Pkg; Pkg.instantiate()'
```
For the python code, there is a `requirements.txt` file in `bench/`. To reproduce this environment, (optionally) activate a virtual environment then install from the requirements file
```
(venv) $ pip install -r requirements.txt
```

```@setup bench
using CSV
using DataFrames
using StatsPlots
benchdir(args...) = joinpath("..", ".." ,"bench", args...);
```

## ADI Reduction

These benchmarks show the duration to fully reduce ADI data for various algorithms. The data used are $\beta$ Pictoris and HR8799 from [HCIDatasets.jl](https://github.com/JuliaHCI/HCIDatasets.jl).


```@example bench
adi_data = CSV.File(benchdir("adi_benchmarks.csv")) |> DataFrame |> sort!
cube_labels = @. ifelse(adi_data[:N] == 622261, "Beta Pictoris", "HR8799")
insertcols!(adi_data, 4, :cube => cube_labels)
adi_groups = groupby(adi_data, :framework)
```

```@example bench
cube_groups = groupby(adi_data, :cube)
plot(
    @df(cube_groups[1], groupedbar(:alg, :time, group=:framework, yscale=:log10)),
    @df(cube_groups[2], groupedbar(:alg, :time, group=:framework)),
    size=(700, 350),
    leg=:topleft,
    ylabel="time (s)",
    title=["Beta Pictoris" "HR8799"]
)
```


## Detection Maps

This benchmark measures the duration to produce a signal-to-noise ratio (S/N) map. Rather than test exact cubes, these benchmarks test randomly generated frames of various sizes. The FWHM is fixed at 5.

```@example bench
snrmap_data = CSV.File(benchdir("snrmap_benchmarks.csv")) |> DataFrame |> sort!
snrmap_groups = groupby(snrmap_data, :framework)
```

```@example bench
@df snrmap_data scatter(
    :N,
    :time,
    group=:framework,
    msw=0,
    ms=6,
    xlabel="number of pixels",
    ylabel="time (s)"
)
```

## Contrast Curves

Finally, this benchmark measures the duration to generate a contrast curve for analyzing the algorithmic throughput of an ADI algorithm. For both benchmarks 3 azimuthal branches are used for throughput injections and a FWHM of 8. A Gaussian PSF function is evaluated in a `(21, 21)` grid for the injections. The data used are $\beta$ Pictoris and HR8799 from [HCIDatasets.jl](https://github.com/JuliaHCI/HCIDatasets.jl).

```@example bench
contrast_data = CSV.File(benchdir("contrast_benchmarks.csv")) |> DataFrame |> sort!
cube_labels = @. ifelse(contrast_data[:N] == 622261, "Beta Pictoris", "HR8799")
insertcols!(contrast_data, 4, :cube => cube_labels)
contrast_groups = groupby(contrast_data, :framework)
```

```@example bench
@df contrast_data groupedbar(
    :cube,
    :time,
    group=:framework,
    leg=:topleft,
    ylabel="time (s)"
)
```
