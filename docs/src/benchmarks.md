# Benchmarks

The large scale image-processing required for ADI algorithms can lead to concerns about runtime efficiency. To this end, ADI.jl (and the associated JuliaHCI packages) are developed with performance in mind. These packages do not aim to be as fast as possible; rather they focus on being as fast as *is convenient* (for the users and the devs).

### Meta Information

The benchmarks here can be found in the [`bench/`](https://github.com/JuliaHCI/ADI.jl/blob/master/bench/) folder organized into Julia files. The benchmarks utilize BenchmarkTools.jl, PyCall.jl with `virtualenv`, and CSV.jl for accuracy, reproducibility, and organization.

<details>
<summary>**System Information**</summary>
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
</details>

## ADI.jl vs. VIP-HCI

The [Vortex Imaging Pipeline](https://github.com/vortex-exoplanet/vip) (VIP) is the inspiration for ADI.jl. It is one of the major Python HCI packages and it offers many more features than ADI.jl. Some of the common uses for both packages include full-frame ADI processing, S/N maps, and contrast curves.

```@example bench
using CSV
using DataFrames
using StatsPlots
```

### ADI Reduction

These benchmarks show the duration to fully reduce ADI data for various algorithms. The data used are $\beta$ Pictoris and HR8799 from [HCIDatasets.jl](https://github.com/JuliaHCI/HCIDatasets.jl).

```@example bench
adi_data = CSV.File("assets/adi_benchmarks.csv") |> DataFrame |> sort!
```

```@example bench
groups = groupby(adi_data, :framework)
plot(
    @df groups[1] barplot(:alg, :time, group=:N),
    @df groups[2] barplot(:alg, :time, group=:N),
    legendtitle="Num. Pixels",
    size=(900, 500),
    dpi=75,
    ylabel="time (s)",
    title=["ADI.jl" "VIP"]
)
```


### Detection Maps

This benchmark measures the duration to produce a signal-to-noise ratio (S/N) map. Rather than test exact cubes, these benchmarks test randomly generated frames of various sizes. The FWHM is fixed at 5.

```@example bench
snrmap_data = CSV.File("assets/snrmap_benchmarks.csv") |> DataFrame |> sort!
```

```@example bench
@df snrmap_data scatter(
    :N,
    :time,
    group=:framework,
    msw=0,
    xlabel="number of pixels
    ylabel="time (s)"
)
```

### Contrast Curves

Finally, this benchmark measures the duration to generate a contrast curve for analyzing the algorithmic throughput of an ADI algorithm. For both benchmarks 3 azimuthal branches are used for throughput injections and a FWHM of 8. A Gaussian PSF function is used for the injections. The data used are $\beta$ Pictoris and HR8799 from [HCIDatasets.jl](https://github.com/JuliaHCI/HCIDatasets.jl).

```@example bench
contrast_data = CSV.File("assets/contrast_benchmarks.csv") |> DataFrame |> sort!
```

```@example bench
@df contrast_data scatter(
    :N,
    :time,
    group=:framework,
    msw=0,
    xlabel="number of pixels
    ylabel="time (s)"
)
```
