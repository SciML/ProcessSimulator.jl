using ProcessSimulator
using Test

@testset "public API" begin
    public_names = (
        :Reaction,
        :MaterialSource,
        :Port,
        :MaterialStream,
        :SimpleAdiabaticCompressor,
        :SimpleIsobaricHeatExchanger,
        :CSTR,
    )

    @static if isdefined(Base, :ispublic)
        @test all(name -> Base.ispublic(ProcessSimulator, name), public_names)
    else
        @test all(name -> isdefined(ProcessSimulator, name), public_names)
    end
    @test all(name -> !isnothing(Base.Docs.doc(ProcessSimulator, name)), public_names)

    source = ProcessSimulator.MaterialSource("helium";
        Mw = 0.004,
        molar_density = (p, T, xᵢ; kwargs...) -> p / (8.314 * T),
        VT_enthalpy = (ϱ, T, xᵢ) -> 20.0 * T,
        VT_internal_energy = (ϱ, T, xᵢ) -> 12.0 * T,
        VT_entropy = (ϱ, T, xᵢ) -> log(T),
    )

    @test source.Mw == [0.004]
    @test nameof(ProcessSimulator.Port(source; name = :inlet)) == :inlet
    @test nameof(ProcessSimulator.MaterialStream(source; name = :stream)) == :stream
    @test nameof(ProcessSimulator.SimpleAdiabaticCompressor(source; name = :compressor)) ==
        :compressor
    @test nameof(ProcessSimulator.SimpleIsobaricHeatExchanger(source; name = :exchanger)) ==
        :exchanger
    @test nameof(ProcessSimulator.CSTR(source; name = :reactor)) == :reactor
end
