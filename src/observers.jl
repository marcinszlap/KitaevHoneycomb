# Observer API is in ITensorMPS (not ITensors) in your setup
mutable struct TruncErrObserver <: ITensorMPS.AbstractObserver
    last_truncerr::Float64
end
TruncErrObserver() = TruncErrObserver(NaN)

function ITensorMPS.measure!(obs::TruncErrObserver; truncerr=nothing, kwargs...)
    if truncerr !== nothing
        obs.last_truncerr = float(truncerr)
    end
    return nothing
end
