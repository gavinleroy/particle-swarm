module RUN_PSO



#####
# FFI helpers
#####


# Changes the major order of the array memory layout.
# Row -> Column and vice-versa
# This functions works on *any* sized 2-dimensional array
# and not just square matrices.
change_major_order(X::AbstractArray, sizes...=size(X)...) =
    permutedims(reshape(X, sizes...), length(sizes):-1:1)

column_major_to_row(M::AbstractArray)::AbstractVector =
    collect(reshape(change_major_order(M), length(M)))

array_to_column_major(A::AbstractArray, sizes...)::AbstractArray =
    collect(change_major_order(reshape(A, sizes)))


#####
# FFI Functions:
#####

struct PsoParams 
    f::Ptr{Cvoid}
    intertia::Real
    social::Real
    cognition::Real
    local_refinement_box_size::Real
    min_minimizer_distance::Real
    dimensions::Integer
    pop_size::Integer
    time_max::Integer
    n_trials::Integer
    bounds_low::Vector{<:Real}
    bounds_high::Vector{<:Real}
    vmin::Vector{<:Real}
    vmax::Vector{<:Real}
    initial_positions::Vector{<:Real}
end


# Abstract pso object
struct Pso
end


function pso_alloc()
    ccall(
        (:alloc_pso_data_constant_inertia, :libpso),
        Ptr{Pso},
        ()
    )
end



function pso_init(params::PsoParams, pso::Ptr{Pso})
    ccall(
        (:pso_constant_inertia_init, :libpso),
        Cvoid,
        (
            Ptr{Pso}, Ptr{Cvoid},
            Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
            Cint, Cint, Cint, Cint,
            Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}
        ),
        pso, params.f,
        params.intertia, params.social, params.cognition,
        params.local_refinement_box_size, params.min_minimizer_distance,
        params.dimensions, params.pop_size, params.time_max, params.n_trials,
        params.bounds_low, params.bounds_high, params.vmin, params.vmax,
        params.initial_positions
    )
end

function pso_first_steps(pso::Ptr{Pso})
    ccall(
        (:pso_constant_inertia_first_steps, :libpso),
        Cvoid,
        (Ptr{Pso},),
        pso
    )
end

function pso_mainloop(pso::Ptr{Pso})
    ccall(
        (:pso_constant_inertia_loop, :libpso),
        Cvoid,
        (Ptr{Pso},),
        pso
    )
end


function pso_current_positions(pso::Ptr{Pso}, popsize, dimensions)
    pos = ccall(
        (:get_current_particle_positions, :libpso),
        Ptr{Cdouble},
        (Ptr{Pso},),
        pso
    )
    # columns are vectors
    pos_array = unsafe_wrap(Array, pos, (dimensions, popsize); own = false)
    pos_array
end

function pso_current_velocities(pso::Ptr{Pso}, popsize, dimensions)
    vel = ccall(
        (:get_current_particle_velocities, :libpso),
        Ptr{Cdouble},
        (Ptr{Pso},),
        pso
    )
    vel_array = unsafe_wrap(Array, vel, (dimensions, popsize); own = false)
    copy(vel_array)
end



function surrogate_eval(pso::Ptr{Pso}, x::Vector{<:Real})
    ccall(
        (:surrogate_eval, :libpso),
        Cdouble,
        (Ptr{Pso},Ptr{Cdouble}),
        pso, x
    )
end



end #module