

include("run_pso_ffi.jl")

using Plots

DIMENSIONS=2


function my_f(x::Real, y::Real)
    return cos(x) * sin(y) + 1e-3 * (x - 1) * (x - 1) + 1e-3 * y * y
end

function my_f(x::Array{<:Real})
    return my_f(x[1], x[2])
end

function my_f(a_ptr::Ptr{Cdouble})::Cdouble
    x = unsafe_wrap(Array, a_ptr, DIMENSIONS; own = false)
    my_f(x)
end




T_MAX = 1000

params = RUN_PSO.PsoParams(
    @cfunction(my_f, Cdouble, (Ptr{Cdouble},)),
    1., 0.1, 0.2,
    1., 1.,
    DIMENSIONS, 5, T_MAX, 5,
    [-10., -10.], [10., 10.],
    [-10., -10.], [10., 10.],
    [-1., 8., 5., -3., 5., 6., 7., 3., -9., -2.]
)


x = LinRange(-10, 10, 40)
y = LinRange(-10, 10, 40)


pso = RUN_PSO.pso_alloc()

RUN_PSO.pso_init(params, pso)


RUN_PSO.pso_first_steps(pso)



surrogate = (x,y) -> RUN_PSO.surrogate_eval(pso, [x,y])

for i in 1:T_MAX

    RUN_PSO.pso_mainloop(pso)

    cur_pos = RUN_PSO.pso_current_positions(pso, params.pop_size, params.dimensions)
    cur_vel = RUN_PSO.pso_current_velocities(pso, params.pop_size, params.dimensions)

    p1 = plot()

    heatmap!(p1, x, y, my_f, alpha=0.6)
    scatter!(p1, cur_pos[1,:], cur_pos[2,:])

    p2 = plot()

    heatmap!(p2, x, y, surrogate, alpha=0.6)
    scatter!(p2, cur_pos[1,:], cur_pos[2,:])

    p = plot(p1, p2)

    savefig(p, "plots/t$(lpad(i, 6, '0')).png")

    for k in 1:params.pop_size
        xk = cur_pos[:, k]
        sur_x = RUN_PSO.surrogate_eval(pso, xk)
        fx = my_f(xk)
        println("x_$k($i): $xk")
        println("sur(x)= $sur_x")
        println("f(x)= $fx")
        println("v_$k($i): $(cur_vel[:, k])")
    end

end

