function deffTHz(cry)
    if cry == 0 # LN
        deff_ = 168e-12
    elseif cry == 2 # ZnTe
        deff_ = 0
    elseif cry == 3 # GaP
        deff_ = 0
    elseif cry == 4 # GaAs
        deff_ = 2 / sqrt(3) * 42.35e-12
    elseif cry == 7 # ZnSe
        deff_ = 0
    end
    return deff_
end

function neo(lambda, T, cry)
    if cry == 4 #GaAs Skauli et al. 2003 0.97-17 um
        l = lambda * 1e6
        a0 = 4.372514
        a = [5.466742 0.0242996 1.957522]
        b = [0.4431307 0.8746453 36.9166]

        n = real(sqrt.(Complex.(a0 .+ 1 .+ a[1] * l .^ 2 ./ (l .^ 2 .- b[1]^2) + a[2] * l .^ 2 ./ (l .^ 2 .- b[2]^2) + a[(3)] * l .^ 2 ./ (l .^ 2 .- b[(3)]^2))))
    end
    return n
end

function ngp(lambda, T, cry)
    lambda1 = lambda * 1e6
    if cry == 4
        @variables l

        a0 = 4.372514
        a = [5.466742 0.0242996 1.957522]
        b = [0.4431307 0.8746453 36.9166]

        n0 = real.(sqrt.(a0 + 1 + a[(1)] * l .^ 2 ./ (l .^ 2 - b[(1)]^2) + a[(2)] * l .^ 2 ./ (l .^ 2 - b[(2)]^2) + a[(3)] * l .^ 2 ./ (l .^ 2 - b[(3)]^2)))

        a = n0 - l * Symbolics.derivative(n0, l)
        #l = lambda1;
        ng = Symbolics.value(substitute(a, l => lambda1))

    end
    return ng
end

function nTHzo(omega, T, cry)
    if cry == 4

        nTHz = real.(sqrt.(er(omega, T, cry)))

    end
end

function diffegy_conv(z, A_kompozit, omega, T, k_omega, k_OMEGA, k_omegaSH, khi_eff, dnu, domega, k_omega0, omega0, gamma, cry)
    n2 = n2value(cry)
    deff_ = 1 * deff(cry)
    c = 3e8    #%m/s
    e0 = 8.854e-12
    abszorpcio = aTHzo(omega, T, cry)
    abszorpcio[abszorpcio.>1e5] .= 1e5
    ATHz = A_kompozit[:, 1]
    Aop = A_kompozit[:, 2]
    ASH = A_kompozit[:, 3]
    NN = length(Aop)
    At = ifft(Aop .* exp.(-0im * (k_omega) * z) * 2 * pi * dnu * length(omega))
    n2pm = fft(1im * e0 * omega0 * neo(2 * pi * 3e8 / omega0, T, cry) * n2 / 2 * abs.(At) .^ 2 .* At) / dnu / 2 / pi / length(omega) .* exp.(0im .* k_omega .* z)
    t1 = @spawn begin
        #      println(threadid())
        temp11 = conv(reverse(conj(Aop) .* exp.(1im .* k_omega .* z)), (Aop .* exp.(-1im * k_omega .* z)))
        temp11 = temp11[NN:end] .* exp.(1im .* k_OMEGA .* z) .* (-1 .* 1im .* khi_eff .* omega .^ 2 / 2 / c^2 ./ k_OMEGA) .* domega - 1 .* abszorpcio / 2 .* ATHz
        temp11[1] = 0
        temp1 = temp11
    end

    t2 = @spawn begin
        #       println(threadid())
        temp21 = conv(reverse(conj(ATHz) .* exp.(1im .* k_OMEGA .* z)), Aop .* exp.(-1im .* k_omega .* z))
        temp21 = temp21[NN:end] .* exp.(1im .* k_omega .* z)
        temp22 = conv(Aop .* exp.(-1im .* k_omega .* z), ATHz .* exp.(-1im .* k_OMEGA .* z))
        temp22 = temp22[1:NN] .* exp.(1im .* k_omega .* z)
        temp20 = -1 * n2pm - 1 * 1im * khi_eff .* omega .^ 2 / 2 / c^2 ./ k_omega .* (temp21 + temp22) .* domega
        temp20[1] = 0
        temp23 = conv(reverse(conj(Aop) .* exp.(1im .* k_omega .* z .* cos(gamma) .^ 2)), ASH .* exp.(-1im .* k_omegaSH .* z .* cos(gamma) .^ 2)) .* domega
        temp23 = -1 * cos(gamma) .* 1im .* deff_ .* omega .^ 2 / c^2 ./ k_omega .* temp23[NN:end] .* exp.(1im .* k_omega .* z .* cos(gamma) .^ 2)
        temp24 = temp20 + temp23
        temp24[1] = 0
        temp2 = temp24
    end
    t3 = @spawn begin
        #        println(threadid())
        temp31 = conv(Aop .* exp.(-1im .* k_omega .* z * cos(gamma)^2), Aop .* exp.(-1im .* k_omega .* z * cos(gamma)^2)) * domega
        temp31 = -1 * cos(gamma) * 1im * deff_ .* omega .^ 2 / 2 / c^2 ./ k_omega .* temp31[1:NN] .* exp.(1im .* k_omegaSH .* z .* cos(gamma) .^ 2)
        temp31[1] = 0
        temp3 = temp31
    end
    return cat(fetch(t1), fetch(t2), fetch(t3), dims=2)
end

function RK4_M(f, step, t0, y0, t_final)
    T = t0:step:t_final
    #println(T);
    Y = y0
    for ii in 2:length(T)
        k1 = f(T[ii-1], Y)
        k2 = f(T[ii-1] + step / 2, Y .+ k1 * step / 2)
        k3 = f(T[ii-1] + step / 2, Y .+ k2 * step / 2)
        k4 = f(T[ii-1] + step, Y .+ k3 * step)
        Y .+= 1 / 6 * (k1 .+ 2 * k2 .+ 2 * k3 .+ k4) * step
    end
    return T, Y
end

function n2value(cry)
    if cry == 4 # GaAs
        n2_ = 5.9e-18
    end
    return n2_
end

function deff(cry)
    if cry == 4 #% GaAs
        deff_ = 65.6e-12
    end
    return deff_
end

function aTHzo(omega, T, cry)
    if cry == 4
        alpha = -2 .* omega / 3e8 .* imag(sqrt.(er(omega, T, cry)))
    end
    return alpha
end
function er(omega, T, cry)
    nu = omega / 2 / pi / 3e8 * 0.01
    if cry == 4 #GaAs
        if T == 300 #ord
            e_inf = 11
            nu_L = 292.1
            nu_T = 268.7
            G = 2.4
            nu = omega / 2 / pi / 3e8 * 1e-2

            er_ = e_inf * (1 .+ (nu_L^2 .- nu_T^2) ./ (nu_T^2 .- nu .^ 2 .+ 1im * G * nu))
        end
    end
    return er_
end

function DataBaseWriter(FID, z, Aop, Iop, ATHz, ETHz, ASH, ISH)
    FID[string(Int(floor(z*1e6)))*"/Aop"] = collect(abs.(Aop))
    FID[string(Int(floor(z*1e6)))*"/Eop"] = collect(Iop)
    FID[string(Int(floor(z*1e6)))*"/ATHz"] = collect(abs.(ATHz))
    FID[string(Int(floor(z*1e6)))*"/ETHz"] = collect(ETHz)
    FID[string(Int(floor(z*1e6)))*"/ASH"] = collect(abs.(ASH))
    FID[string(Int(floor(z*1e6)))*"/ESH"] = collect(ISH)
end

function DataBaseEnder(FID, z, t, nu, effic, efficSH)
    FID["z"] = collect(transpose(z))
    FID["effic"] = collect(transpose(effic))
    FID["efficSH"] = collect(transpose(efficSH))
    FID["t"] = collect(transpose(t))
    FID["nu"] = collect(transpose(nu))
end

