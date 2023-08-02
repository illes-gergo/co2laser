using Plots
using Symbolics
using FFTW
using DSP
using Base.Threads
using HDF5, DelimitedFiles

include("fuggvenyek_full.jl")
@time begin
    plotlyjs()
    cry = 4
    T = 300
    c = 3e8
    lambda0 = 1.4e-6
    N = 4e4
    tau = 1e-12
    I0 = 60e13
    khi_eff = 2 * deffTHz(cry)
    e0 = 8.854e-12
    nu0 = 1.5e12
    deltanu = nu0

    dz = 0.5e-5
    z_vegso = 60e-3
    z = 0:dz:z_vegso
    omega0 = 2 * pi * c / lambda0

    elochirp = 0 * 1 * z_vegso / 2

    omegaMAX = 5 * 2 * omega0

    dt = 2 * pi / omegaMAX
    t = (0:N-1) * dt
    t = t .- t[end] / 2

    domega = omegaMAX / N
    dnu = domega / 2 / pi

    omega = (0:N-1) * domega
    nu = omega / 2 / pi

    deltaOmega = 2 * sqrt(2 * log(2)) / tau

    lambda = 2 * pi * c ./ omega
    lambda[1] = lambda[2]
    ngp0 = ngp(lambda0, T, cry)
    np0 = neo(lambda0, T, cry)
    ngpSH = ngp(lambda0 / 2, T, cry)
    npSH = neo(lambda0 / 2, T, cry)

    nTHz = nTHzo(2 * pi * nu0, T, cry)
    vfTHz = c ./ nTHz

    gamma = acos.(ngp0 / nTHz)

    A0t = sqrt(2 * I0 / neo(lambda0, T, cry) / e0 / c)

    n_omega = neo(lambda, T, cry)
    k_OMEGA = real(omega .* nTHzo(omega, T, cry) / c)
    k_OMEGA0 = real(omega .* nTHzo(2 * pi * nu0, T, cry) / c)

    ddk_omega = -ngp0 .^ 2 / omega0 / c / np0 * tan(gamma)^2
    k_omega = real(1 / cos(gamma) .* (omega .* n_omega / c + 1 * (omega .- omega0) .^ 2 / 2 .* ddk_omega))
    # + 1e5;
    ddk_omegaSH = -ngpSH .^ 2 / omega0 / 2 / c / npSH * tan(gamma)^2
    k_omegaSH = real(1 / cos(gamma) .* (omega .* n_omega / c + 1 * (omega .- 2 * omega0) .^ 2 / 2 .* ddk_omegaSH))
    # + 1e5;
    k_omega0 = real.(1 ./ cos(gamma) .* (omega .* ngp0 / c))
    k_omegaSH0 = real.(1 ./ cos(gamma) .* (omega .* ngpSH / c))

    A0 = sqrt(2 * I0 / neo(lambda0, T, cry) / e0 / c) * tau / (2 * sqrt(2 * pi * log(2)))

    Aot = A0t * exp.(-2 * log(2) * t .^ 2 / tau^2) .* exp.(1im * omega0 * t)
    Aop = fft(Aot) .* exp.(1im * (k_omega - k_omega0) * elochirp) * dt / 2 / pi
    Aop0 = copy(Aop)

    ATHz = zeros(size(Aop))
    ASH = zeros(size(Aop))

    pF = sum(abs.(Aop) .^ 2) * np0

    FI = 4 * nTHz .^ 2 ./ (1 + nTHz) .^ 2
    FA = 2 * nTHz ./ (1 + nTHz)

    A_komp = cat(ATHz, Aop0, ASH, dims=2)
    A_komp0 = copy(A_komp)
    v6_fgv(z, A_kompozit) = diffegy_conv(z, A_kompozit, omega, T, k_omega, k_OMEGA, k_omegaSH, khi_eff, dnu, domega, k_omega0, omega0, gamma, cry)

    effic = zeros(size(z))
    efficSH = zeros(size(effic))

    vg = domega ./ diff(n_omega .* omega ./ c ./cos(gamma))
    dk2 = diff(1 ./vg)
    gvd = dk2./domega
    display(plot(lambda .* 1e6, [vg; 0] ./ c, xlims=(1, 2), #= ylims=(0.294, 0.305) =#))
    display(plot(lambda .* 1e6, [gvd;0;0] ./ c, xlims=(1, 2),  #= ylims=[-3,1].*1e-32 =#))

    vg2 = domega ./ diff(k_omega)
    dk22 = diff(1 ./vg2)
    gvd2 = dk22./domega
    display(plot(lambda .* 1e6, [vg2; 0] ./ c, xlims=(1, 2),#=  ylims=(0.2, 0.4) =#))
    display(plot(lambda .* 1e6, [gvd2;0;0] ./ c, xlims=(1, 2),  #= ylims=[-0.75,-0.38].*1e-31 =#))

    display(plot(lambda .* 1e6, [vg.-vg2; 0] ./ c, xlims=(1, 2), #= ylims=(0.294, 0.305) =#))
    writedlm("vg-no-tilt.txt",[lambda*1e6;;[vg;0]])
    writedlm("gvd-no-tilt.txt",[lambda*1e6;;[gvd;0;0]])

    writedlm("vg-with-tilt.txt",[lambda*1e6;;[vg2;0]])
    writedlm("gvd-with-tilt.txt",[lambda*1e6;;[gvd2;0;0]])
end
#display(plot(z, effic))
println("Végeztem!")
