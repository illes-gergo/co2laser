using Plots
using Symbolics
using FFTW
using DSP
using Base.Threads
using HDF5

include("fuggvenyek_full.jl")
@time begin
    cry = 4
    T = 300
    c = 3e8
    lambda0 = 10.6e-6
    N = 4e4
    tau = 2e-12
    I0 = 80e13
    khi_eff = 2 * deffTHz(cry)
    e0 = 8.854e-12
    nu0 = 1.5e12
    deltanu = nu0

    dz = 0.5e-5
    z_vegso = 8e-3
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
    FID = h5open("DB_full", "w")
    let A_loop = A_komp
        for ii in eachindex(z)[2:end]
            (z2, A_loop) = RK4_M(v6_fgv, dz, z[ii-1], A_loop, z[ii])
            local ATHz = A_loop[:, 1]
            local Aop = A_loop[:, 2]
            local ASH = A_loop[:, 3]
            global effic[ii] = sum(abs.(ATHz) .^ 2 .* FI) / pF
            global efficSH[ii] = sum(abs.(ASH) .^ 2 .* npSH) / pF

            local Iop = np0 * e0 * c / 2 * abs.((ifft(Aop .* exp.(-1im * (k_omega - k_omega0) * z[ii]))) * omegaMAX) .^ 2
            local ETHz = real.((ifft(FA .* ATHz .* exp.(-1im * (k_OMEGA - k_OMEGA0) * z[ii])))) * omegaMAX
            local ISH = abs.((ifft(ASH .* exp.(-1im * (k_omegaSH - k_omegaSH0) * z[ii]))) * omegaMAX) .^ 2
            DataBaseWriter(FID, z[ii], Aop, Iop, ATHz, ETHz, ASH, ISH)
            #println(ii)
        end
    end
end
DataBaseEnder(FID, z, t, nu, effic, efficSH);
close(FID);

#display(plot(z, effic))
println("Végeztem!")
