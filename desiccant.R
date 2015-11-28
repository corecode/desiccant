calcMs <- function(W, T, parms) {
    with (as.list(parms), {
        Psat <- ifelse(T > 100,
                       10**(8.14019-1810.94/(244.485+T)),
                       10**(8.07131-1730.63/(233.426+T)))
        Ms <- 0.622*RH(W) * Psat / (P - 0.378*RH(W) * Psat)
        return(Ms)
    })
}

desiccant <- function(time, state, parms, N) {
    with (as.list(parms), {
        ## params:
        ## rhoe: density of (inlet) air
        ## rhob: density of bed
        ## v: velocity
        ## Re: Reynold's number
        ## r: bed radius
        ## P: inlet pressure
        ## N: number of Z steps
        ## L: length (N * dz)
        ## Me_in: inlet air water mass fraction
        ## Te_in: inlet air temperature
        ## Hads(W): heat of adsorption
        ## RH(W): isotherm relative humidity

        print(time)

        Wavg <- state[1:N]
        Ts <- state[N+1:N]
        Te <- state[2*N+1:N]
        Me <- state[3*N+1:N]

        dz <- L/N
        A <- r**2 * pi                  # area
        p <- 2*r * pi                   # perimeter

        Ga <- rhoe * v
        mdotG <- A * Ga
        KGeff <- 0.704*Ga*Re**-0.51
        Ms <- calcMs(Wavg, Ts, parms)
        cpl <- 1884                     # J/kg/K
        ce <- cpl * Me + 1004 * (1-Me)
        cb <- 4186 * Wavg + 921
        hc <- 0.683*Ga*Re**-0.51 * ce

        dWavg <- -KGeff*p/A/rhob * (Ms-Me)
        dTs <- p/A/rhob/cb * (hc * (Te-Ts) - Hads(Wavg)*KGeff * (Ms-Me))

        dMe_dz <- KGeff*p/mdotG * (Ms-Me) * (1-Me)
        dMe_diff <- diff(c(Me_in, Me))
        dMe <- dMe_dz*dz - dMe_diff
        dTe_dz <- -p/ce/mdotG * (hc+cpl*KGeff * (Ms-Me)) * (Te-Ts)
        dTe_diff <- diff(c(Te_in, Te))
        dTe <- dTe_dz*dz - dTe_diff

        return(list(c(Wavg=dWavg, Ts=dTs, Te=dTe, Me=dMe)))
    })
}


params <- list(
    rhoe=1.18,
    rhob=400.6,
    v=0.45,
    Re=109.47,
    r=0.0651,
    P=15,
    L=77.5e-3,
    Me_in=0.0097,
    Te_in=23.7,
    Hads=function(W) ifelse(W<=0.15, -300*W+2095, 2050),
    RH=function(W) ifelse(W<=0.07, 1.235*W + 267.99*W**2 - 3170.7*W**3 + 10087.16*W**4, 0.3316 + 3.18*W)
)

initial <- list(
    Wavg=0.0088,
    T=23.7
)

initial$M <- calcMs(initial$Wavg, initial$T, params)

N <- 100
state <- c(rep(initial$Wavg, N), rep(initial$T, N), rep(initial$T, N), rep(initial$M, N))
