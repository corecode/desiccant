library(deSolve)
library(reshape2)
library(ggplot2)

calcMs <- function(W, T, parms) {
    with (as.list(parms), {
        Psat <- ifelse(T > 100,
                       10**(8.14019-1810.94/(244.485+T)),
                       10**(8.07131-1730.63/(233.426+T)))/51.71493
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
        ## L: length (N * dz)
        ## Me_in: inlet air water mass fraction
        ## Te_in: inlet air temperature
        ## Hads(W): heat of adsorption
        ## RH(W): isotherm relative humidity

        N <- length(state)/2
        Wavg <- state[1:N]
        Ts <- state[N+1:N]

        #print(list(Wavg, Ts))

        A <- r**2 * pi                  # area
        p <- 2*r * pi                   # perimeter

        Ga <- rhoe * v
        mdotG <- A * Ga
        KGeff <- 14.53*Ga*Re**-0.41
        Ms <- calcMs(Wavg, Ts, parms)
        names(Ms) <- NULL
        cpl <- 1884                     # J/kg/K
        cb <- 4186 * Wavg + 921

        ce <- function(Me) cpl * Me + 1004 * (1-Me)
        hc <- function(Me) 0.683*Ga*Re**-0.51 * ce(Me)

        dzSolve <- function(t, state, parms) {
            Mez <- state[1]
            Tez <- state[2]
            Msz <- Ms[1+t*(N-1)/L]
            Tsz <- Ts[1+t*(N-1)/L]

            dMe_dz <- KGeff*p/mdotG * (Msz-Mez) * (1-Mez)
            dTe_dz <- -p/ce(Mez)/mdotG * (hc(Mez) + cpl*KGeff * (Msz-Mez)) * (Tez-Tsz)

            #print(c(as.integer(t*(N-1)/L), Mez, Tez, dMe_dz, dTe_dz))

            return(list(c(Me=dMe_dz, Te=dTe_dz)))
        }
        z_state <- ode(c(Me_in, Te_in), seq(0, L, length.out=N), dzSolve, parms=parms, method=c("rk4"))
        Me <- z_state[,2]
        Te <- z_state[,3]

        dWavg <- -KGeff*p/A/rhob * (Ms-Me)
        dTs <- p/A/rhob/cb * (hc(Me) * (Te-Ts) - Hads(Wavg)*KGeff * (Ms-Me))

        return(list(c(dWavg, dTs), c(Ms=Ms, Me=Me, Te=Te)))
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
    Me_in=0.008,
    Te_in=23.7,
    Wavg_0=0.00,
    T_0=23.7,
    Hads=function(W) ifelse(W<=0.15, -300*W+2095, 2050)*1000,
    RH=function(W) ifelse(W<=0.07, 1.235*W + 267.99*W**2 - 3170.7*W**3 + 10087.16*W**4, 0.3316 + 3.18*W)
)

params$M_0<- calcMs(params$Wavg_0, params$T_0, params)

solveDesiccant <- function(times, N, params) {
    state <- c(Wavg=rep(params$Wavg_0, N), Ts=rep(params$T_0, N))

    res <- ode.1D(state, times, desiccant, params, nspec=2, method=c("euler"))
    attr(res, "params") <- params
    res
}

plotDesiccant <- function(res) {
    data <- melt(as.data.frame(res), "time")
    data <- within(data, {
        variable <- as.character(variable)
        matches <- regmatches(variable, regexec("([^0-9]*)(.*)", variable))
        variable <- factor(sapply(matches, function(e) e[2]))
        z <- as.integer(sapply(matches, function(e) e[3]))
        rm(matches)
    })

    N <- attr(res, "dimens")
    ggplot(subset(data, z == N), aes(time, value)) + geom_line() + facet_wrap(~variable, scales="free")
}

adsorbedWater <- function(res) {
    params <- attr(res, "params")
    N <- attr(res, "dimens")
    boxWeight <- params$rhob*params$L*params$r**2*pi/N
    adsorbedMass <- res[,1+1:N]*boxWeight
    adsorbedMass
}
