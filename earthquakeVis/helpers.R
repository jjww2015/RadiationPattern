# calculate radiation patterns
bystep <- 6
rp <- function(u,v,sig,phi,delta,theta,eta) {
    # input 
    # u = displacement * area
    # v = volume explosion
    # sig = Poisson ratio
    # phi = strike in degree
    # delta = dip in degree
    # theta = tensile angle in degree
    # eta = rake angle in degree
    
    # convert degree to radian
    phi <- phi*pi/180.
    delta <- delta*pi/180.
    theta <- theta*pi/180.
    eta <- eta*pi/180.
    r <- 2*sig/(2*sig -1)
    
    # calculate radiation pattern
    azimuths <- seq(0,360,by = bystep)
    pols <- c(seq(0,180, by = bystep),180)
    
    x <- matrix(, nrow = length(azimuths), ncol = length(pols))
    y <- matrix(, nrow = length(azimuths), ncol = length(pols))
    z <- matrix(, nrow = length(azimuths), ncol = length(pols))
    rad <- matrix(, nrow = length(azimuths), ncol = length(pols))
    
    for(i in seq_along(azimuths)) {
        for(j in seq_along(pols)) {
            az <- azimuths[i] *pi/180.
            to <- pols[j] * pi/180.
            rad[i,j] <- -(2*cos(delta)*cos(theta)*sin(delta)*sin(eta) - 2*(((cos(az)*cos(delta)*cos(eta) - (2*sin(az)*sin(delta)^2 - sin(az))*sin(eta))*cos(phi) + (cos(delta)*cos(eta)*sin(az) + (2*cos(az)*sin(delta)^2 - cos(az))*sin(eta))*sin(phi))*r*cos(theta) + 2*(cos(delta)*cos(phi)*sin(az)*sin(delta) - cos(az)*cos(delta)*sin(delta)*sin(phi))*r*sin(theta))*cos(to)*sin(to) + 2*(((2*cos(az)^2 - 1)*cos(delta)*cos(phi)^2*sin(delta)*sin(eta) + (2*cos(az)^2 - 1)*cos(eta)*cos(phi)*sin(delta)*sin(phi) - (cos(az)^2 + 1)*cos(delta)*sin(delta)*sin(eta) - (2*cos(az)*cos(delta)*cos(phi)*sin(az)*sin(delta)*sin(eta)*sin(phi) - 2*cos(az)*cos(eta)*cos(phi)^2*sin(az)*sin(delta) + cos(az)*cos(eta)*sin(az)*sin(delta))*r)*cos(theta) - (2*r*cos(az)*cos(phi)*sin(az)*sin(delta)^2*sin(phi) - (2*cos(az)^2 - 1)*cos(phi)^2*sin(delta)^2 + (cos(az)^2 + 1)*sin(delta)^2 - 1)*sin(theta))*sin(to)^2 + (2*sin(delta)^2 + r - 2)*sin(theta))*u - (3*r - 2)*v
            #rsv <- (((cos(az)*cos(delta)*cos(eta) - (2*sin(az)*sin(delta)^2 - sin(az))*sin(eta))*cos(phi) + (cos(delta)*cos(eta)*sin(az) + (2*cos(az)*sin(delta)^2 - cos(az))*sin(eta))*sin(phi))*r*cos(theta) + 2*(cos(delta)*cos(phi)*sin(az)*sin(delta) - cos(az)*cos(delta)*sin(delta)*sin(phi))*r*sin(theta) - 2*(((2*sin(az)^2 - 1)*cos(delta)*sin(delta)*sin(eta)*sin(phi)^2 - (2*sin(az)^2 - 1)*cos(eta)*cos(phi)*sin(delta)*sin(phi) - (sin(az)^2 + 1)*cos(delta)*sin(delta)*sin(eta) - (2*cos(az)*cos(delta)*cos(phi)*sin(az)*sin(delta)*sin(eta)*sin(phi) + 2*cos(az)*cos(eta)*sin(az)*sin(delta)*sin(phi)^2 - cos(az)*cos(eta)*sin(az)*sin(delta))*r)*cos(theta) - (2*r*cos(az)*cos(phi)*sin(az)*sin(delta)^2*sin(phi) - (2*sin(az)^2 - 1)*sin(delta)^2*sin(phi)^2 + (sin(az)^2 + 1)*sin(delta)^2 - 1)*sin(theta))*cos(to)*sin(to) - 2*(((cos(az)*cos(delta)*cos(eta) - (2*sin(az)*sin(delta)^2 - sin(az))*sin(eta))*cos(phi) + (cos(delta)*cos(eta)*sin(az) + (2*cos(az)*sin(delta)^2 - cos(az))*sin(eta))*sin(phi))*r*cos(theta) + 2*(cos(delta)*cos(phi)*sin(az)*sin(delta) - cos(az)*cos(delta)*sin(delta)*sin(phi))*r*sin(theta))*sin(to)^2)*u
            #rsh <- -((((cos(delta)*cos(eta)*sin(az) + (2*cos(az)*sin(delta)^2 - cos(az))*sin(eta))*cos(phi) - (cos(az)*cos(delta)*cos(eta) - (2*sin(az)*sin(delta)^2 - sin(az))*sin(eta))*sin(phi))*r*cos(theta) - 2*(cos(az)*cos(delta)*cos(phi)*sin(delta) + cos(delta)*sin(az)*sin(delta)*sin(phi))*r*sin(theta))*cos(to) + ((4*cos(az)*cos(delta)*sin(az)*sin(delta)*sin(eta)*sin(phi)^2 - 4*cos(az)*cos(eta)*cos(phi)*sin(az)*sin(delta)*sin(phi) - 2*cos(az)*cos(delta)*sin(az)*sin(delta)*sin(eta) - (2*(2*cos(az)^2 - 1)*cos(delta)*cos(phi)*sin(delta)*sin(eta)*sin(phi) + 2*(2*cos(az)^2 - 1)*cos(eta)*sin(delta)*sin(phi)^2 - (2*cos(az)^2 - 1)*cos(eta)*sin(delta))*r)*cos(theta) - 2*((2*cos(az)^2 - 1)*r*cos(phi)*sin(delta)^2*sin(phi) - 2*cos(az)*sin(az)*sin(delta)^2*sin(phi)^2 + cos(az)*sin(az)*sin(delta)^2)*sin(theta))*sin(to))*u
            #rs <- sqrt(rsv*rsv + rsh*rsh)
            x[i,j] <- abs(rad[i,j])*sin(to)*cos(az)
            y[i,j] <- abs(rad[i,j])*sin(to)*sin(az)
            z[i,j] <- abs(rad[i,j])*cos(to)
        }
    }
    return(list(x=x,y=y,z=z,r=rad))
}

rsv <- function(u,v,sig,phi,delta,theta,eta) {
    # input 
    # u = displacement * area
    # v = volume explosion
    # sig = Poisson ratio
    # phi = strike in degree
    # delta = dip in degree
    # theta = tensile angle in degree
    # eta = rake angle in degree
    
    # convert degree to radian
    phi <- phi*pi/180.
    delta <- delta*pi/180.
    theta <- theta*pi/180.
    eta <- eta*pi/180.
    r <- 2*sig/(2*sig -1)
    
    # calculate radiation pattern
    azimuths <- seq(0,360,by =bystep)
    pols <- c(seq(0,180, by = bystep),180)
    
    x <- matrix(, nrow = length(azimuths), ncol = length(pols))
    y <- matrix(, nrow = length(azimuths), ncol = length(pols))
    z <- matrix(, nrow = length(azimuths), ncol = length(pols))
    rad <- matrix(, nrow = length(azimuths), ncol = length(pols))
    
    for(i in seq_along(azimuths)) {
        for(j in seq_along(pols)) {
            az <- azimuths[i] *pi/180.
            to <- pols[j] * pi/180.
            #rp <- -(2*cos(delta)*cos(theta)*sin(delta)*sin(eta) - 2*(((cos(az)*cos(delta)*cos(eta) - (2*sin(az)*sin(delta)^2 - sin(az))*sin(eta))*cos(phi) + (cos(delta)*cos(eta)*sin(az) + (2*cos(az)*sin(delta)^2 - cos(az))*sin(eta))*sin(phi))*r*cos(theta) + 2*(cos(delta)*cos(phi)*sin(az)*sin(delta) - cos(az)*cos(delta)*sin(delta)*sin(phi))*r*sin(theta))*cos(to)*sin(to) + 2*(((2*cos(az)^2 - 1)*cos(delta)*cos(phi)^2*sin(delta)*sin(eta) + (2*cos(az)^2 - 1)*cos(eta)*cos(phi)*sin(delta)*sin(phi) - (cos(az)^2 + 1)*cos(delta)*sin(delta)*sin(eta) - (2*cos(az)*cos(delta)*cos(phi)*sin(az)*sin(delta)*sin(eta)*sin(phi) - 2*cos(az)*cos(eta)*cos(phi)^2*sin(az)*sin(delta) + cos(az)*cos(eta)*sin(az)*sin(delta))*r)*cos(theta) - (2*r*cos(az)*cos(phi)*sin(az)*sin(delta)^2*sin(phi) - (2*cos(az)^2 - 1)*cos(phi)^2*sin(delta)^2 + (cos(az)^2 + 1)*sin(delta)^2 - 1)*sin(theta))*sin(to)^2 + (2*sin(delta)^2 + r - 2)*sin(theta))*u - (3*r - 2)*v
            rad[i,j] <- (((cos(az)*cos(delta)*cos(eta) - (2*sin(az)*sin(delta)^2 - sin(az))*sin(eta))*cos(phi) + (cos(delta)*cos(eta)*sin(az) + (2*cos(az)*sin(delta)^2 - cos(az))*sin(eta))*sin(phi))*r*cos(theta) + 2*(cos(delta)*cos(phi)*sin(az)*sin(delta) - cos(az)*cos(delta)*sin(delta)*sin(phi))*r*sin(theta) - 2*(((2*sin(az)^2 - 1)*cos(delta)*sin(delta)*sin(eta)*sin(phi)^2 - (2*sin(az)^2 - 1)*cos(eta)*cos(phi)*sin(delta)*sin(phi) - (sin(az)^2 + 1)*cos(delta)*sin(delta)*sin(eta) - (2*cos(az)*cos(delta)*cos(phi)*sin(az)*sin(delta)*sin(eta)*sin(phi) + 2*cos(az)*cos(eta)*sin(az)*sin(delta)*sin(phi)^2 - cos(az)*cos(eta)*sin(az)*sin(delta))*r)*cos(theta) - (2*r*cos(az)*cos(phi)*sin(az)*sin(delta)^2*sin(phi) - (2*sin(az)^2 - 1)*sin(delta)^2*sin(phi)^2 + (sin(az)^2 + 1)*sin(delta)^2 - 1)*sin(theta))*cos(to)*sin(to) - 2*(((cos(az)*cos(delta)*cos(eta) - (2*sin(az)*sin(delta)^2 - sin(az))*sin(eta))*cos(phi) + (cos(delta)*cos(eta)*sin(az) + (2*cos(az)*sin(delta)^2 - cos(az))*sin(eta))*sin(phi))*r*cos(theta) + 2*(cos(delta)*cos(phi)*sin(az)*sin(delta) - cos(az)*cos(delta)*sin(delta)*sin(phi))*r*sin(theta))*sin(to)^2)*u
            #rsh <- -((((cos(delta)*cos(eta)*sin(az) + (2*cos(az)*sin(delta)^2 - cos(az))*sin(eta))*cos(phi) - (cos(az)*cos(delta)*cos(eta) - (2*sin(az)*sin(delta)^2 - sin(az))*sin(eta))*sin(phi))*r*cos(theta) - 2*(cos(az)*cos(delta)*cos(phi)*sin(delta) + cos(delta)*sin(az)*sin(delta)*sin(phi))*r*sin(theta))*cos(to) + ((4*cos(az)*cos(delta)*sin(az)*sin(delta)*sin(eta)*sin(phi)^2 - 4*cos(az)*cos(eta)*cos(phi)*sin(az)*sin(delta)*sin(phi) - 2*cos(az)*cos(delta)*sin(az)*sin(delta)*sin(eta) - (2*(2*cos(az)^2 - 1)*cos(delta)*cos(phi)*sin(delta)*sin(eta)*sin(phi) + 2*(2*cos(az)^2 - 1)*cos(eta)*sin(delta)*sin(phi)^2 - (2*cos(az)^2 - 1)*cos(eta)*sin(delta))*r)*cos(theta) - 2*((2*cos(az)^2 - 1)*r*cos(phi)*sin(delta)^2*sin(phi) - 2*cos(az)*sin(az)*sin(delta)^2*sin(phi)^2 + cos(az)*sin(az)*sin(delta)^2)*sin(theta))*sin(to))*u
            #rs <- sqrt(rsv*rsv + rsh*rsh)
            x[i,j] <- abs(rad[i,j])*sin(to)*cos(az)
            y[i,j] <- abs(rad[i,j])*sin(to)*sin(az)
            z[i,j] <- abs(rad[i,j])*cos(to)
        }
    }
    return(list(x=x,y=y,z=z,r=rad))
}

rsh <- function(u,v,sig,phi,delta,theta,eta) {
    # input 
    # u = displacement * area
    # v = volume explosion
    # sig = Poisson ratio
    # phi = strike in degree
    # delta = dip in degree
    # theta = tensile angle in degree
    # eta = rake angle in degree
    
    # convert degree to radian
    phi <- phi*pi/180.
    delta <- delta*pi/180.
    theta <- theta*pi/180.
    eta <- eta*pi/180.
    r <- 2*sig/(2*sig -1)
    
    # calculate radiation pattern
    azimuths <- seq(0,360,by =bystep)
    pols <- c(seq(0,180, by = bystep),180)
    
    x <- matrix(, nrow = length(azimuths), ncol = length(pols))
    y <- matrix(, nrow = length(azimuths), ncol = length(pols))
    z <- matrix(, nrow = length(azimuths), ncol = length(pols))
    rad <- matrix(, nrow = length(azimuths), ncol = length(pols))
    
    for(i in seq_along(azimuths)) {
        for(j in seq_along(pols)) {
            az <- azimuths[i] *pi/180.
            to <- pols[j] * pi/180.
            #rp <- -(2*cos(delta)*cos(theta)*sin(delta)*sin(eta) - 2*(((cos(az)*cos(delta)*cos(eta) - (2*sin(az)*sin(delta)^2 - sin(az))*sin(eta))*cos(phi) + (cos(delta)*cos(eta)*sin(az) + (2*cos(az)*sin(delta)^2 - cos(az))*sin(eta))*sin(phi))*r*cos(theta) + 2*(cos(delta)*cos(phi)*sin(az)*sin(delta) - cos(az)*cos(delta)*sin(delta)*sin(phi))*r*sin(theta))*cos(to)*sin(to) + 2*(((2*cos(az)^2 - 1)*cos(delta)*cos(phi)^2*sin(delta)*sin(eta) + (2*cos(az)^2 - 1)*cos(eta)*cos(phi)*sin(delta)*sin(phi) - (cos(az)^2 + 1)*cos(delta)*sin(delta)*sin(eta) - (2*cos(az)*cos(delta)*cos(phi)*sin(az)*sin(delta)*sin(eta)*sin(phi) - 2*cos(az)*cos(eta)*cos(phi)^2*sin(az)*sin(delta) + cos(az)*cos(eta)*sin(az)*sin(delta))*r)*cos(theta) - (2*r*cos(az)*cos(phi)*sin(az)*sin(delta)^2*sin(phi) - (2*cos(az)^2 - 1)*cos(phi)^2*sin(delta)^2 + (cos(az)^2 + 1)*sin(delta)^2 - 1)*sin(theta))*sin(to)^2 + (2*sin(delta)^2 + r - 2)*sin(theta))*u - (3*r - 2)*v
            #rsv <- (((cos(az)*cos(delta)*cos(eta) - (2*sin(az)*sin(delta)^2 - sin(az))*sin(eta))*cos(phi) + (cos(delta)*cos(eta)*sin(az) + (2*cos(az)*sin(delta)^2 - cos(az))*sin(eta))*sin(phi))*r*cos(theta) + 2*(cos(delta)*cos(phi)*sin(az)*sin(delta) - cos(az)*cos(delta)*sin(delta)*sin(phi))*r*sin(theta) - 2*(((2*sin(az)^2 - 1)*cos(delta)*sin(delta)*sin(eta)*sin(phi)^2 - (2*sin(az)^2 - 1)*cos(eta)*cos(phi)*sin(delta)*sin(phi) - (sin(az)^2 + 1)*cos(delta)*sin(delta)*sin(eta) - (2*cos(az)*cos(delta)*cos(phi)*sin(az)*sin(delta)*sin(eta)*sin(phi) + 2*cos(az)*cos(eta)*sin(az)*sin(delta)*sin(phi)^2 - cos(az)*cos(eta)*sin(az)*sin(delta))*r)*cos(theta) - (2*r*cos(az)*cos(phi)*sin(az)*sin(delta)^2*sin(phi) - (2*sin(az)^2 - 1)*sin(delta)^2*sin(phi)^2 + (sin(az)^2 + 1)*sin(delta)^2 - 1)*sin(theta))*cos(to)*sin(to) - 2*(((cos(az)*cos(delta)*cos(eta) - (2*sin(az)*sin(delta)^2 - sin(az))*sin(eta))*cos(phi) + (cos(delta)*cos(eta)*sin(az) + (2*cos(az)*sin(delta)^2 - cos(az))*sin(eta))*sin(phi))*r*cos(theta) + 2*(cos(delta)*cos(phi)*sin(az)*sin(delta) - cos(az)*cos(delta)*sin(delta)*sin(phi))*r*sin(theta))*sin(to)^2)*u
            rad[i,j] <- -((((cos(delta)*cos(eta)*sin(az) + (2*cos(az)*sin(delta)^2 - cos(az))*sin(eta))*cos(phi) - (cos(az)*cos(delta)*cos(eta) - (2*sin(az)*sin(delta)^2 - sin(az))*sin(eta))*sin(phi))*r*cos(theta) - 2*(cos(az)*cos(delta)*cos(phi)*sin(delta) + cos(delta)*sin(az)*sin(delta)*sin(phi))*r*sin(theta))*cos(to) + ((4*cos(az)*cos(delta)*sin(az)*sin(delta)*sin(eta)*sin(phi)^2 - 4*cos(az)*cos(eta)*cos(phi)*sin(az)*sin(delta)*sin(phi) - 2*cos(az)*cos(delta)*sin(az)*sin(delta)*sin(eta) - (2*(2*cos(az)^2 - 1)*cos(delta)*cos(phi)*sin(delta)*sin(eta)*sin(phi) + 2*(2*cos(az)^2 - 1)*cos(eta)*sin(delta)*sin(phi)^2 - (2*cos(az)^2 - 1)*cos(eta)*sin(delta))*r)*cos(theta) - 2*((2*cos(az)^2 - 1)*r*cos(phi)*sin(delta)^2*sin(phi) - 2*cos(az)*sin(az)*sin(delta)^2*sin(phi)^2 + cos(az)*sin(az)*sin(delta)^2)*sin(theta))*sin(to))*u
            #rs <- sqrt(rsv*rsv + rsh*rsh)
            x[i,j] <- abs(rad[i,j])*sin(to)*cos(az)
            y[i,j] <- abs(rad[i,j])*sin(to)*sin(az)
            z[i,j] <- abs(rad[i,j])*cos(to)
        }
    }
    return(list(x=x,y=y,z=z,r=rad))
}

rs <- function(u,v,sig,phi,delta,theta,eta) {
    # input 
    # u = displacement * area
    # v = volume explosion
    # sig = Poisson ratio
    # phi = strike in degree
    # delta = dip in degree
    # theta = tensile angle in degree
    # eta = rake angle in degree
    
    # convert degree to radian
    phi <- phi*pi/180.
    delta <- delta*pi/180.
    theta <- theta*pi/180.
    eta <- eta*pi/180.
    r <- 2*sig/(2*sig -1)
    
    # calculate radiation pattern
    azimuths <- seq(0,360,by =bystep)
    pols <- c(seq(0,180, by = bystep),180)
    
    x <- matrix(, nrow = length(azimuths), ncol = length(pols))
    y <- matrix(, nrow = length(azimuths), ncol = length(pols))
    z <- matrix(, nrow = length(azimuths), ncol = length(pols))
    rad <- matrix(, nrow = length(azimuths), ncol = length(pols))
    
    for(i in seq_along(azimuths)) {
        for(j in seq_along(pols)) {
            az <- azimuths[i] *pi/180.
            to <- pols[j] * pi/180.
            #rp <- -(2*cos(delta)*cos(theta)*sin(delta)*sin(eta) - 2*(((cos(az)*cos(delta)*cos(eta) - (2*sin(az)*sin(delta)^2 - sin(az))*sin(eta))*cos(phi) + (cos(delta)*cos(eta)*sin(az) + (2*cos(az)*sin(delta)^2 - cos(az))*sin(eta))*sin(phi))*r*cos(theta) + 2*(cos(delta)*cos(phi)*sin(az)*sin(delta) - cos(az)*cos(delta)*sin(delta)*sin(phi))*r*sin(theta))*cos(to)*sin(to) + 2*(((2*cos(az)^2 - 1)*cos(delta)*cos(phi)^2*sin(delta)*sin(eta) + (2*cos(az)^2 - 1)*cos(eta)*cos(phi)*sin(delta)*sin(phi) - (cos(az)^2 + 1)*cos(delta)*sin(delta)*sin(eta) - (2*cos(az)*cos(delta)*cos(phi)*sin(az)*sin(delta)*sin(eta)*sin(phi) - 2*cos(az)*cos(eta)*cos(phi)^2*sin(az)*sin(delta) + cos(az)*cos(eta)*sin(az)*sin(delta))*r)*cos(theta) - (2*r*cos(az)*cos(phi)*sin(az)*sin(delta)^2*sin(phi) - (2*cos(az)^2 - 1)*cos(phi)^2*sin(delta)^2 + (cos(az)^2 + 1)*sin(delta)^2 - 1)*sin(theta))*sin(to)^2 + (2*sin(delta)^2 + r - 2)*sin(theta))*u - (3*r - 2)*v
            rsv <- (((cos(az)*cos(delta)*cos(eta) - (2*sin(az)*sin(delta)^2 - sin(az))*sin(eta))*cos(phi) + (cos(delta)*cos(eta)*sin(az) + (2*cos(az)*sin(delta)^2 - cos(az))*sin(eta))*sin(phi))*r*cos(theta) + 2*(cos(delta)*cos(phi)*sin(az)*sin(delta) - cos(az)*cos(delta)*sin(delta)*sin(phi))*r*sin(theta) - 2*(((2*sin(az)^2 - 1)*cos(delta)*sin(delta)*sin(eta)*sin(phi)^2 - (2*sin(az)^2 - 1)*cos(eta)*cos(phi)*sin(delta)*sin(phi) - (sin(az)^2 + 1)*cos(delta)*sin(delta)*sin(eta) - (2*cos(az)*cos(delta)*cos(phi)*sin(az)*sin(delta)*sin(eta)*sin(phi) + 2*cos(az)*cos(eta)*sin(az)*sin(delta)*sin(phi)^2 - cos(az)*cos(eta)*sin(az)*sin(delta))*r)*cos(theta) - (2*r*cos(az)*cos(phi)*sin(az)*sin(delta)^2*sin(phi) - (2*sin(az)^2 - 1)*sin(delta)^2*sin(phi)^2 + (sin(az)^2 + 1)*sin(delta)^2 - 1)*sin(theta))*cos(to)*sin(to) - 2*(((cos(az)*cos(delta)*cos(eta) - (2*sin(az)*sin(delta)^2 - sin(az))*sin(eta))*cos(phi) + (cos(delta)*cos(eta)*sin(az) + (2*cos(az)*sin(delta)^2 - cos(az))*sin(eta))*sin(phi))*r*cos(theta) + 2*(cos(delta)*cos(phi)*sin(az)*sin(delta) - cos(az)*cos(delta)*sin(delta)*sin(phi))*r*sin(theta))*sin(to)^2)*u
            rsh <- -((((cos(delta)*cos(eta)*sin(az) + (2*cos(az)*sin(delta)^2 - cos(az))*sin(eta))*cos(phi) - (cos(az)*cos(delta)*cos(eta) - (2*sin(az)*sin(delta)^2 - sin(az))*sin(eta))*sin(phi))*r*cos(theta) - 2*(cos(az)*cos(delta)*cos(phi)*sin(delta) + cos(delta)*sin(az)*sin(delta)*sin(phi))*r*sin(theta))*cos(to) + ((4*cos(az)*cos(delta)*sin(az)*sin(delta)*sin(eta)*sin(phi)^2 - 4*cos(az)*cos(eta)*cos(phi)*sin(az)*sin(delta)*sin(phi) - 2*cos(az)*cos(delta)*sin(az)*sin(delta)*sin(eta) - (2*(2*cos(az)^2 - 1)*cos(delta)*cos(phi)*sin(delta)*sin(eta)*sin(phi) + 2*(2*cos(az)^2 - 1)*cos(eta)*sin(delta)*sin(phi)^2 - (2*cos(az)^2 - 1)*cos(eta)*sin(delta))*r)*cos(theta) - 2*((2*cos(az)^2 - 1)*r*cos(phi)*sin(delta)^2*sin(phi) - 2*cos(az)*sin(az)*sin(delta)^2*sin(phi)^2 + cos(az)*sin(az)*sin(delta)^2)*sin(theta))*sin(to))*u
            rad[i,j] <- sqrt(rsv*rsv + rsh*rsh)
            x[i,j] <- abs(rad[i,j])*sin(to)*cos(az)
            y[i,j] <- abs(rad[i,j])*sin(to)*sin(az)
            z[i,j] <- abs(rad[i,j])*cos(to)
        }
    }
    return(list(x=x,y=y,z=z,r=rad))
    
}