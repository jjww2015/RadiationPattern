An interactive 3D earthquake energy radiation pattern calculator
===
author: Student
date: August 23 2016
transition: linear
font-import: http://fonts.googleapis.com/css?family=Slabo 27px
font-family: 'Slabo 27px'
width: 1920
height: 1080
<img src=pwave2.png height = 300>
</img>
<img src=swave2.png height = 300>
</img>

Background: Motivation
===
incremental: true

1. Seismic wave energy radiation is **import** to seismologists.
2. The computation of the radiation pattern is, however, **non-trivial**.
3. A small tool as **convenient** as a calculator would be greatly appreciated in the seismology community.
4. **Real-time** update of the pattern allows seismologist to interact and exam the spatial patten. 
5. As far as the author concerns, such a tool is yet **unavailable** in the public domain, thus the author decided to implement this tool for sharing.

Introduction: Earthquake Source
========================================================
incremental: true
left: 35%
<img src=Picture3.png height = 600>
</img>

***
The earquake source can be defined by the fault plane:

- The fault is defined by **strike** ($\phi$) and **dip** ($\delta$) in (a).
- The moving direction is defined by **tensile opening** 
angle ($\theta$) and **rake** angle ($\eta$), as shown in (b).
- The discontinuity across the fault plane is controled by **$u$**.
- The volumn explosion can be determined by **$v$** (not displayed).

The radiation pattern is calculated by the following function:


```r
rp <- function(u,v,sig,phi,delta,theta,eta) {
    # input 
    # u = displacement * area
    # v = volume explosion
    # sig = Poisson ratio
    # phi = strike in degree
    # delta = dip in degree
    # theta = tensile angle in degree
    # eta = rake angle in degree
    ... # detail implementation is omitted
    
return(list(x=x,y=y,z=z,r=rad))
}
```


The Shiny App
========================================================

The web tool can be access here: https://jjww2015.shinyapps.io/EarthquakeRadiationPatternVis/

<img src =shinyapp.png width = 1000>
</img>
* From top to bottom and left to right, the UI panel defines the following parameters
   1. discontinuity $u$,
   2. volume explosion $v$,
   3. poisson ratio $\sigma$,
   4. strike, dip, tensile, and rake angle in degree.

Earthquake Energy Radiation Pattern Visualization
========================================================
<font size = 6>
code snippet:

```r
# Only part of code is shown
input <- list(u = 1.0, v=0.0, sigma = 0.25, phi = 0.0, delta = 45.0, theta = 0.0, eta = 0.0)
surf3D(radP$x, radP$y, radP$z, colvar = radP$r, lighting = TRUE, shade = 1.0, 
       phi = 40, theta = 120, main = "P-wave", bty = "b2",xlab = "North", ylab = "East", zlab = "Down")
surf3D(radS$x, radS$y, radS$z, colvar = radS$r, lighting = TRUE, shade = 1.0, 
       phi = 40, theta = 120, main = "S-wave", bty = "b2",xlab = "North", ylab = "East", zlab = "Down")
```
</font>

<img src="radpatr-figure/unnamed-chunk-3-1.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" width="800px" style="display: block; margin: auto;" />

For interactive view please visit [shinyio](https://jjww2015.shinyapps.io/EarthquakeRadiationPatternVis/). Library RGL was used for the web tool. 

<h2 align = right>The End</h2>
