# This file contains interpolation functions that operate on the FD stencil.

# Set weighting coefficients
function set_weigths(
    wNW :: Matrix{Float64},
    wNE :: Matrix{Float64},
    wSW :: Matrix{Float64},
    wSE :: Matrix{Float64},
)
    # Get sizes
    nx, ny = size(wNW)

    # Set weights
    for idy in 1:ny
        for idx in 1:nx

            # North West
            (idx == 1 && idy == 1                             ) ? wNW[idx, idy] = 4.0 : nothing # Left corner
            (idx >= 2 && idx <=nx-1 && idy == 1               ) ? wNW[idx, idy] = 2.0 : nothing # Left edge
            (idx == 1               && idy >= 2 && idy <= ny-1) ? wNW[idx, idy] = 2.0 : nothing # Top edge
            (idx == nx              && idy >= 1               ) ? wNW[idx, idy] = 0.0 : nothing # Bottom edge
            (idx >= 1               && idy == ny              ) ? wNW[idx, idy] = 0.0 : nothing # Right edge

            # North East
            (idx == 1 && idy == ny                            ) ? wNE[idx, idy] = 4.0 : nothing # Right corner
            (idx >= 2 && idx <=nx-1 && idy == ny              ) ? wNE[idx, idy] = 2.0 : nothing # Right edge
            (idx == 1               && idy >= 2 && idy <= ny-1) ? wNE[idx, idy] = 2.0 : nothing # Top edge
            (idx == nx              && idy >= 1               ) ? wNE[idx, idy] = 0.0 : nothing # Bottom edge
            (idx >= 1               && idy == 1               ) ? wNE[idx, idy] = 0.0 : nothing # Left edge

            # South West
            (idx == nx && idy == 1                             ) ? wSW[idx, idy] = 4.0 : nothing # Left corner
            (idx >= 2  && idx <=nx-1 && idy == 1               ) ? wSW[idx, idy] = 2.0 : nothing # Left edge
            (idx == nx               && idy >= 2 && idy <= ny-1) ? wSW[idx, idy] = 2.0 : nothing # Bottom edge
            (idx == 1                && idy >= 1               ) ? wSW[idx, idy] = 0.0 : nothing # Top edge
            (idx >= 1                && idy == ny              ) ? wSW[idx, idy] = 0.0 : nothing # Right edge

            # South East
            (idx == nx && idy == ny                           ) ? wSE[idx, idy] = 4.0 : nothing # Right corner
            (idx >= 2 && idx <=nx-1 && idy == ny              ) ? wSE[idx, idy] = 2.0 : nothing # Right edge
            (idx == nx              && idy >= 2 && idy <= ny-1) ? wSE[idx, idy] = 2.0 : nothing # Bottom edge
            (idx == 1               && idy >= 1               ) ? wSE[idx, idy] = 0.0 : nothing # Top edge
            (idx >= 1               && idy == 1               ) ? wSE[idx, idy] = 0.0 : nothing # left edge
        end
    end

    # Return
    return nothing
end

# Interpolate from centers to vertices
function c2v(
    Ac   :: Matrix{Float64},
    Av   :: Matrix{Float64},
    AvNW :: Matrix{Float64},
    AvNE :: Matrix{Float64},
    AvSW :: Matrix{Float64},
    AvSE :: Matrix{Float64},
    wNW  :: Matrix{Float64},
    wNE  :: Matrix{Float64},
    wSW  :: Matrix{Float64},
    wSE  :: Matrix{Float64}
)
    # Get sizes
    nvx, nvy = size(Av)

    # Interpolate
    for idy in 1:nvy
        for idx in 1:nvx
            (idx >= 1 && idx <= nvx-1 && idy >= 1 && idy <= nvy-1) ? AvNW[idx, idy] = Ac[idx,   idy  ] : nothing # NW
            (idx >= 1 && idx <= nvx-1 && idy >= 2 && idy <= nvy  ) ? AvNE[idx, idy] = Ac[idx,   idy-1] : nothing # NE
            (idx >= 2 && idx <= nvx   && idy >= 1 && idy <= nvy-1) ? AvSW[idx, idy] = Ac[idx-1, idy  ] : nothing # SW
            (idx >= 2 && idx <= nvx   && idy >= 2 && idy <= nvy  ) ? AvSE[idx, idy] = Ac[idx-1, idy-1] : nothing # SE
            Av[idx, idy] = 0.25 * (wNW[idx, idy] * AvNW[idx, idy] + wNE[idx, idy] * AvNE[idx, idy] + wSW[idx, idy] * AvSW[idx, idy] + wSE[idx, idy] * AvSE[idx, idy])
        end
    end

    # Return
    return nothing
end

# Interpolate from vertices to centers
function v2c(
    Av :: Matrix{Float64},
    Ac :: Matrix{Float64}
)
    # Get sizes
    ncx, ncy = size(Ac)
    nvx, nvy = size(Av)
    
    # Interpolate
    for idy in 1:ncy
        for idx in 1:ncx
            Ac[idx, idy] = 0.25 * (Av[idx, idy] + Av[idx+1, idy] + Av[idx, idy+1] + Av[idx+1, idy+1])
        end
    end

    # Return
    return nothing
end

# Main function
function main()

    # Initialize sizes
    ncx, ncy = 6, 6
    nvx, nvy = ncx+1, ncy+1

    # Initialize Arrays
    Ac   =  ones(ncx, ncy)      # Cell center array
    Bc   = zeros(ncx, ncy)      # Cell center array
    Av   = zeros(nvx, nvy)      # Vertex array
    AvNW = zeros(nvx, nvy)      # Average tables for different locations
    AvNE = zeros(nvx, nvy)
    AvSW = zeros(nvx, nvy)
    AvSE = zeros(nvx, nvy)
    wNW  =  ones(nvx, nvy)      # Weights based on different locations
    wNE  =  ones(nvx, nvy)
    wSW  =  ones(nvx, nvy)
    wSE  =  ones(nvx, nvy)

    # Set weigthing coefficients
    set_weigths(wNW, wNE, wSW, wSE)

    # Fill test table in the cell center with some values
    Ac[[1, end], [1, end]] .= 3.0
    Ac[1       ,  2:end-1] .= 2.0
    Ac[end     ,  2:end-1] .= 2.0
    Ac[2:end-1 ,        1] .= 2.0
    Ac[2:end-1 ,      end] .= 2.0

    # Interpolate those values linearly to the vertices
    c2v(Ac, Av, AvNW, AvNE, AvSW, AvSE, wNW, wNE, wSW, wSE)
    v2c(Av, Bc)

    # Display results
    display(Ac)
    display(Av)
    display(Bc)

    # Return
    return nothing
end

# Run main function
main();
