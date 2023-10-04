# Solving 2D porous convection equation with central finite differences and fully implicit time integration
using Pkg, CairoMakie, Printf, FFMPEG
if isfile("Project.toml") && isfile("Manifest.toml")
    Pkg.activate(".")
end
fontsize_theme = Theme(fontsize = 25)
set_theme!(fontsize_theme)
# Define Function
@views function porousConvection_2D()
    # Physics
    lx, ly     = 40.0, 20.0            # Length in x and y
    k_ηf       = 1.0                   # Diffusion coefficient
    αρgx,αρgy  = 0.0,1.0
    αρg        = sqrt(αρgx^2+αρgy^2)
    ΔT         = 200.0
    ϕ          = 0.1
    Ra         = 1000
    λ_ρCp      = 1/Ra*(αρg*k_ηf*ΔT*ly/ϕ) # Ra = αρg*k_ηf*ΔT*ly/λ_ρCp/ϕ
    T_Bot      = ΔT/2.0                # BC at domain bottom
    T_Top      = -ΔT/2.0                #BC at domain top
    # Numerics
    nt         = 500 #ncx^2/5              # No. of time steps to compute
    re_D       = 4π
    CFL        = 1.0/sqrt(2.1)
    ncx, ncy   = 127,63               # Number of cells in x and y
    ϵtol       = 1e-8                 # Residual tolerance
    maxiter    = 20*max(ncx,ncy)      # Maximum no. iterations
    ncheck     = ceil(Int,0.25ncx)    # Convergence check frequency
    nvis       = 5                    # Visualisation frequency
    qstp       = 4                    # Quiver density: higher values = less dense
    aspRat     = lx/ly                # Model aspect ratio
    # Misc
    printFig   = true                # Printing switch
    writeMov   = true                # Write movie switch
    figName    = "./Solvers/PorousConvection/doc/png/PC_impl_BConFlux"    # Figure name
    movName    = "./Solvers/PorousConvection/doc/PC_impl_BConFlux_mov.mov"        # Movie name
    # Derived Numerics
    dx        = lx/ncx                # Spatial step size in x
    dy        = ly/ncy                # Spatial step size in y
    dt_diff   = min(dx,dy)^2/λ_ρCp/4.1 # Thermal diffusion time step
    θ_dτ_D    = max(lx,ly)/re_D/CFL/min(dx,dy)
    β_dτ_D    = (re_D*k_ηf)/(CFL*min(dx,dy)*max(lx,ly))
    xc        = LinRange(-lx/2.0+dx/2, lx/2.0-dx/2, ncx) # Coordinate array
    yc        = LinRange(-ly+dy/2,         0.0, ncy) # Coordinate array
    # Initialisation
    T         = @. ΔT*exp(-xc^2 - (yc'+ly/2.0)^2); T_ini = copy(T); T_old = copy(T)
    Tx_avg    = zeros(Float64, ncx-1, ncy  )
    Ty_avg    = zeros(Float64, ncx  , ncy-1)
    dTdt      = zeros(Float64, ncx  , ncy  )
    Pf        = zeros(Float64, ncx  , ncy  )
    qDx       = zeros(Float64, ncx+1, ncy  )
    qDy       = zeros(Float64, ncx  , ncy+1)
    qDx_avg   = zeros(Float64, ncx  , ncy  )
    qDy_avg   = zeros(Float64, ncx  , ncy  )
    qTx       = zeros(Float64, ncx+1, ncy  )
    qTy       = zeros(Float64, ncx  , ncy+1)
    r_Pf      = zeros(Float64, ncx  , ncy  )
    r_T       = zeros(Float64, ncx  , ncy  )
    # Plot initial configuration
    fig1      = Figure(figure_padding=30)                 # Plotting
    ax1       = Axis(fig1[1, 1], title="Porous convection 2D", aspect=aspRat, xlabel=L"\textit{xc}", ylabel=L"\textit{yc}", limits=(minimum(xc), maximum(xc), minimum(yc), maximum(yc)))
    heatmap!(ax1, xc, yc, T_ini', colormap=:heat)
    Colorbar(fig1[2,1], vertical=false, colormap=:heat, label=L"\textit{T}", limits=(-ΔT/2, ΔT/2))
    fig1
    # Time loop
    for it=1:nt
        T_old .= T
        # time stepping
        dt = if it == 1
            0.1*min(dx,dy)/(αρg*ΔT*k_ηf)
        else
            min(5.0*min(dx,dy)/(αρg*ΔT*k_ηf),ϕ*min(dx/maximum(abs.(qDx)), dy/maximum(abs.(qDy)))/2.1)
        end
        re_T    = π + sqrt(π^2 + ly^2/λ_ρCp/dt)
        θ_dτ_T  = max(lx,ly)/re_T/CFL/min(dx,dy)
        β_dτ_T  = (re_T*λ_ρCp)/(CFL*min(dx,dy)*max(lx,ly))
        # Calculate fluid pressure (implicit scheme)
        iter = 1; err_D = 2ϵtol; err_T = 2ϵtol; err_T_ini = 1.0; err_D_ini = 1.0
        while max(err_D, err_T) >= ϵtol && iter <= maxiter
            # ---------------------------------- #
            # Fluid pressure calculation         #
            # ---------------------------------- #
            Tx_avg .= (T[2:end,:] .+ T[1:end-1,:])./2.0             # Some averaging related to staggered grid
            Ty_avg .= (T[:,2:end] .+ T[:,1:end-1])./2.0
            #@printf("Tx_avg = %.2e; Ty_avg = %.2e\n",maximum(abs.(Tx_avg)),maximum(abs.(Ty_avg)))
            # Compute Darcy fluxes
            qDx[2:end-1,:].-=   1.0./(1.0 + θ_dτ_D).*(qDx[2:end-1,:] + k_ηf.*diff(Pf, dims=1)./dx .- αρgx.*Tx_avg) 
            qDy[:,2:end-1].-=   1.0./(1.0 + θ_dτ_D).*(qDy[:,2:end-1] + k_ηf.*diff(Pf, dims=2)./dy .- αρgy.*Ty_avg)
            # Fluid pressure update
            Pf            .-=   (diff(qDx, dims=1)./dx + diff(qDy, dims=2)./dy)./β_dτ_D
            # ---------------------------------- #
            # Temperature calculation            #
            # ---------------------------------- #
            # Boundary conditions
            qTx[[1,end],:]  .= 0.0;                                  # Left and right
            # qTy[:,1  ]      .= -λ_ρCp.*(2*(T_old[:,1  ].-T_Bot))./dy     # Bottom
            # qTy[:,end]      .= -λ_ρCp.*(2*(T_Top.-T_old[:,end]))./dy     # Top
            qTy[:,1]  .-= 1.0./(1.0 + θ_dτ_T).*(qTy[:,1] + λ_ρCp.*2*(T[:,1  ].-T_Bot)./dy) # Super important: the bounds must be update in the implicit fashion. An explicit update as in the lines above will not work!!!!!!!
            qTy[:,end].-= 1.0./(1.0 + θ_dτ_T).*(qTy[:,end] + λ_ρCp.*2*(T_Top.-T[:,end])./dy)
            # Compute fluxes
            qTx[2:end-1,:      ] .-= 1.0./(1.0 + θ_dτ_T).*(qTx[2:end-1,:      ] + λ_ρCp.*diff(T, dims=1)./dx)
            qTy[:      ,2:end-1] .-= 1.0./(1.0 + θ_dτ_T).*(qTy[:      ,2:end-1] + λ_ρCp.*diff(T, dims=2)./dy)
            # Temperature material time derivative
            dTdt                   .= (T-T_old)./dt 
            dTdt[2:end  ,:      ] .+= max.(0.0,qDx[2:end-1,:      ]).*diff(T,dims=1)./dx./ϕ
            dTdt[1:end-1,:      ] .+= min.(0.0,qDx[2:end-1,:      ]).*diff(T,dims=1)./dx./ϕ
            dTdt[:      ,2:end  ] .+= max.(0.0,qDy[:      ,2:end-1]).*diff(T,dims=2)./dy./ϕ
            dTdt[:      ,1:end-1] .+= min.(0.0,qDy[:      ,2:end-1]).*diff(T,dims=2)./dy./ϕ
            T                     .-= (dTdt .+ diff(qTx, dims=1)./dx .+ diff(qTy, dims=2)./dy)./(1.0/dt + β_dτ_T)
            # Convergence check
            if iter % ncheck == 0
                qDx_avg .= (qDx[2:end,:] .+ qDx[1:end-1,:]); qDy_avg .= (qDy[:,2:end] .+ qDy[:,1:end-1])
                r_Pf    .= diff(qDx, dims=1)./dx + diff(qDy, dims=2)./dy 
                r_T     .= dTdt .+ diff(qTx,dims=1)./dx .+ diff(qTy,dims=2)./dy
                err_D    = maximum(abs.(r_Pf))
                err_T    = maximum(abs.(r_T))
                #@printf("  iter/ncx=%.1f, err_D=%1.3e, err_T=%1.3e\n",iter/ncx,err_D,err_T)
            end
            # Store errors and update iteration count
            if iter==2 
                r_Pf     .= diff(qDx, dims=1)./dx + diff(qDy, dims=2)./dy 
                r_T      .= dTdt .+ diff(qTx,dims=1)./dx .+ diff(qTy,dims=2)./dy
                err_D_ini = maximum(abs.(r_Pf)); err_T_ini = maximum(abs.(r_T))
            end
            iter += 1
        end
        err_D_final = err_D; err_T_final = err_T
        # Print to stdout
        @printf("\nTime = %d\n\n", it)
        @printf("T : Initial residual = %1.3e; Final residual = %1.3e; iter/ncx = %.1f\n",err_T_ini,err_T_final,iter/ncx)
        @printf("Pf: Initial residual = %1.3e; Final residual = %1.3e; iter/ncx = %.1f\n",err_D_ini,err_D_final,iter/ncx)
        # Calculate temperature (explicit scheme)
        dt_adv = ϕ*min(dx/maximum(abs.(qDx)), dy/maximum(abs.(qDy)))/2.1
        dt     = min(dt_diff,dt_adv) # Time step calculation
        # Visualisation
        if it % nvis == 0
            sleep(1.0)
            empty!(ax1)
            heatmap!(ax1, xc, yc, T, colormap=:heat)
            arrows!(ax1, xc[1:qstp:end], yc[1:qstp:end], qDx_avg[1:qstp:end,1:qstp:end], qDy_avg[1:qstp:end,1:qstp:end], arrowsize=10, lengthscale=0.5, normalize=true)
            display(fig1)
            # Print figure
            if printFig
                save(figName*@sprintf("_%05d",it)*".png",fig1)
            end
        end
    end

    # Write movie
    if writeMov
        FFMPEG.ffmpeg_exe(`-framerate 12 -f image2 -pattern_type glob -i $(figName)_'*'.png -vf "scale=1920:960" -c:v libx264 -pix_fmt yuv420p -y "$(movName)"`)    
    end
end

# Call function
porousConvection_2D()