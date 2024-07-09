# Some structures and functions that help me nicely organize my figures

# ------------------------ #
#|      Cairo Makie       |#
# ------------------------ #

struct makie_plot_options
    fig_res     :: Float64
    line_width  :: Float64    
    font_size   :: Float64
    marker_size :: Float64
end

function makie_plot_options(;   fig_res     :: Float64 = 4.0,
                                line_width  :: Float64 = 2.5,
                                font_size   :: Float64 = 18.0,
                                marker_size :: Float64 = 12.0)
        # Return structure
        return makie_plot_options(fig_res, line_width, font_size, marker_size)
end