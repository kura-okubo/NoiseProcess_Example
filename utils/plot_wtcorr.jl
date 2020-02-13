using SeisIO, SeisNoise, JLD2, Plots, Printf

export sort_pairs, get_corrtype

function plot_wtcorr(C::CorrData; figdir::String="./fig_wtcorr")

	Plots.pyplot()
	
	if !ispath(figdir); mkdir(figdir); end

	freqband = C.misc["freqband"]
	Nfreqband = length(freqband) - 1 

	fontsize = 12

	T, N = size(C.corr)
	tvec = -C.maxlag:1/C.fs:C.maxlag
	Plots.plot(bg = :white)


	for i=1:Nfreqband
		tr = C.corr[:, i]
		normamp = maximum(abs.(tr))
		tr ./= normamp

		Plots.plot!(tr .+ 2*(float(i)-1), tvec,
				linewidth = 1,
				linecolor = :black,
				legend=false,
				xlabel = "short-time window cc functions", ylabel = "Time lag",
				xticks = false,
				xtickfontsize=fontsize,
				ytickfontsize=fontsize)
	end
	Plots.plot!(size=(800, 600), show=false)
	Plots.savefig(figdir*"/raw_cc.png")

	# plot wavelet transform
	Plots.plot3d(bg = :white)
	for ifreq = Nfreqband:-1:1
		for j = 1:12
			tr = C.misc["wtcorr"][:, j, ifreq]
			normamp = maximum(abs.(tr))
			tr ./= normamp

			yp = ones(length(tvec), 1)
			Plots.plot3d!(tr .+ 2*(float(j)-1), yp .* (float(ifreq) - 1),  tvec,
					linewidth = 1,
					linecolor = :auto,
					legend=false,
					xlabel = "short-time window cc functions",
					ylabel = "Frequency bands",
					yticks = false,
					zlabel = "Time lag",
					xtickfontsize=fontsize,
					ytickfontsize=fontsize)
		end
	end
	Plots.plot3d!(camera=(15, 40), size=(800, 600), show=false)
	Plots.savefig(figdir*"/wt_cc.png")
end
