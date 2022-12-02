## --- Import packages we'll need

    using StatGeochem, Statistics
    using Plots; gr(); default(fmt = :svg)

## --- Read data files

    kernel_sigma = 330 # Myr

    abbott = importdataset("data/Abbott.csv", '|', importas=:Dict)
    abbott["Age_sigma"] = kernel_sigma * ones(size(abbott["Age"]))

    herzberg = importdataset("data/Herzberg.csv",'|', importas=:Dict)
    herzberg["Age_sigma"] = kernel_sigma * ones(size(herzberg["Age"]))

    ganne = importdataset("data/Ganne_FeTi.csv",'|', importas=:Dict)
    ganne["Age_sigma"] = kernel_sigma * ones(size(ganne["Age"]))

    condieHM = importdataset("data/CondieHM.csv",'|', importas=:Dict)
    condieHM["Age_sigma"] = kernel_sigma * ones(size(condieHM["Age"]))

    condieDM = importdataset("data/CondieDM.csv",'|', importas=:Dict)
    condieDM["Age_sigma"] = kernel_sigma * ones(size(condieDM["Age"]))

    condieEM = importdataset("data/CondieEM.csv",'|', importas=:Dict)
    condieEM["Age_sigma"] = kernel_sigma * ones(size(condieEM["Age"]))

    ks2012 = importdataset("data/KellerSchoene.csv",'|', importas=:Dict)

    kimura = Dict()
    kimura["Age"]  = [68.051, 144.932, 345.154, 454.212, 497.127, 742.077, 802.843, 801.042, 1283.93, 1578.87, 1743.42, 1873.8, 2027.54, 2463.79, 2463.86, 2465.69, 2665.88, 2662.36, 2762.47, 2758.93, 2937.62, 3030.64, 3030.67, 3291.55, 3398.95, 3398.98, 3470.46]
    kimura["T"] = [1389.64, 1396.78, 1368.1, 1368.07, 1379.99, 1406.16, 1371.56, 1348.89, 1629.08, 1550.28, 1663.55, 1438.07, 1412.98, 1443.88, 1558.39, 1626.38, 1541.64, 1635.87, 1625.11, 1687.14, 1533.22, 1609.53, 1660.83, 1417.42, 1642.84, 1684.59, 1627.31]


## --- Plot Herzberg, Kimura and Abbott data and bootstrap-resampled averages

    h1 = plot( xlabel="Age (Ma)", ylabel="T (C)", framestyle=:box,
        xlims=(0,4000), ylims=(1200,1700))

    # Plot dots
    plot!(h1, ganne["Age"], ganne["T_p"], st=:scatter, msalpha=0, label="")
    plot!(h1, abbott["Age"], abbott["T_l"], # yerror=2*abbott["T_l_sigma"],
        st=:scatter, msalpha=0, label="")
    plot!(h1, herzberg["Age"], herzberg["T_p"], st=:scatter, color=:blue, msalpha=0, label="")
    plot!(h1, kimura["Age"], kimura["T"], st=:scatter, color=:darkblue, msalpha=0, label="")

    nresamplings = 10000
    nbins = 10
    tmin = 0
    tmax = 4000

    c = Dict(); m = Dict(); e = Dict();

    # Calculate binned means and uncertainties
    # (c = bincenters, m = mean, el = lower 95% CI, eu = upper 95% CI)
    k = invweight_age([ganne["Age"]; herzberg["Age"]])
    p = 1.0 ./ ((k .* median(5.0 ./ k)) .+ 1.0) # Keep rougly one-fith of the abbott in each resampling
    (c["primelts"], m["primelts"], e["primelts"]) = bin_bsr([ganne["Age"]; herzberg["Age"]],
        [ganne["T_p"]; herzberg["T_p"]], tmin, tmax, nbins, [ganne["Age_sigma"]; herzberg["Age_sigma"]], nresamplings, p)
    plot!(h1, c["primelts"], m["primelts"], yerror=2*e["primelts"], label="", color=:darkblue, mscolor=:auto)


    # Calculate binned means and uncertainties
    # (c = bincenters, m = mean, el = lower 95% CI, eu = upper 95% CI)
    k = invweight_age(abbott["Age"])
    p = 1.0 ./ ((k .* median(5.0 ./ k)) .+ 1.0) # Keep rougly one-fith of the abbott in each resampling
    (c["abbott"], m["abbott"], e["abbott"]) = bin_bsr(abbott["Age"],abbott["T_l"], tmin, tmax, nbins, abbott["Age_sigma"],nresamplings,p)
    plot!(h1, c["abbott"], m["abbott"], yerror=2*e["abbott"], label="", color=:darkred, mscolor=:auto)

    # Plot results
    display(h1)
    savefig(h1,"AbbottHerzbergKimuraT.pdf")


## --- Plot Condie data and bootstrap-resampled averages

    # Calculate binned means and uncertainties
    # (c = bincenters, m = mean, el = lower 95% CI, eu = upper 95% CI)
    nbins = 10

    h2 = plot( xlabel="Age (Ma)", ylabel="T (C)", framestyle=:box, xlims=(0,4000), ylims=(1300,1650))
    plot!(h2, condieEM["Age"], condieEM["T_g"], st=:scatter, color=:red, msalpha=0, label="", markersize=2)
    k = invweight_age(condieEM["Age"]); p = 1.0 ./ ((k .* median(5.0 ./ k)) .+ 1.0) # Keep rougly one-fith of the abbott in each resampling
    (c["condieEM"], m["condieEM"], e["condieEM"]) = bin_bsr(condieEM["Age"], condieEM["T_g"], tmin, tmax, nbins, condieEM["Age_sigma"], nresamplings, p)
    plot!(h2, c["condieEM"], m["condieEM"], yerror=2*e["condieEM"], label="", color=:red, mscolor=:auto)

    plot!(h2, condieDM["Age"], condieDM["T_g"], st=:scatter, color=:cyan, msalpha=0, label="", markersize=2)
    k = invweight_age(condieDM["Age"]); p = 1.0 ./ ((k .* median(5.0 ./ k)) .+ 1.0) # Keep rougly one-fith of the abbott in each resampling
    (c["condieDM"], m["condieDM"], e["condieDM"]) = bin_bsr(condieDM["Age"], condieDM["T_g"], tmin, tmax, nbins, condieDM["Age_sigma"], nresamplings, p)
    plot!(h2, c["condieDM"], m["condieDM"], yerror=2*e["condieDM"], label="", color=:cyan, mscolor=:auto)

    plot!(h2, condieHM["Age"], condieHM["T_g"], st=:scatter, color=:blue, msalpha=0, label="", markersize=2)
    k = invweight_age(condieHM["Age"]); p = 1.0 ./ ((k .* median(5.0 ./ k)) .+ 1.0) # Keep rougly one-fith of the abbott in each resampling
    (c["condieHM"], m["condieHM"], e["condieHM"]) = bin_bsr(condieHM["Age"], condieHM["T_g"], tmin, tmax, nbins, condieHM["Age_sigma"], nresamplings, p)
    plot!(h2, c["condieHM"], m["condieHM"], yerror=2*e["condieHM"], label="", color=:blue, mscolor=:auto)

    display(h2)
    savefig(h2,"CondieDMEMHM.pdf")

    k = invweight_age([condieHM["Age"]; condieDM["Age"]]); p = 1.0 ./ ((k .* median(5.0 ./ k)) .+ 1.0)
    (c["condie"], m["condie"], e["condie"]) = bin_bsr([condieHM["Age"]; condieDM["Age"]], [condieHM["T_g"]; condieDM["T_g"]], tmin, tmax, nbins, [condieHM["Age_sigma"]; condieDM["Age_sigma"]], nresamplings, p)


## --- Plot relative cooling curves

    h3 = plot(xlabel="Age (Ma)", ylabel="ΔT (C)", framestyle=:box, xlims=(0,4000), ylims=(-30, 230), legend=:topleft, fg_color_legend=:white)
    # plot!(h3, c["condieHM"], m["condieHM"] .- m["condieHM"][1], yerror=2*e["condieHM"], label="condieHM", color=:orange, mscolor=:auto)
    # plot!(h3, c["condieDM"], m["condieDM"] .- m["condieDM"][1], yerror=2*e["condieDM"], label="condieDM", color=:red, mscolor=:auto)
    plot!(h3, c["condie"], m["condie"] .- m["condie"][1], yerror=2*e["condie"], label="condie", color=:orangered, mscolor=:auto)
    plot!(h3, c["abbott"], m["abbott"] .- m["abbott"][1], yerror=2*e["abbott"], label="abbott", color=:purple, mscolor=:auto)
    plot!(h3, ks2012["Age"], ks2012["T_p"] .- ks2012["T_p"][1], yerror=2*ks2012["T_p_sigma"], label="ks2012", color=:darkred, mscolor=:auto)
    plot!(h3, c["primelts"], m["primelts"] .- m["primelts"][1], yerror=2*e["primelts"], label="primelts", color=:blue, mscolor=:auto)
    display(h3)
    savefig(h3,"Trel.pdf")


## --- Calculate boostrap-resampled averages in 100 Myr bins

    nbins = 39
    tmin = 0
    tmax = 3900

    k = invweight_age([ganne["Age"]; herzberg["Age"]]); p = 1.0 ./ ((k .* median(5.0 ./ k)) .+ 1.0)
    (c["primelts"], m["primelts"], e["primelts"]) = bin_bsr([ganne["Age"]; herzberg["Age"]], [ganne["T_p"]; herzberg["T_p"]], tmin, tmax, nbins, [ganne["Age_sigma"]; herzberg["Age_sigma"]], nresamplings, p)

    k = invweight_age(abbott["Age"]); p = 1.0 ./ ((k .* median(5.0 ./ k)) .+ 1.0)
    (c["abbott"], m["abbott"], e["abbott"]) = bin_bsr(abbott["Age"],abbott["T_l"], tmin, tmax, nbins, abbott["Age_sigma"],nresamplings,p)

    k = invweight_age(condieDM["Age"]); p = 1.0 ./ ((k .* median(5.0 ./ k)) .+ 1.0)
    (c["condieDM"], m["condieDM"], e["condieDM"]) = bin_bsr(condieDM["Age"], condieDM["T_g"], tmin, tmax, nbins, condieDM["Age_sigma"], nresamplings, p)

    k = invweight_age(condieHM["Age"]); p = 1.0 ./ ((k .* median(5.0 ./ k)) .+ 1.0)
    (c["condieHM"], m["condieHM"], e["condieHM"]) = bin_bsr(condieHM["Age"], condieHM["T_g"], tmin, tmax, nbins, condieHM["Age_sigma"], nresamplings, p)

    k = invweight_age(condieDM["Age"]); p = 1.0 ./ ((k .* median(5.0 ./ k)) .+ 1.0)
    (c["condieDM"], m["condieDM"], e["condieDM"]) = bin_bsr(condieDM["Age"], condieDM["T_g"], tmin, tmax, nbins, condieDM["Age_sigma"], nresamplings, p)

    k = invweight_age([condieHM["Age"]; condieDM["Age"]]); p = 1.0 ./ ((k .* median(5.0 ./ k)) .+ 1.0)
    (c["condie"], m["condie"], e["condie"]) = bin_bsr([condieHM["Age"]; condieDM["Age"]], [condieHM["T_g"]; condieDM["T_g"]], tmin, tmax, nbins, [condieHM["Age_sigma"]; condieDM["Age_sigma"]], nresamplings, p)


    Age = c["primelts"]
    T_ave = ones(nbins)
    T_ave_sigma = ones(nbins)
    for i=1:nbins
        T = [m["condieDM"][i], m["condieHM"][i], m["abbott"][i], m["primelts"][i]]
        T_sigma = [e["condieDM"][i], e["condieHM"][i], e["abbott"][i], e["primelts"][i]]
        (T_ave[i], T_ave_sigma[i], mswd) = awmean(T, T_sigma)
    end
    plot(Age,T_ave,yerror=2*T_ave_sigma,mscolor=:auto,label="",framestyle=:box)
    savefig("TrelNon-KS2012.pdf")


## --- Plot weighted-mean Trel

    Age = c["primelts"]
    T_ave = ones(nbins)
    T_ave_sigma = ones(nbins)
    for i=1:nbins
        T = [m["condieDM"][i], m["condieHM"][i], m["abbott"][i], m["primelts"][i], ks2012["T_p"][i]]
        T_sigma = [e["condieDM"][i], e["condieHM"][i], e["abbott"][i], e["primelts"][i], ks2012["T_p_sigma"][i]]
        (T_ave[i], T_ave_sigma[i], mswd) = awmean(T, T_sigma)
    end
    plot(Age,T_ave,yerror=2*T_ave_sigma,mscolor=:auto,label="",framestyle=:box)
    plot!(ks2012["Age"], ks2012["T_p"], yerror=2*ks2012["T_p_sigma"], label="", color=:darkred, mscolor=:auto)
    savefig("TrelCombined.pdf")


## --- Export weighted-mean Trel to csv

    CoolingAverage=Dict()
    CoolingAverage["Age"] = Age
    CoolingAverage["T"] = T_ave
    CoolingAverage["T_sigma"] = T_ave_sigma
    CoolingAverage["Trel"] = T_ave .- T_ave[1]
    CoolingAverage["Trel_sigma"] = T_ave_sigma
    CoolingAverage["elements"] = collect(keys(CoolingAverage))
    exportdataset(CoolingAverage, "data/CoolingAverage.csv", ',')

## --- End of File
