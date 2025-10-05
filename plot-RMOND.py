from getdist import plots

g=plots.get_subplot_plotter(chain_dir=r'/mnt/f5e499d6-0bd8-4c53-b406-34024468918c/kazemrezazadeh/Kazem Rezazadeh/0001 Study and Research/0002 Research/0001 Theoretical physics/0025 Relativistic MOND Theory from Modified Entropic Gravity/0002 Computational file/0005 MCMC-RMOND/MCMC-RMOND')
g.settings.fontfamily = "Times"
g.settings.axes_fontsize = 16
g.settings.axes_labelsize = 18
g.settings.colorbar_axes_fontsize = 18
g.settings.fontsize = 18
g.settings.legend_fontsize = 18
roots = ['chains']
params = ['x1', 'x2', 'x5', 'x6', 'x7', 'x9']
g.triangle_plot(roots, params, filled=True)
g.export()
