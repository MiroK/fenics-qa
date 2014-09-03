def compute_norm_at_vertices(w, target):
    wx, wy  = w.split(deepcopy=True)            # copy? Can we avoid it?
    wnorm2 = wx.vector() * wx.vector()  \
             + wy.vector() * wy.vector()        # doesn't work in parallel
    wnorm = np.sqrt(wnorm2.get_local())         # sqrt from numpy - avoid?
    target.vector().set_local(wnorm)
