def callback(xs, ys, zs, np, ax, scatter):
    def x(i):
        scatter.set_offsets(np.c_(xs[i], ys[i], zs[i]))
        return None

    return x