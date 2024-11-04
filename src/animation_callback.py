def callback(xs, ys, zs, np, ax):
    def x(i):
        return [ax.scatter(xs[i], ys[i], zs[i], c="red")]

    return x