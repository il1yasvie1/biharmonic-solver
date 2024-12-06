import matplotlib.pyplot as plt


def plot_deformation(y, fname, equal_aspect=False):
    y1, y2, y3 = y.dat.data.T
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d', azim=-60, elev=30)
    trisurf = ax.plot_trisurf(
        y1, y2, y3,
        cmap='viridis',
        facecolors=plt.cm.viridis(y3),
        edgecolor='k',
        linewidth=0.2,
        antialiased=True,
        alpha=0.9,
        shade=True
    )
    ax.grid(False)
    ax.set_xlabel('y1')
    ax.set_ylabel('y2')
    ax.set_zlabel('y3')
    if equal_aspect:
        ax.set_aspect('equal')
    plt.savefig(fname, dpi=300, bbox_inches='tight')
