"""Still pictures of a system in a matplotlib figure, for reports and books.

Each atom is a ball: a darker rim, the body and a highlight, drawn from the back
to the front so that nearer atoms cover farther ones, with far atoms faded towards
the background. The view, colours and sizes are those of the interactive player.
matplotlib is imported only when something is drawn.
"""
import numpy as np

# The player's colours: species in order, except water's oxygen (red) and hydrogen (white).
PALETTE = np.array([[41, 145, 184], [214, 122, 48], [64, 166, 102], [153, 102, 191],
                    [204, 77, 89], [128, 140, 153], [230, 232, 235]])/255
ELEMENTS = {'O': 4, 'H': 6}


def _colours(species):
    names = list(dict.fromkeys(species))
    kind = {name: ELEMENTS.get(name, names.index(name) % 6) for name in names}
    return PALETTE[[kind[name] for name in species]]


def _rotation(box, projection):
    """Rows: screen right, screen up, towards the viewer (the player's starting views)."""
    if projection == 'orthographic':
        # Face-on to the xz plane with z up, or with z across for a box long in z (a slab).
        if box[2] == max(box) and box[2] > 1.5*min(box):
            return np.array([[0., 0, 1], [1, 0, 0], [0, 1, 0]])
        return np.array([[1., 0, 0], [0, 0, 1], [0, -1, 0]])
    if max(box) > 1.5*min(box):
        # A long box, such as a slab: look across its long axis, which then runs across the view.
        toward = np.array([1, .35, 0]) if np.argmax(box) == 2 else np.array([0, .35, 1])
        toward /= np.linalg.norm(toward)
        right = np.array([toward[2], 0, -toward[0]])/np.hypot(toward[0], toward[2])
        return np.array([right, np.cross(toward, right), toward])
    cy, sy, cp, sp = np.cos(.35), np.sin(.35), np.cos(.3), np.sin(.3)
    return np.array([[cy, 0, -sy], [-sy*sp, cp, -cy*sp], [sy*cp, sp, cy*cp]])


def draw_atoms(positions, box, species, radii, *, bonds=(), walled=(False, False, False), framed=True,
               ax=None, projection='perspective', background='white'):
    """Draw atoms into ``ax`` (a new figure if None) and return the Axes.

    positions (N, D) and box (D,) are in the system's length unit; radii has one
    display radius per atom. ``walled`` marks the axes with walls; ``framed=False``
    leaves out the box (an open cluster).
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.collections import EllipseCollection, LineCollection
        from matplotlib.colors import to_rgb
    except ImportError as exc:
        raise ImportError('Drawing needs matplotlib; it comes with fys2160-simple-md[notebook].') from exc
    if projection not in ('perspective', 'orthographic'):
        raise ValueError('projection must be "perspective" or "orthographic".')
    positions = np.array(positions, dtype=float)
    box = np.asarray(box, dtype=float)
    N, D = positions.shape
    radii = np.broadcast_to(np.asarray(radii, dtype=float), (N,))
    periodic = ~np.asarray(walled[:D], dtype=bool)
    # Whole molecules: each bonded partner at its image nearest the first atom.
    for i, j in bonds:
        jump = np.where(periodic, box*np.rint((positions[j]-positions[i])/box), 0.) if framed else 0.
        positions[j] -= jump
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 5) if D == 2 else (5.5, 5))
    edges = []
    if D == 2:
        xy, r, shade = positions, radii, np.ones(N)
        if framed:
            corners = np.array([[0, 0], [box[0], 0], [box[0], box[1]], [0, box[1]]])
            for a in range(4):
                b = (a+1) % 4
                wall = walled[0] if a % 2 else walled[1]  # sides 1 and 3 are x = const
                ax.plot(*np.array([corners[a], corners[b]]).T, color='0.15' if wall else '0.6',
                        lw=2 if wall else .8, ls='-' if wall else (0, (4, 3)), zorder=1, solid_capstyle='projecting')
        order = np.arange(N)
    else:
        R = _rotation(box, projection)
        centre = box/2
        p = (positions-centre) @ R.T
        distance = 2.2*np.linalg.norm(box) if framed else 4*np.ptp(positions, axis=0).max()
        scale = distance/(distance-p[:, 2]) if projection == 'perspective' else np.ones(N)
        xy, r = p[:, :2]*scale[:, None], radii*scale
        order = np.argsort(p[:, 2], kind='stable')  # far first
        z = (p[:, 2]-p[:, 2].min())/max(np.ptp(p[:, 2]), 1e-12)
        shade = .55+.45*z  # far atoms fade towards the background
        if framed:
            corners = np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)], dtype=float)*box
            q = (corners-centre) @ R.T
            qs = distance/(distance-q[:, 2]) if projection == 'perspective' else np.ones(8)
            q2 = q[:, :2]*qs[:, None]
            back, front = [], []
            for a in range(8):
                for d in range(3):
                    b = a ^ (1 << (2-d))
                    if a < b:
                        (front if q[a, 2]+q[b, 2] > 1e-9 else back).append([q2[a], q2[b]])
            ax.add_collection(LineCollection(back, colors='0.65', linewidths=.8, zorder=1))
            ax.add_collection(LineCollection(front, colors='0.45', linewidths=.8, zorder=4))
            edges = q2
    if len(bonds):
        ax.add_collection(LineCollection([[xy[i], xy[j]] for i, j in bonds], colors='0.55',
                                         linewidths=1.2, zorder=2))
    base = _colours(species)
    fade = shade[:, None]
    behind = np.array(to_rgb(background))
    rim = fade*.55*base+(1-fade)*behind
    body = fade*.85*base+(1-fade)*behind
    light = fade*(.55*base+.45)+(1-fade)*behind
    o = order
    # Per atom a rim disc, the body and a highlight up and to the left, interleaved so
    # that a nearer atom covers all three layers of a farther one.
    offsets = np.stack([xy[o], xy[o], xy[o]+np.c_[-.3*r[o], .3*r[o]]], axis=1).reshape(-1, 2)
    sizes = np.stack([2*r[o], 1.72*r[o], .76*r[o]], axis=1).reshape(-1)
    colours = np.stack([rim[o], body[o], light[o]], axis=1).reshape(-1, 3)
    ax.add_collection(EllipseCollection(sizes, sizes, np.zeros(len(sizes)), units='xy', offsets=offsets,
                                        offset_transform=ax.transData, facecolors=colours, linewidths=0, zorder=3))
    pad = 1.2*r.max()
    lo, hi = xy.min(axis=0)-pad, xy.max(axis=0)+pad
    if len(edges):
        lo, hi = np.minimum(lo, edges.min(axis=0)), np.maximum(hi, edges.max(axis=0))
    elif D == 2 and framed:
        lo, hi = np.minimum(lo, 0), np.maximum(hi, box)
    ax.set_xlim(lo[0], hi[0])
    ax.set_ylim(lo[1], hi[1])
    ax.set_aspect('equal')
    if D == 2:
        ax.set_xlabel('x')
        ax.set_ylabel('y')
    else:
        ax.set_axis_off()
    return ax


def display_radii(species, model, radius=None):
    """One display radius per atom: ``radius`` (a number or one per species), else the player's sizes."""
    from .visualization import WATER_DISPLAY
    names = [str(name) for name in dict.fromkeys(species)]
    if radius is None:
        radius = {name: WATER_DISPLAY.get(name, 1.)/2 for name in names} if model == 'water-rigid' else .5
    if isinstance(radius, dict):
        missing = [name for name in names if name not in radius]
        if missing:
            raise ValueError(f'radius has no value for species {missing[0]!r}.')
        values = np.array([float(radius[name]) for name in species])
    else:
        values = np.full(len(species), float(radius))
    if not np.isfinite(values).all() or np.any(values <= 0):
        raise ValueError('radius must be positive and finite.')
    return values
