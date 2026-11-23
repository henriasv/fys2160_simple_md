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


# Named views: (elevation, azimuth) in degrees, as in matplotlib's 3D plots, with z up.
VIEWS = {'angled': (20., -60.), 'front': (0., -90.), 'side': (0., 0.), 'top': (90., -90.)}


def _look(elevation, azimuth):
    """Rows: screen right, screen up, towards the viewer, for a camera at that
    elevation above the xy plane and azimuth from the x axis, with z up."""
    e, a = np.radians(elevation), np.radians(azimuth)
    toward = np.array([np.cos(e)*np.cos(a), np.cos(e)*np.sin(a), np.sin(e)])
    right = np.array([-np.sin(a), np.cos(a), 0.])
    up = np.array([-np.sin(e)*np.cos(a), -np.sin(e)*np.sin(a), np.cos(e)])
    return np.array([right, up, toward])


def _rotation(box, view, projection):
    """The view's rows: screen right, screen up, towards the viewer."""
    if view is None:
        if box[2] == max(box) and box[2] > 1.5*min(box):
            # A box long in z, such as a slab, lies across the picture as in the player:
            # z to the right, y up, seen from x (a little from above in perspective).
            if projection == 'orthographic':
                return np.array([[0., 0, 1], [0, 1, 0], [-1, 0, 0]])
            e = np.radians(20)
            toward = np.array([-np.cos(e), np.sin(e), 0.])
            return np.array([[0., 0, 1], np.cross(toward, [0., 0, 1]), toward])
        view = 'front' if projection == 'orthographic' else 'angled'
    if isinstance(view, str):
        if view not in VIEWS:
            raise ValueError(f'view must be one of {", ".join(VIEWS)}, or (elevation, azimuth) in degrees.')
        return _look(*VIEWS[view])
    try:
        elevation, azimuth = (float(angle) for angle in view)
    except (TypeError, ValueError):
        raise ValueError('view must be a name or (elevation, azimuth) in degrees.') from None
    return _look(elevation, azimuth)


def _arrows(ax, R, background):
    """Small x, y, z arrows near the lower left corner, foreshortened as the axes
    are, on a soft patch of the background colour. An axis along the line of sight
    is a circle: with a dot when it points towards the viewer, a cross away from it."""
    from matplotlib.patches import Circle
    from matplotlib.transforms import ScaledTranslation
    figure = ax.figure
    # Inches from the corner of the Axes, so the arrows keep their size in any figure.
    at = figure.dpi_scale_trans+ScaledTranslation(0, 0, ax.transAxes)
    centre, length = np.array([.42, .42]), .27
    ax.add_patch(Circle(centre, .38, transform=at, facecolor=background, edgecolor='none', alpha=.8,
                        zorder=5, clip_on=False))
    for k, name in enumerate('xyz'):
        v = R[:2, k]
        size = np.hypot(*v)
        if size < .2:
            ax.plot(*centre, marker='o', ms=7, mfc='none', mec='0.3', mew=.9, transform=at, clip_on=False, zorder=6)
            ax.plot(*centre, marker='.' if R[2, k] > 0 else 'x', ms=4 if R[2, k] > 0 else 3.5, color='0.3',
                    mew=.9, transform=at, clip_on=False, zorder=6)
            label = centre+[-.13, -.13]
        else:
            ax.annotate('', xy=centre+length*v, xycoords=at, xytext=centre, textcoords=at, zorder=6,
                        arrowprops=dict(arrowstyle='-|>', color='0.3', lw=.9, shrinkA=0, shrinkB=0, mutation_scale=7))
            label = centre+(length*size+.08)*v/size
        ax.text(*label, name, transform=at, ha='center', va='center', fontsize=8, color='0.25', zorder=6)


def draw_atoms(positions, box, species, radii, *, bonds=(), walled=(False, False, False), framed=True,
               ax=None, view=None, projection='perspective', arrows=True, background='white', tilt=None):
    """Draw atoms into ``ax`` (a new figure if None) and return the Axes.

    positions (N, D) and box (D,) are in the system's length unit; radii has one
    display radius per atom. ``walled`` marks the axes with walls; ``framed=False``
    leaves out the box (an open cluster). ``view`` and ``projection`` as in
    ``System.draw``; ``arrows`` adds x, y, z arrows to a 3D picture. ``tilt``
    (xy, xz, yz) leans the box edges as in LAMMPS.
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
        R = _rotation(box, view, projection)
        xy_, xz_, yz_ = (0., 0., 0.) if tilt is None else tilt
        cell = np.array([[box[0], 0, 0], [xy_, box[1], 0], [xz_, yz_, box[2]]])  # rows: the edges a, b, c
        centre = cell.sum(axis=0)/2
        p = (positions-centre) @ R.T
        distance = 2.2*np.linalg.norm(box) if framed else 4*np.ptp(positions, axis=0).max()
        scale = distance/(distance-p[:, 2]) if projection == 'perspective' else np.ones(N)
        xy, r = p[:, :2]*scale[:, None], radii*scale
        order = np.argsort(p[:, 2], kind='stable')  # far first
        z = (p[:, 2]-p[:, 2].min())/max(np.ptp(p[:, 2]), 1e-12)
        shade = .55+.45*z  # far atoms fade towards the background
        if framed:
            corners = np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)], dtype=float) @ cell
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
        if arrows:
            _arrows(ax, R, behind)
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
