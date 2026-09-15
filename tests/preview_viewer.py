"""Write a standalone crystal/depth-check page: python preview_viewer.py /tmp/depth.html."""
import json
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
cells, rho = 5, 1.1
length = (4*cells**3/rho)**(1/3)
basis = [(0,0,0), (0,.5,.5), (.5,0,.5), (.5,.5,0)]
x = [[(v+.25+b[d])*length/cells for d,v in enumerate((i,j,k))]
     for i in range(cells) for j in range(cells) for k in range(cells) for b in basis]
data = dict(frames=[dict(x=x,box=[length]*3,step=0,time=0)],zoom=1.4,
            projection='orthographic',note='Dense FCC crystal: drag or use arrow keys to rotate.')
page = ('<!doctype html><meta charset="utf-8"><title>MD sphere-depth check</title>'
        '<body style="font:14px system-ui"><div id="view"></div>'
        '<button id="check">Check sphere depth</button><pre id="results"></pre><script>'
        + (root/'src/fys2160_md/viewer.js').read_text()
        + '\ncreateMDView(document.getElementById("view"),'+json.dumps(data)+');\n'
        + Path(__file__).with_name('viewer_depth_browser.js').read_text()+'</script>')
Path(sys.argv[1]).write_text(page)
