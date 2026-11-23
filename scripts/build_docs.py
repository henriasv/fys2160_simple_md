"""Execute the short showcase and export static HTML with working saved players."""
import argparse
import html
from pathlib import Path
import shutil
import tempfile
import nbformat
from nbclient import NotebookClient
import markdown
from pygments import highlight
from pygments.lexers import PythonLexer
from pygments.formatters import HtmlFormatter

ROOT=Path(__file__).resolve().parents[1]

def export(nb):
    formatter=HtmlFormatter(cssclass="code")
    cells=[]
    for cell in nb.cells:
        if cell.cell_type=="markdown":
            cells.append('<section class="markdown">'+markdown.markdown(cell.source,extensions=["fenced_code","tables"])+"</section>")
        elif cell.cell_type=="code":
            cells.append('<section class="cell">'+highlight(cell.source,PythonLexer(),formatter))
            for output in cell.get("outputs",[]):
                data=output.get("data",{})
                if "text/html" in data:
                    text=data["text/html"]
                    cells.append('<div class="output">'+("".join(text) if isinstance(text,list) else text)+"</div>")
                elif "image/png" in data:
                    cells.append('<img class="plot" alt="Recorded simulation plot" src="data:image/png;base64,'+data["image/png"]+'">')
                elif output.output_type=="stream":
                    cells.append('<pre class="text-output">'+html.escape(output.text)+"</pre>")
                elif "text/plain" in data:
                    cells.append('<pre class="text-output">'+html.escape(data["text/plain"])+"</pre>")
                elif output.output_type=="error":
                    raise RuntimeError("Notebook contains an error output")
            cells.append("</section>")
    css="""
    *{box-sizing:border-box}body{margin:0;font:16px/1.65 system-ui,sans-serif;color:#203c43;background:#fff}
    header{padding:18px 5%;border-bottom:1px solid #dce8e4;color:#18756a;font-weight:650;font-size:13px}
    main{max-width:1060px;margin:auto;padding:22px 5% 70px}h1{font-size:32px;line-height:1.2;letter-spacing:-.03em;color:#123e40}h2{margin-top:42px;font-size:23px}
    .markdown{max-width:820px}.cell{margin:22px 0}.code{padding:14px 18px;background:#f2f6f6;border:1px solid #dfE8e7;border-radius:7px;overflow:auto;font-size:13px;line-height:1.6}
    pre{margin:0}.output{margin:18px 0;max-width:100%;overflow:auto}.plot{max-width:100%;height:auto}.text-output{margin:14px 0;font-size:14px;white-space:pre-wrap}a{color:#087f73}
    """+formatter.get_style_defs(".code")
    return '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Recorded MD notebook</title><style>'+css+'</style></head><body><header>FYS2160 SIMPLE MD · EXECUTED NOTEBOOK · BROWSER PLAYBACK ONLY</header><main>'+"\n".join(cells)+"</main></body></html>"

def draw_gallery():
    """Pictures for the drawing section of the API page, made with the current code."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from fys2160_md import FCC, IceIh, Random, RigidWater
    target=ROOT/"docs/images";target.mkdir(parents=True,exist_ok=True)
    crystal=FCC(N=256, rho=0.8)
    fig, axes=plt.subplots(2, 3, figsize=(12, 8))
    views=[("crystal.draw()", {}), ('view="front", projection="orthographic"', dict(view="front", projection="orthographic")),
           ('view="side"', dict(view="side")), ('view="top", projection="orthographic"', dict(view="top", projection="orthographic")),
           ("view=(30, -30)", dict(view=(30, -30))), ("view=(10, 200)", dict(view=(10, 200)))]
    for ax, (title, options) in zip(axes.flat, views):
        crystal.draw(ax=ax, **options);ax.set_title(title, fontsize=10)
    fig.tight_layout();fig.savefig(target/"draw_views.png", dpi=90);plt.close(fig)
    ice=IceIh(N=768, density=0.92, molecule="water", bond=RigidWater(oh_length=1.0, hoh_angle=109.47), seed=1)
    Lx, Ly, Lz=ice.box.lengths
    fig, axes=plt.subplots(1, 3, figsize=(13, 4.2), gridspec_kw=dict(width_ratios=[1, 1, 1.6]))
    Random(N=300, box=(40, 40), boundary="walls", min_distance=1.0, seed=1).draw(ax=axes[0])
    axes[0].set_title("flat, with walls", fontsize=10)
    Random(species={"A": 60, "B": 40}, rho=0.05, min_distance=1.0, seed=2).draw(ax=axes[1], radius={"A": 0.5, "B": 0.7})
    axes[1].set_title('a mixture, radius={"A": 0.5, "B": 0.7}', fontsize=10)
    ice.with_box([Lx, Ly, 4*Lz]).draw(ax=axes[2])
    axes[2].set_title("a water slab: z lies across", fontsize=10)
    fig.tight_layout();fig.savefig(target/"draw_systems.png", dpi=90);plt.close(fig)
    print("Drawing gallery made.")

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--execute",action="store_true");args=parser.parse_args()
    for name in ("quickstart", "isobar", "gravity", "water"):
        path=ROOT/"examples"/(name+".ipynb");nb=nbformat.read(path,as_version=4)
        if args.execute:
            import matplotlib.font_manager  # Build font cache before recording cells.
            with tempfile.TemporaryDirectory(prefix="fys2160-docs-") as work:
                NotebookClient(nb,timeout=600,kernel_name="python3",resources={"metadata":{"path":work}}).execute()
            nbformat.write(nb,path)
        if not any(c.get("outputs") for c in nb.cells):
            raise RuntimeError(f"Run with --execute first to create {name}'s recorded outputs.")
        target=ROOT/"docs/examples/rendered";target.mkdir(parents=True,exist_ok=True)
        (target/(name+".html")).write_text(export(nb))
        downloads=ROOT/"docs/downloads";downloads.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(path,downloads/path.name)
        print(f"Static {name} notebook exported with embedded trajectories; no kernel required.")
    draw_gallery()

if __name__=="__main__":main()
