"""Optional Jupyter bridge. The rendering code is bundled with the package."""
from pathlib import Path
import anywidget
import traitlets


class MDWidget(anywidget.AnyWidget):
    _esm = Path(__file__).with_name('viewer.js').read_text() + '''
function render({model,el}) {
  const view=createMDView(el,{...model.get('configuration'),frames:[model.get('frame')],live:true});
  const frameChanged=()=>view.setFrame(model.get('frame'));
  const statusChanged=()=>view.setStatus(model.get('status'));
  model.on('change:frame',frameChanged);
  model.on('change:status',statusChanged);
  statusChanged();
  return ()=>{
    model.off('change:frame',frameChanged);
    model.off('change:status',statusChanged);
    view.dispose();
  };
}
export default {render};
'''
    configuration = traitlets.Dict().tag(sync=True)
    frame = traitlets.Dict().tag(sync=True)
    status = traitlets.Unicode('running').tag(sync=True)
