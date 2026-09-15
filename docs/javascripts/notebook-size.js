// Let the static notebook participate in the page's normal scroll flow.
function fitRecordedNotebooks(){
  document.querySelectorAll('iframe.notebook-frame').forEach(frame=>{
    if(frame.dataset.autoHeight)return;
    frame.dataset.autoHeight='true';
    function attach(){
      if(frame._sizeObserver)frame._sizeObserver.disconnect();
      const body=frame.contentDocument?.body;
      if(!body)return;
      const fit=()=>{frame.style.height=(Math.ceil(body.getBoundingClientRect().height)+2)+'px';};
      frame._sizeObserver=new ResizeObserver(fit);
      frame._sizeObserver.observe(body);fit();
    }
    frame.addEventListener('load',attach);attach();
  });
}
if(document.readyState==='loading')document.addEventListener('DOMContentLoaded',fitRecordedNotebooks);
else fitRecordedNotebooks();
