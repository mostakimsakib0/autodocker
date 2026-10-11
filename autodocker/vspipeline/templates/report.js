/* theme */
(function(){ var t=document.getElementById('theme-toggle');
  var cur=localStorage.getItem('vs-theme');
  if(cur){ document.documentElement.setAttribute('data-theme',cur); }
  if(t){ t.onclick=function(){ var n=document.documentElement.getAttribute('data-theme')==='dark'?'light':'dark';
    document.documentElement.setAttribute('data-theme',n); localStorage.setItem('vs-theme',n); }; }
})();
/* table */
(function(){ var s=document.getElementById('results-search');
  if(s){ s.oninput=function(){ var q=s.value.toLowerCase();
    document.querySelectorAll('#all-results tbody tr').forEach(function(tr){
      tr.style.display = tr.textContent.toLowerCase().indexOf(q)>=0 ? '' : 'none'; }); }; }
  document.querySelectorAll('#all-results th.sortable').forEach(function(th){
    th.onclick=function(){ var tb=th.closest('table').querySelector('tbody');
      var i=Array.prototype.indexOf.call(th.parentNode.children,th);
      var asc=!th.classList.contains('sort-asc');
      th.parentNode.querySelectorAll('th').forEach(function(c){c.classList.remove('sort-asc','sort-desc');});
      th.classList.add(asc?'sort-asc':'sort-desc');
      Array.prototype.slice.call(tb.children).sort(function(a,b){
        var x=a.children[i].textContent.trim(), y=b.children[i].textContent.trim();
        var nx=parseFloat(x), ny=parseFloat(y);
        if(!isNaN(nx)&&!isNaN(ny)) return asc?nx-ny:ny-nx;
        return asc? x.localeCompare(y): y.localeCompare(x); }).forEach(function(r){tb.appendChild(r);});
    }; });
})();
/* chart */
(function(){ var tip=document.getElementById('chart-tip');
  if(!tip){ tip=document.createElement('div'); tip.id='chart-tip'; tip.className='chart-tip';
    document.body.appendChild(tip); }
  function bind(){ document.querySelectorAll('[data-tip]').forEach(function(el){
    el.addEventListener('mousemove',function(e){ tip.style.display='block';
      tip.textContent=el.getAttribute('data-tip'); tip.style.left=(e.clientX+12)+'px';
      tip.style.top=(e.clientY+12)+'px'; });
    el.addEventListener('mouseleave',function(){ tip.style.display='none'; }); }); }
  bind();
  window.addEventListener('load',function(){ if(window.selectHit){
    document.querySelectorAll('.scatter-pt').forEach(function(pt){
      pt.addEventListener('click',function(){ window.selectHit(pt.getAttribute('data-name')); }); }); } });
})();
