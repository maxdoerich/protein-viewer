/* Snapshot button for every viewer on the page (they all share the .gvv markup): a PNG of the current frame,
 * LONG pixels on its long side whatever the size on screen. Load it after the viewer script.
 *
 * Mol* draws the protein at that size on a transparent background (its own screenshot helper). The scene's
 * background and text boxes (card, legend, notes) are DOM, so they are rasterised at the same size by drawing the
 * boxes, with their computed styles, into an SVG <foreignObject> image. The layers go back in the order the page
 * stacks them: background + boxes, protein, note lines. Buttons are not part of any layer.
 */
(function () {
  'use strict';

  var LONG = 3840;
  var BOXES = '.gvv-card, .gvv-legend, .gvv-side, .pv-note';
  var ICON = '<svg class="gvv-icon" viewBox="0 0 24 24" aria-hidden="true" focusable="false" fill="none" ' +
             'stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round">' +
             '<path d="M4 8h3l1.5-2.5h7L17 8h3v11H4Z" /><circle cx="12" cy="13.5" r="3.5" /></svg>';

  /** the page's web font as an @font-face with the file inlined: an image of the boxes cannot fetch it itself */
  var fontP;
  function fontCss(host) {
    return fontP || (fontP = (async function () {
      var fam = getComputedStyle(host).fontFamily.split(',')[0].replace(/["']/g, '').trim(), out = '';
      for (var sh of document.styleSheets) {
        var rules;
        try { rules = sh.cssRules; } catch (e) { continue; }       // a sheet from another origin is unreadable
        for (var r of rules) {
          if (!(r instanceof CSSFontFaceRule) || r.style.getPropertyValue('font-family').replace(/["']/g, '').trim() !== fam) continue;
          var m = /url\(["']?([^"')]+)/.exec(r.style.getPropertyValue('src'));
          if (!m) continue;
          var blob = await (await fetch(new URL(m[1], sh.href || location.href))).blob();
          var uri = await new Promise(function (ok) {
            var f = new FileReader();
            f.onload = function () { ok(f.result); };
            f.readAsDataURL(blob);
          });
          out += '@font-face{font-family:"' + fam + '";font-weight:' + r.style.getPropertyValue('font-weight') +
                 ';src:url(' + uri + ')}';
        }
      }
      return out;
    })().catch(function (e) { console.warn('[gv-snapshot] font not embedded, text falls back to the system font', e); return ''; }));
  }

  /** copy what the page computed for a, and for its descendants, onto the clone b: it will be styled by no sheet */
  function inline(a, b) {
    var s = getComputedStyle(a), t = '';
    for (var i = 0; i < s.length; i++) t += s[i] + ':' + s.getPropertyValue(s[i]) + ';';
    b.style.cssText = t;
    for (var c = a.firstElementChild, d = b.firstElementChild; c; c = c.nextElementSibling, d = d.nextElementSibling) inline(c, d);
  }

  /** the visible elements matching `sel` as an image of W x H pixels, drawn in the stage's w x h CSS pixels */
  async function raster(host, stage, sel, bg, W, H) {
    var sr = stage.getBoundingClientRect(), w = stage.clientWidth, h = stage.clientHeight;
    var root = document.createElement('div');
    root.style.cssText = 'position:relative;overflow:hidden;width:' + w + 'px;height:' + h + 'px';
    if (bg) {
      var bs = getComputedStyle(stage);
      root.style.backgroundImage = bs.backgroundImage;
      root.style.backgroundColor = bs.backgroundColor;
    }
    var n = 0;
    Array.prototype.forEach.call(host.querySelectorAll(sel), function (e) {
      var s = getComputedStyle(e), r = e.getBoundingClientRect();
      // on a phone the text box is under the scene, not on it
      if (s.position !== 'absolute' || s.display === 'none' || s.visibility === 'hidden' || !r.width || !r.height) return;
      var c = e.cloneNode(true);
      inline(e, c);
      Object.assign(c.style, { position: 'absolute', left: r.left - sr.left + 'px', top: r.top - sr.top + 'px',
        right: 'auto', bottom: 'auto', boxSizing: 'border-box', width: r.width + 'px', height: r.height + 'px',
        maxWidth: 'none', maxHeight: 'none', margin: '0', transform: 'none' });
      root.appendChild(c);
      n++;
    });
    if (!n && !bg) return null;
    var css = (await fontCss(host)).replace(/&/g, '&amp;').replace(/</g, '&lt;');
    var svg = '<svg xmlns="http://www.w3.org/2000/svg" width="' + W + '" height="' + H + '" viewBox="0 0 ' + w + ' ' + h + '">' +
      '<style>' + css + '</style><foreignObject width="' + w + '" height="' + h + '">' +
      new XMLSerializer().serializeToString(root) + '</foreignObject></svg>';
    var img = new Image();
    img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svg);
    await img.decode();
    return img;
  }

  async function snapshot(host) {
    var plugin = host.__pv.plugin, stage = host.querySelector('.gvv-stage');
    var k = LONG / Math.max(stage.clientWidth, stage.clientHeight);      // the stage's own proportions: same framing
    var W = Math.round(stage.clientWidth * k), H = Math.round(stage.clientHeight * k);

    var shot = plugin.helpers.viewportScreenshot, was = shot.values;
    shot.behaviors.values.next(Object.assign({}, was, { transparent: true, axes: { name: 'off', params: {} },
      resolution: { name: 'custom', params: { width: W, height: H } } }));
    shot.behaviors.relativeCrop.next({ x: 0, y: 0, width: 1, height: 1 });
    // the outline is a width in device pixels (scale x pixel ratio, in the image as on screen): widen it with the
    // image, or it turns into a hairline; the ratio is divided out, as the screen's is no part of this image
    var c3d = plugin.canvas3d, pp = c3d.props.postprocessing, ol = pp.outline;
    if (ol.name === 'on') {
      c3d.setProps({ postprocessing: Object.assign({}, pp, { outline: { name: 'on',
        params: Object.assign({}, ol.params, { scale: ol.params.scale * Math.min(k / c3d.webgl.pixelRatio, 5) }) } }) });
    }
    try { await shot.getImageDataUri(); } finally {
      shot.behaviors.values.next(was);
      c3d.setProps({ postprocessing: pp });
    }

    var out = document.createElement('canvas'), ctx = out.getContext('2d');
    out.width = W; out.height = H;
    ctx.drawImage(await raster(host, stage, BOXES, true, W, H), 0, 0);
    ctx.drawImage(shot.canvas, 0, 0, W, H);
    var lines = await raster(host, stage, '.pv-leaders', false, W, H);       // developer mode's note lines
    if (lines) ctx.drawImage(lines, 0, 0);

    var blob = await new Promise(function (ok) { out.toBlob(ok, 'image/png'); });
    var a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = 'protein-snapshot.png';
    a.click();
    setTimeout(function () { URL.revokeObjectURL(a.href); }, 10000);
  }

  function mount(host) {
    var tools = host.querySelector('.gvv-tools');
    if (!tools || tools.querySelector('[data-tool="snapshot"]')) return;
    var b = document.createElement('button');
    b.className = 'gvv-tool';
    b.type = 'button';
    b.innerHTML = ICON;
    b.setAttribute('data-tool', 'snapshot');
    b.setAttribute('aria-label', 'Save image');
    b.title = 'Save image (' + LONG + ' px, PNG)';
    b.disabled = !host.__pv;                                  // until Mol* is up
    host.addEventListener('pv:ready', function () { b.disabled = false; });
    b.addEventListener('click', function () {
      b.disabled = true;
      snapshot(host).catch(function (e) { console.error('[gv-snapshot]', e); })
        .then(function () { b.disabled = false; });
    });
    tools.insertBefore(b, tools.querySelector('[data-tool="expand"]'));   // before full screen, or last
  }

  function init() { Array.prototype.forEach.call(document.querySelectorAll('.gvv'), mount); }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
  else init();
})();
