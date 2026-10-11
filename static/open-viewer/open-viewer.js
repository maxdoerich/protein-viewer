/* Open-a-structure viewer: any PDB entry (by ID) or local .pdb / .cif / .mmcif / .bcif file, shown the way the other
 * viewers show a protein: same stage, text box, legend and buttons (gv-viewer.css, loaded first), the family's cartoon
 * colours, and the developer mode (static/dev/dev-mode.js), whose edits are kept per structure.
 *
 * Mounted onto every [data-open-viewer] element. Mol* 5.13 is vendored in ../molstar/. Only the polymer is drawn
 * (protein and nucleic acid as cartoon); water and ligands are left out. The last PDB ID fetched is remembered and
 * fetched again on the next visit, so the developer mode's edits (text, notes, colours) have their structure back;
 * an uploaded file is not kept.
 */
(function () {
  'use strict';

  var HERE = (document.currentScript && document.currentScript.src || '').replace(/[^/]*$/, '');
  var MOLSTAR = HERE + '../molstar/';
  var PDB_ID = /^[0-9][A-Za-z0-9]{3}$/;
  var FORMAT = { pdb: 'pdb', ent: 'pdb', cif: 'mmcif', mmcif: 'mmcif', bcif: 'mmcif' };
  var LAST = 'pv-open:last';
  var LEGEND = [['--gv-helix', 'α-helix'], ['--gv-strand', 'β-strand'], ['--gv-loop', 'loop']];

  var ICON = {
    reset: '<path d="M4 12a8 8 0 1 0 2.4-5.7" /><path d="M4 4v4.5h4.5" />',
    spin: '<ellipse cx="12" cy="12" rx="9" ry="4" /><path d="M12 3v18" stroke-dasharray="2 2.5" />' +
          '<path d="m16.5 14.6 2.6 1-.8 2.7" />',
    expand: '<path d="M4 9V4h5M15 4h5v5M20 15v5h-5M9 20H4v-5" />',
    shrink: '<path d="M9 4v5H4M20 9h-5V4M15 20v-5h5M4 15h5v5" />'
  };
  function icon(name) {
    return '<svg class="gvv-icon" viewBox="0 0 24 24" aria-hidden="true" focusable="false" fill="none" ' +
           'stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round">' +
           ICON[name] + '</svg>';
  }
  function el(tag, cls, html) {
    var e = document.createElement(tag);
    if (cls) e.className = cls;
    if (html != null) e.innerHTML = html;
    return e;
  }
  function esc(s) { return String(s).replace(/[&<>"]/g, function (c) { return '&#' + c.charCodeAt(0) + ';'; }); }
  function stillMotion() {
    var m = document.documentElement.getAttribute('data-motion') || 'auto';
    return m === 'reduced' || (m !== 'full' && matchMedia('(prefers-reduced-motion: reduce)').matches);
  }
  var colourCtx = document.createElement('canvas').getContext('2d');
  /** a --gv-* custom property (resolved to its token) as 0xRRGGBB */
  function colour(host, name) {
    var v = getComputedStyle(host).getPropertyValue(name).trim();
    colourCtx.fillStyle = '#000';
    colourCtx.fillStyle = v;
    var s = colourCtx.fillStyle;
    if (s.charAt(0) === '#') return parseInt(s.slice(1, 7), 16);
    var m = s.match(/[\d.]+/g) || [0, 0, 0];
    return (Math.round(m[0]) << 16) | (Math.round(m[1]) << 8) | Math.round(m[2]);
  }
  function loadScript(src) {
    return new Promise(function (ok, fail) {
      if (window.molstar) return ok();
      var s = document.createElement('script');
      s.src = src; s.async = true;
      s.onload = function () { ok(); };
      s.onerror = function () { fail(new Error('could not load ' + src)); };
      document.head.appendChild(s);
    });
  }
  function loadStyle(href) {
    if (document.querySelector('link[data-molstar]')) return;
    var l = document.createElement('link');
    l.rel = 'stylesheet'; l.href = href; l.setAttribute('data-molstar', '');
    document.head.appendChild(l);
  }
  function webglOk() {
    try {
      var c = document.createElement('canvas');
      return !!(c.getContext('webgl2') || c.getContext('webgl'));
    } catch (e) { return false; }
  }

  function mount(host) {
    var state = { spin: false };

    /* ---- markup: the family's, plus the bar to fetch / upload. The bar goes before the viewer, not in it: the text
       box is positioned against the viewer's top-left, so anything above the stage would sit under it. ---- */
    host.innerHTML = '';
    host.classList.add('gvv', 'ov');
    var bar = el('form', 'ov-bar',
      '<input class="ov-id" placeholder="PDB ID, e.g. 1CRN" maxlength="4" size="10" autocomplete="off" spellcheck="false" aria-label="PDB ID">' +
      '<button class="ov-btn ov-btn--go">Fetch</button>' +
      '<label class="ov-btn">Upload file<input type="file" accept=".pdb,.ent,.cif,.mmcif,.bcif" hidden></label>' +
      '<span class="ov-msg" role="status"></span>');
    var stage = el('div', 'gvv-stage');
    var canvas = el('div', 'gvv-canvas');
    var card = el('div', 'gvv-card');
    card.setAttribute('aria-live', 'polite');
    var tools = el('div', 'gvv-tools');
    var chips = el('div', 'gvv-chips');
    var legend = el('ul', 'gvv-legend');
    var hint = el('div', 'gvv-hint', 'Use <kbd>Ctrl</kbd> + scroll to zoom');
    var veil = el('div', 'gvv-veil',
      '<span class="gvv-pulse" aria-hidden="true"></span><span class="gvv-veil-text">Preparing the 3D model…</span>');
    stage.append(canvas, tools, chips, legend, hint, veil);
    host.append(stage, card);
    host.before(bar);
    var idInput = bar.querySelector('.ov-id'), fileInput = bar.querySelector('input[type=file]'),
        msg = bar.querySelector('.ov-msg');

    function toolButton(name, label, fn) {
      var b = el('button', 'gvv-tool', icon(name));
      b.type = 'button';
      b.setAttribute('aria-label', label);
      b.setAttribute('data-tool', name);
      b.title = label;
      b.addEventListener('click', fn);
      tools.appendChild(b);
      return b;
    }
    toolButton('reset', 'Reset view', function () { frame(700); });
    var bSpin = toolButton('spin', 'Rotate', function () { setSpin(!state.spin); });
    bSpin.setAttribute('aria-pressed', 'false');
    var bFull = null;
    if (host.requestFullscreen) {
      bFull = toolButton('expand', 'Full screen', function () {
        if (document.fullscreenElement === host) document.exitFullscreen();
        else host.requestFullscreen().catch(function () {});
      });
    }

    var hintTimer;
    stage.addEventListener('wheel', function (e) {
      if (e.ctrlKey || e.metaKey || document.fullscreenElement === host) return;
      e.stopPropagation();
      hint.classList.add('is-on');
      clearTimeout(hintTimer);
      hintTimer = setTimeout(function () { hint.classList.remove('is-on'); }, 1400);
    }, { capture: true, passive: true });
    if (/Mac|iPhone|iPad/.test(navigator.platform)) hint.innerHTML = 'Use <kbd>⌘</kbd> + scroll to zoom';

    /* ---- the text box and legend: what dev-mode edits (it keys an edit by the text painted here) ---- */
    function paint(eyebrow, title, text, keys) {
      card.innerHTML = '<span class="gvv-eyebrow">' + esc(eyebrow) + '</span><strong class="gvv-title">' + esc(title) +
        '</strong><span class="gvv-text">' + esc(text) + '</span>';
      legend.innerHTML = '';
      (keys || []).forEach(function (k) {
        legend.appendChild(el('li', 'gvv-key',
          '<span class="gvv-swatch" style="--gvv-swatch: var(' + k[0] + ')"></span><span class="gvv-key-label">' + k[1] + '</span>'));
      });
    }
    function paintEmpty() {
      paint('Open a structure', 'No structure yet',
        'Enter a PDB ID above and press Fetch, or upload a .pdb or .cif file. It is drawn as a cartoon in these colours; ' +
        'turn on Developer mode to change its text, colours and display.', LEGEND);
    }
    function say(t, err) { msg.textContent = t; msg.classList.toggle('is-error', !!err); }
    paintEmpty();

    function fail(t) {
      veil.classList.add('is-error');
      veil.querySelector('.gvv-veil-text').textContent = t;
      veil.hidden = false;
    }
    if (!webglOk()) return fail('This 3D figure needs WebGL, which this browser has switched off.');

    /* ---- Mol*, set up as the other viewers do ---- */
    var viewer, plugin, lib, COLORS, built = {}, focus = null;

    var booted, inView = true;
    /** Mol* is fetched (5 MB) on first need, as in the other viewers: when the figure nears the viewport */
    function boot() { return booted || (booted = bootMolstar()); }
    /** no drawing while the figure is scrolled away: a spinning model would otherwise keep the GPU busy */
    function applyView() { if (plugin) { if (inView) plugin.canvas3d.resume(); else plugin.canvas3d.pause(true); } }

    async function bootMolstar() {
      loadStyle(MOLSTAR + 'molstar.css');
      await loadScript(MOLSTAR + 'molstar.js');
      lib = window.molstar.lib;
      viewer = await window.molstar.Viewer.create(canvas, {
        extensions: [], layoutIsExpanded: false, layoutShowControls: false, layoutShowRemoteState: false,
        layoutShowSequence: false, layoutShowLog: false, layoutShowLeftPanel: false, collapseLeftPanel: true,
        viewportShowExpand: false, viewportShowControls: false, viewportShowSettings: false,
        viewportShowSelectionMode: false, viewportShowAnimation: false, viewportShowTrajectoryControls: false,
        viewportShowReset: false, viewportShowScreenshotControls: false, viewportShowToggleFullscreen: false,
        viewportFocusBehavior: 'disabled', volumeStreamingDisabled: true
      });
      plugin = viewer.plugin;
      COLORS = {};
      ['--gv-helix', '--gv-strand', '--gv-loop', '--gv-outline', '--gv-highlight', '--gv-bg']
        .forEach(function (n) { COLORS[n] = colour(host, n); });
      styleCanvas();
      veil.hidden = true;
      host.classList.add('is-ready');
      host.__pv = { plugin: plugin, lib: lib, built: built };       // what static/dev/dev-mode.js works on
      host.dispatchEvent(new CustomEvent('pv:ready'));
      if (!inView) applyView();

      if (window.ResizeObserver) new ResizeObserver(function () { viewer.handleResize(); }).observe(stage);
      document.addEventListener('fullscreenchange', function () {
        var full = document.fullscreenElement === host;
        if (bFull) {
          bFull.innerHTML = icon(full ? 'shrink' : 'expand');
          bFull.setAttribute('aria-label', full ? 'Leave full screen' : 'Full screen');
          bFull.title = bFull.getAttribute('aria-label');
        }
        setTimeout(function () { viewer.handleResize(); frame(0); }, 60);
      });
    }

    function styleCanvas() {
      var p = plugin.canvas3d.props;
      plugin.canvas3d.setProps({
        transparentBackground: true,
        checkeredTransparentBackground: false,
        camera: Object.assign({}, p.camera, {
          manualReset: true,
          helper: Object.assign({}, p.camera.helper, { axes: { name: 'off', params: {} } }) }),
        renderer: Object.assign({}, p.renderer, {
          backgroundColor: COLORS['--gv-bg'],
          highlightColor: COLORS['--gv-highlight'],
          selectColor: COLORS['--gv-highlight'],
          highlightStrength: 0.45,
          celSteps: 4,
          ambientIntensity: 0.55,
          exposure: 1.05
        }),
        postprocessing: Object.assign({}, p.postprocessing, {
          occlusion: { name: 'on', params: {
            samples: 32, multiScale: { name: 'off', params: {} }, radius: 5, bias: 0.8, blurKernelSize: 15,
            blurDepthBias: 0.5, resolutionScale: 1, color: COLORS['--gv-outline'], transparentThreshold: 0.4 } },
          outline: { name: 'on', params: {
            scale: 1, threshold: 0.33, color: COLORS['--gv-outline'], includeTransparent: true } }
        }),
        cameraFog: { name: 'on', params: { intensity: 12 } }
      });
    }

    /* ---- camera: the structure's bounding sphere, seen along -z, clear of the text box ---- */
    function frame(ms) {
      if (!plugin || !focus) return;
      var Vec3 = lib.math.LinearAlgebra.Vec3, cam = plugin.canvas3d.camera;
      var dir = Vec3.create(0, 0, -1), up = Vec3.create(0, 1, 0), right = Vec3.cross(Vec3(), dir, up);
      var target = Vec3.clone(focus.center), radius = focus.radius;
      var aspect = stage.clientWidth / Math.max(stage.clientHeight, 1);
      if (aspect < 1) radius /= aspect;                       // portrait: the width limits
      var over = getComputedStyle(card).position === 'absolute';   // from 768 up the text box sits over the scene
      radius *= over ? 1.3 : 1.1;
      if (over) {
        Vec3.scaleAndAdd(target, target, right, -0.3 * radius);
        Vec3.scaleAndAdd(target, target, up, 0.08 * radius);
      }
      var snap = cam.getInvariantFocus(target, radius, up, dir);
      snap.radiusMax = focus.radius * 4 + 50;
      cam.setState(snap, stillMotion() ? 0 : ms);
      plugin.canvas3d.requestDraw();
    }

    function setSpin(on) {
      state.spin = on;
      bSpin.setAttribute('aria-pressed', on ? 'true' : 'false');
      if (!plugin) return;
      var tb = plugin.canvas3d.props.trackball;
      plugin.canvas3d.setProps({ trackball: Object.assign({}, tb, {
        animate: on ? { name: 'spin', params: { speed: 0.25, axis: [0, -1, 0] } } : { name: 'off', params: {} } }) });
    }

    /* ---- one structure at a time: the polymer as cartoon, coloured by secondary structure as GvpA is ---- */
    var LOOK = { celShaded: true, material: { metalness: 0, roughness: 1, bumpiness: 0 } };
    function cartoon() {
      var c = COLORS, ss = { alphaHelix: c['--gv-helix'], threeTenHelix: c['--gv-helix'], piHelix: c['--gv-helix'],
        betaStrand: c['--gv-strand'], betaTurn: c['--gv-loop'], coil: c['--gv-loop'], bend: c['--gv-loop'],
        turn: c['--gv-loop'], dna: c['--gv-loop'], rna: c['--gv-loop'], carbohydrate: c['--gv-loop'] };
      return { type: 'cartoon', typeParams: Object.assign({ sizeFactor: 0.32 }, LOOK),
        color: 'secondary-structure', colorParams: { saturation: 0, lightness: 0, colors: { name: 'custom', params: ss } } };
    }

    /** `source` gives { data, format } once Mol* is up: a download cell, or the file's own bytes */
    async function show(what, source) {
      if (show.busy) return;
      show.busy = true;
      say('Loading ' + what + '…');
      veil.querySelector('.gvv-veil-text').textContent = 'Loading ' + what + '…';
      veil.classList.add('is-soft');
      veil.hidden = false;
      host.classList.add('is-busy');
      try {
        await boot();
        await plugin.clear();
        Object.keys(built).forEach(function (k) { delete built[k]; });
        focus = null;
        var s = await source();
        var traj = await plugin.builders.structure.parseTrajectory(s.data, s.format);
        var model = await plugin.builders.structure.createModel(traj);
        var struct = await plugin.builders.structure.createStructure(model);
        var poly = await plugin.builders.structure.tryCreateComponentStatic(struct, 'polymer');
        if (!poly) throw new Error('no protein or nucleic-acid chains in ' + what);
        built.polymer = (await plugin.builders.structure.representation.addRepresentation(poly, cartoon())).ref;
        var sphere = poly.data.boundary.sphere;
        focus = { center: sphere.center, radius: Math.max(sphere.radius, 5) };
        frame(0);
        var m = struct.data.model;
        paint(what, m.entryId || m.label || what,
          struct.data.polymerResidueCount + ' residues, ' + struct.data.elementCount + ' atoms. Drag to turn it.', LEGEND);
        say(what);
      } catch (e) {
        console.error('[open-viewer]', e);
        say('Could not load ' + what + '. Check the ID, or that the file is a PDB / mmCIF file.', true);
        paintEmpty();
      } finally {
        show.busy = false;
        veil.hidden = true; veil.classList.remove('is-soft');
        host.classList.remove('is-busy');
      }
    }

    function fetchId(id) {
      id = id.toLowerCase();
      show('PDB ' + id.toUpperCase(), async function () {
        var data = await plugin.builders.data.download(
          { url: 'https://www.ebi.ac.uk/pdbe/entry-files/download/' + id + '.bcif', isBinary: true, label: id.toUpperCase() },
          { state: { isGhost: true } });
        try { localStorage.setItem(LAST, id); } catch (e) { /* private window: no memory, no harm */ }
        return { data: data, format: 'mmcif' };
      });
    }
    function openFile(f) {
      var ext = f.name.split('.').pop().toLowerCase();
      show(f.name, async function () {
        if (!FORMAT[ext]) throw new Error('unsupported extension .' + ext);
        var bin = ext === 'bcif';
        var data = await plugin.builders.data.rawData(
          { data: bin ? new Uint8Array(await f.arrayBuffer()) : await f.text(), label: f.name }, { state: { isGhost: true } });
        return { data: data, format: FORMAT[ext] };
      });
    }

    bar.addEventListener('submit', function (e) {
      e.preventDefault();
      var v = idInput.value.trim();
      if (!PDB_ID.test(v)) return say('A PDB ID is four characters: a digit, then three letters or digits.', true);
      fetchId(v);
    });
    fileInput.addEventListener('change', function () {
      var f = fileInput.files[0];
      fileInput.value = '';                                   // so picking the same file again still fires
      if (f) openFile(f);
    });

    var last = null, started = false;
    try { last = localStorage.getItem(LAST); } catch (e) { /* none */ }
    if (last && PDB_ID.test(last)) idInput.value = last.toUpperCase();
    function start() {
      if (started) return;
      started = true;
      boot().then(function () { if (last && PDB_ID.test(last)) fetchId(last); }, function (err) {
        console.error('[open-viewer]', err);
        fail('The 3D viewer could not be loaded. Reload the page to try again.');
      });
    }
    if ('IntersectionObserver' in window) {
      new IntersectionObserver(function (es) {
        inView = es[es.length - 1].isIntersecting;
        if (inView) start();
        applyView();
      }, { rootMargin: '600px 0px' }).observe(host);
    } else start();
  }

  function init() { Array.prototype.forEach.call(document.querySelectorAll('[data-open-viewer]'), mount); }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})();
