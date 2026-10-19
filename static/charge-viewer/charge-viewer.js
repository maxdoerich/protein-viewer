/* GvpA + GvpC charge viewer: four GvpA and one GvpC repeat of PDB 8GBS (Anabaena flos-aquae), their surfaces coloured
 * by Coulombic potential, the same scene as the vault's ChimeraX image. One chip: "Remove GvpC" / "Add GvpC" takes
 * GvpC off the wall and puts it back.
 *
 * Mounted onto every [data-charge-viewer] element; layers on gv-viewer.css (loaded first), so stage, text box, legend,
 * buttons and developer mode (static/dev/dev-mode.js) are the family's. Mol* 5.13 is vendored in ../molstar/.
 * Data: gvpc-charge.cif (B-factor = the potential at each atom's surface, kcal/(mol·e)) and charge-viewer.json (net
 * charge, charged residues, camera), both written by proteins/gvpc_charge/tools/prep_viewer.py. data-base="…/" loads
 * them from another folder.
 */
(function () {
  'use strict';

  var HERE = (document.currentScript && document.currentScript.src || '').replace(/[^/]*$/, '');
  var MOLSTAR = HERE + '../molstar/';

  var ICON = {
    reset: '<path d="M4 12a8 8 0 1 0 2.4-5.7" /><path d="M4 4v4.5h4.5" />',
    spin: '<ellipse cx="12" cy="12" rx="9" ry="4" /><path d="M12 3v18" stroke-dasharray="2 2.5" />' +
          '<path d="m16.5 14.6 2.6 1-.8 2.7" />',
    expand: '<path d="M4 9V4h5M15 4h5v5M20 15v5h-5M9 20H4v-5" />',
    shrink: '<path d="M9 4v5H4M20 9h-5V4M15 20v-5h5M4 15h5v5" />',
    add: '<circle cx="12" cy="12" r="8.5" /><path d="M12 8v8M8 12h8" />',
    remove: '<circle cx="12" cy="12" r="8.5" /><path d="M8 12h8" />'
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
  function signed(n) { return n > 0 ? '+' + n : n < 0 ? '−' + Math.abs(n) : '0'; }

  /* the parts of the scene, as MolScript (chains as in gvpc-charge.cif: GvpA A1-A4, GvpC C). Names go in bare:
     a quoted "C" matches nothing */
  var GVPA = '(sel.atom.atom-groups :chain-test (set.has (set A1 A2 A3 A4) atom.auth_asym_id))';
  var GVPC = '(sel.atom.atom-groups :chain-test (= atom.auth_asym_id C))';

  function mount(host) {
    var base = host.getAttribute('data-base') || HERE;
    var state = { spin: false, gvpc: true };

    host.innerHTML = '';
    host.classList.add('gvv', 'cv');
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
    var method = el('p', 'cv-method');
    stage.append(canvas, tools, chips, legend, hint, veil);
    host.append(stage, card, method);

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

    var cGvpc = el('button', 'gvv-chip cv-chip-gvpc');
    cGvpc.type = 'button';
    cGvpc.addEventListener('click', function () { setGvpc(!state.gvpc); });
    chips.appendChild(cGvpc);

    var hintTimer;
    stage.addEventListener('wheel', function (e) {
      if (e.ctrlKey || e.metaKey || document.fullscreenElement === host || host.hasAttribute('data-wheel-zoom')) return;
      e.stopPropagation();
      hint.classList.add('is-on');
      clearTimeout(hintTimer);
      hintTimer = setTimeout(function () { hint.classList.remove('is-on'); }, 1400);
    }, { capture: true, passive: true });
    if (/Mac|iPhone|iPad/.test(navigator.platform)) hint.innerHTML = 'Use <kbd>⌘</kbd> + scroll to zoom';

    function fail(t) {
      veil.classList.add('is-error');
      veil.querySelector('.gvv-veil-text').textContent = t;
      veil.hidden = false;
    }
    if (!webglOk()) return fail('This 3D figure needs WebGL, which this browser has switched off.');

    var viewer, plugin, lib, COLORS, manifest, built = {};

    /* ---- the text box, legend and chips say what is on screen ---- */
    function paint() {
      var g = manifest && manifest.gvpc;
      var pos = g ? g.charged.filter(function (c) { return c.q > 0; }).map(function (c) { return c.res; }) : [];
      var neg = g ? g.charged.filter(function (c) { return c.q < 0; }).map(function (c) { return c.res; }) : [];
      var title, text;
      if (state.gvpc) {
        title = 'GvpC on four GvpA';
        text = 'One 33-residue GvpC repeat (residues ' + (g ? g.residues.join('–') : '85–117') + ') bound to four GvpA ' +
          'of the gas vesicle wall, seen from outside. Each surface is coloured by its own electrostatic potential.' +
          (g ? ' <strong>GvpC net charge: ' + signed(g.net) + '</strong> (' + pos.length + ' positive: ' + pos.join(', ') +
               '; ' + neg.length + ' negative: ' + neg.join(', ') + ').' : '');
      } else {
        title = 'Four GvpA without GvpC';
        text = 'The same four GvpA with GvpC taken off: the face GvpC binds to, coloured by electrostatic potential. ' +
          'Add GvpC to put it back.';
      }
      card.innerHTML = '<span class="gvv-eyebrow">PDB 8GBS · Anabaena flos-aquae</span>' +
        '<strong class="gvv-title">' + title + '</strong><span class="gvv-text">' + text + '</span>';

      legend.innerHTML = '<li class="gvv-key gvv-key--scale"><span class="gvv-key-label">−10</span>' +
        '<span class="gvv-swatch gvv-swatch--scale cv-scale"></span><span class="gvv-key-label">+10</span></li>' +
        '<li class="cv-unit">Coulombic potential, kcal/(mol·e)</li>';

      cGvpc.innerHTML = icon(state.gvpc ? 'remove' : 'add') + '<span>' + (state.gvpc ? 'Remove GvpC' : 'Add GvpC') + '</span>';
      cGvpc.setAttribute('aria-pressed', state.gvpc ? 'true' : 'false');

      method.innerHTML = 'GvpA and the GvpC backbone: PDB 8GBS (Dutka et al. 2023). 8GBS deposits GvpC without side ' +
        'chains; here they are taken from the AlphaFold DB model AF-P09413-F1 and placed residue by residue on the ' +
        'deposited backbone' + (g ? ' (Cα RMSD ' + g.af_rmsd_ca.toFixed(2) + ' Å' + (g.rotamers_replaced.length ?
          '; ' + g.rotamers_replaced.join(', ') + ' as library rotamers where the AlphaFold one ran into GvpA' : '') + ')' : '') +
        ', so GvpC’s charge distribution is a model. Potential: ChimeraX coulombic (AMBER ff14SB charges, ε = 4r, 1.4 Å ' +
        'off the surface); GvpA from the four GvpA, GvpC from GvpC alone.';
    }

    /* ---- Mol*, set up as the other viewers do ---- */
    var booted, inView = true;
    function boot() { return booted || (booted = bootMolstar()); }
    function applyView() { if (plugin) { if (inView) plugin.canvas3d.resume(); else plugin.canvas3d.pause(true); } }

    async function bootMolstar() {
      loadStyle(MOLSTAR + 'molstar.css');
      var mf = fetch(base + 'charge-viewer.json').then(function (r) {
        if (!r.ok) throw new Error('charge-viewer.json: ' + r.status);
        return r.json();
      });
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
      ['--cv-negative', '--cv-zero', '--cv-positive', '--gv-outline', '--gv-highlight', '--gv-bg']
        .forEach(function (n) { COLORS[n] = colour(host, n); });
      styleCanvas();
      manifest = await mf;
      await build();
      paint();
      frame(0);
      veil.hidden = true;
      host.classList.add('is-ready');
      host.__pv = { plugin: plugin, lib: lib, built: built };   // what static/dev/dev-mode.js works on
      host.dispatchEvent(new CustomEvent('pv:ready'));
      if (!inView) applyView();

      if (window.ResizeObserver) new ResizeObserver(function () { viewer.handleResize(); frame(0); }).observe(stage);
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

    /* ---- the scene: GvpA surface, GvpC surface; built once, then shown or hidden ---- */
    var LOOK = { celShaded: true, material: { metalness: 0, roughness: 1, bumpiness: 0 } };
    function potential(alpha) {
      return { type: 'molecular-surface',
        typeParams: Object.assign({ quality: 'high', probeRadius: 1.4, alpha: alpha == null ? 1 : alpha }, LOOK),
        color: 'uncertainty', colorParams: { domain: manifest.range,
          // Mol*'s uncertainty theme runs its scale backwards (low value = last colour), hence positive first
          list: { kind: 'interpolate', colors: [COLORS['--cv-positive'], COLORS['--cv-zero'], COLORS['--cv-negative']] } } };
    }
    var structRef;
    async function component(key, expr, props, hidden) {
      var comp = await plugin.builders.structure.tryCreateComponent(structRef,
        { type: { name: 'script', params: { language: 'mol-script', expression: expr } }, nullIfEmpty: true, label: key }, key);
      if (!comp) throw new Error('nothing in the structure for ' + key);
      var repr = await plugin.builders.structure.representation.addRepresentation(comp, props,
        { initialState: { isHidden: !!hidden } });
      built[key] = repr.ref;
      return comp;
    }
    async function build() {
      var data = await plugin.builders.data.download({ url: base + 'gvpc-charge.cif', label: 'GvpA + GvpC (8GBS)' },
        { state: { isGhost: true } });
      var traj = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
      var model = await plugin.builders.structure.createModel(traj);
      var struct = await plugin.builders.structure.createStructure(model);
      structRef = struct.ref;
      await component('gvpa', GVPA, potential());
      await component('gvpc', GVPC, potential());
    }

    function setVisible(name, on) {
      if (!built[name]) return;
      var cell = plugin.state.data.cells.get(built[name]);
      if (cell && !!cell.state.isHidden === on) plugin.state.data.updateCellState(built[name], { isHidden: !on });
    }
    var alphaNow = { gvpc: 1 };
    function setAlpha(name, a) {
      alphaNow[name] = a;
      return plugin.build().to(built[name]).update(function (old) { old.type.params.alpha = a; return old; }).commit();
    }
    function isShown(name) { var c = built[name] && plugin.state.data.cells.get(built[name]); return !!c && !c.state.isHidden; }
    function nextFrame() { return new Promise(function (r) { requestAnimationFrame(function () { r(); }); }); }
    /** GvpC's surface fades to `to` (never exactly 0: Mol* leaves a fully transparent part out of the bounds) */
    async function fade(to, ms) {
      var from = alphaNow.gvpc;
      if (stillMotion() || from === to) return setAlpha('gvpc', to);
      var t0 = performance.now(), k = 0;
      while (k < 1) {
        await nextFrame();
        k = Math.min(1, (performance.now() - t0) / ms);
        await setAlpha('gvpc', from + (to - from) * k * k * (3 - 2 * k));
      }
    }
    /** brings the scene to `state`, one change at a time: clicks that come quicker than a fade just queue up */
    var queue = Promise.resolve();
    function sync() {
      queue = queue.then(async function () {
        if (!plugin || !manifest) return;
        if (state.gvpc) {
          if (!isShown('gvpc')) { await setAlpha('gvpc', 0.02); setVisible('gvpc', true); }
          await fade(1, 450);
        } else if (isShown('gvpc')) {
          await fade(0.02, 450);
          setVisible('gvpc', false);
        }
      }).catch(function (e) { console.error('[charge-viewer]', e); });
      return queue;
    }
    function setGvpc(on) { state.gvpc = on; paint(); return sync(); }

    /* ---- camera: onto the wall's outer face, GvpC across the picture, clear of the text box ---- */
    function frame(ms) {
      if (!plugin || !manifest) return;
      var Vec3 = lib.math.LinearAlgebra.Vec3, cam = plugin.canvas3d.camera, v = manifest.view;
      var dir = Vec3.create(v.dir[0], v.dir[1], v.dir[2]), up = Vec3.create(v.up[0], v.up[1], v.up[2]);
      var right = Vec3.cross(Vec3(), dir, up);
      var target = Vec3.create(manifest.centre[0], manifest.centre[1], manifest.centre[2]);
      var radius = manifest.radius;
      var aspect = stage.clientWidth / Math.max(stage.clientHeight, 1);
      if (aspect < 1) radius /= aspect;
      var over = getComputedStyle(card).position === 'absolute';   // from 768 up the text box sits over the scene
      radius *= over ? 1.25 : 1.0;
      if (over) {
        Vec3.scaleAndAdd(target, target, right, -0.36 * radius);
        Vec3.scaleAndAdd(target, target, up, -0.04 * radius);
      }
      var snap = cam.getInvariantFocus(target, radius, up, dir);
      snap.radiusMax = manifest.radius * 4;
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

    paint();
    var started = false;
    function start() {
      if (started) return;
      started = true;
      boot().catch(function (err) {
        console.error('[charge-viewer]', err);
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

  function init() { Array.prototype.forEach.call(document.querySelectorAll('[data-charge-viewer]'), mount); }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})();
