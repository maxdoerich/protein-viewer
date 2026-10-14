/* Open-a-structure viewer: any PDB entry (by ID) or local .pdb / .cif / .mmcif / .bcif file, shown the way the other
 * viewers show a protein: same stage, text box, legend and buttons (gv-viewer.css, loaded first), the family's cartoon
 * colours, and the developer mode (static/dev/dev-mode.js), whose edits are kept per structure.
 *
 * Mounted onto every [data-open-viewer] element. Mol* 5.13 is vendored in ../molstar/. Only the polymer is drawn
 * (protein and nucleic acid as cartoon); water and ligands are left out. Every fetch or upload adds a structure to the
 * scene; developer mode lists them (host.__pvObjects, at the end of mount()) to select, align, rename, delete, or
 * copy a selection into a structure of its own. The PDB IDs fetched are remembered and fetched again on the next
 * visit, so the developer mode's edits (text, notes, colours) have their structures back; an uploaded file is not kept.
 * For a scene file (developer mode's Export / Import) the structures travel with it: the PDB ID, or the file's own
 * content, never a path. See host.__pvScene at the end of mount().
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
  function b64(u8) {                                         // chunked: one apply() over a whole file overflows the stack
    for (var s = '', i = 0; i < u8.length; i += 0x8000) s += String.fromCharCode.apply(null, u8.subarray(i, i + 0x8000));
    return btoa(s);
  }
  function unb64(s) { return Uint8Array.from(atob(s), function (c) { return c.charCodeAt(0); }); }
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
      '<label class="ov-btn">Upload file<input type="file" accept=".pdb,.ent,.cif,.mmcif,.bcif" multiple hidden></label>' +
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
        'Enter a PDB ID above and press Fetch, or upload a .pdb or .cif file; every one adds to the scene. They are drawn ' +
        'as cartoons in these colours; turn on Developer mode to align them, change their text, colours and display.', LEGEND);
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
    var viewer, plugin, lib, COLORS, built = {}, skey = {}, focus = null;

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
      host.__pv = { plugin: plugin, lib: lib, built: built, skey: skey };   // what static/dev/dev-mode.js works on
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

    /* ---- the structures in the scene, in the order they were added: each the polymer as cartoon, coloured by secondary
       structure as GvpA is. { id, name, what, title, text, desc, data, root, poly, repr, sphere }: id is what dev-mode
       keeps its edits under (a file name, made unique), name what its list shows and the user may change, desc how a
       scene file keeps it, data the cell to delete to remove it. built[id] is its cartoon, skey[root | poly] its id. ---- */
    var LOOK = { celShaded: true, material: { metalness: 0, roughness: 1, bumpiness: 0 } };
    var objs = [], restoring = false;
    function cartoon() {
      var c = COLORS, ss = { alphaHelix: c['--gv-helix'], threeTenHelix: c['--gv-helix'], piHelix: c['--gv-helix'],
        betaStrand: c['--gv-strand'], betaTurn: c['--gv-loop'], coil: c['--gv-loop'], bend: c['--gv-loop'],
        turn: c['--gv-loop'], dna: c['--gv-loop'], rna: c['--gv-loop'], carbohydrate: c['--gv-loop'] };
      return { type: 'cartoon', typeParams: Object.assign({ sizeFactor: 0.32 }, LOOK),
        color: 'secondary-structure', colorParams: { saturation: 0, lightness: 0, colors: { name: 'custom', params: ss } } };
    }
    function byId(id) { return objs.filter(function (o) { return o.id === id; })[0]; }
    function unique(id) {
      for (var k = id, n = 2; byId(k); n++) k = id + '#' + n;
      return k;
    }
    function forget(map) { Object.keys(map).forEach(function (k) { delete map[k]; }); }
    function reprOf(o) { return plugin.state.data.cells.get(o.repr).obj.data.repr; }
    /** where a structure sits: the representation's matrix, as the other viewers move parts (dev-mode follows it) */
    function poseOf(o) { return reprOf(o).state.transform; }
    function setPose(o, m) { var r = reprOf(o); r.setState({ transform: m }); plugin.canvas3d.update(r); }
    function validPose(p) { return Array.isArray(p) && p.length === 16 && p.every(Number.isFinite) ? p : undefined; }

    /** the PDB IDs, to fetch again on the next visit. ponytail: an uploaded or copied structure is not kept, only
        scene files carry those */
    function save() {
      if (restoring) return;
      try {
        localStorage.setItem(LAST, JSON.stringify(objs.filter(function (o) { return o.desc.pdb; }).map(function (o) {
          return { pdb: o.desc.pdb, name: o.name, pose: Array.from(poseOf(o)) };
        })));
      } catch (e) { /* private window: no memory, no harm */ }
    }
    /** the text box says what was added last */
    function paintLast() {
      var o = objs[objs.length - 1];
      if (o) paint(o.what, o.title, o.text, LEGEND); else paintEmpty();
    }
    /** the camera takes in every structure where it sits now; the memory and dev-mode's list follow */
    function changed(ms) {
      focus = null;
      if (objs.length) {
        var V = lib.math.LinearAlgebra.Vec3, c = V(), cs = objs.map(function (o) {
          return { c: V.transformMat4(V(), o.sphere.center, poseOf(o)), r: o.sphere.radius };
        });
        cs.forEach(function (s) { V.add(c, c, s.c); });
        V.scale(c, c, 1 / cs.length);
        focus = { center: c, radius: Math.max.apply(null, cs.map(function (s) { return V.distance(c, s.c) + s.r; }).concat(5)) };
      }
      frame(ms);
      save();
      host.dispatchEvent(new CustomEvent('pv:objects'));
    }

    /** `source` gives { data, format } once Mol* is up: a download cell, or the file's own bytes. What it reads is added
        to the scene. `desc` is how a scene file keeps it; `extra` is { key: the file name, name, pose }. Resolves to
        whether it worked (a failure is shown on the page, not thrown). */
    async function show(what, source, desc, extra) {
      if (show.busy) return false;
      show.busy = true;
      say('Loading ' + what + '…');
      veil.querySelector('.gvv-veil-text').textContent = 'Loading ' + what + '…';
      veil.classList.add('is-soft');
      veil.hidden = false;
      host.classList.add('is-busy');
      var made;
      try {
        await boot();
        var s = await source();
        made = s.data.ref;
        var traj = await plugin.builders.structure.parseTrajectory(s.data, s.format);
        var model = await plugin.builders.structure.createModel(traj);
        var struct = await plugin.builders.structure.createStructure(model);
        var poly = await plugin.builders.structure.tryCreateComponentStatic(struct, 'polymer');
        if (!poly) throw new Error('no protein or nucleic-acid chains in ' + what);
        var repr = await plugin.builders.structure.representation.addRepresentation(poly, cartoon());
        var m = struct.data.model, sphere = poly.data.boundary.sphere;
        var o = { id: unique(extra.key), what: what, title: m.entryId || m.label || what, desc: desc, data: made,
                  root: struct.ref, poly: poly.ref, repr: repr.ref,
                  sphere: { center: sphere.center, radius: Math.max(sphere.radius, 5) },
                  text: struct.data.polymerResidueCount + ' residues, ' + struct.data.elementCount + ' atoms. Drag to turn it.' };
        o.name = extra.name || desc.pdb || extra.key;
        built[o.id] = o.repr; skey[o.root] = skey[o.poly] = o.id;
        if (extra.pose) setPose(o, extra.pose);
        objs.push(o);
        paintLast();
        changed(0);
        say(objs.length > 1 ? what + ' added: ' + objs.length + ' structures' : what);
        return true;
      } catch (e) {
        console.error('[open-viewer]', e);
        say('Could not load ' + what + '. Check the ID, or that the file is a PDB / mmCIF file.', true);
        if (made) plugin.build().delete(made).commit().catch(function () {});     // what it had built before it failed
        paintLast();
        return false;
      } finally {
        show.busy = false;
        veil.hidden = true; veil.classList.remove('is-soft');
        host.classList.remove('is-busy');
      }
    }

    function fetchId(id, extra) {
      id = id.toLowerCase();
      return show('PDB ' + id.toUpperCase(), async function () {
        var data = await plugin.builders.data.download(
          { url: 'https://www.ebi.ac.uk/pdbe/entry-files/download/' + id + '.bcif', isBinary: true, label: id.toUpperCase() },
          { state: { isGhost: true } });
        return { data: data, format: 'mmcif' };
      }, { pdb: id.toUpperCase() }, Object.assign({ key: id + '.bcif' }, extra));
    }
    /** a file's name and content: text, or bytes for .bcif. The content is kept in the structure's desc, for a scene file */
    function openData(name, data, extra) {
      var ext = name.split('.').pop().toLowerCase(), bin = typeof data !== 'string';
      return show(name, async function () {
        if (!FORMAT[ext]) throw new Error('unsupported extension .' + ext);
        return { data: await plugin.builders.data.rawData({ data: data, label: name }, { state: { isGhost: true } }),
                 format: FORMAT[ext] };
      }, bin ? { file: name, base64: b64(data) } : { file: name, text: data }, Object.assign({ key: name }, extra));
    }
    async function openFile(f) {
      return openData(f.name, /\.bcif$/i.test(f.name) ? new Uint8Array(await f.arrayBuffer()) : await f.text());
    }

    async function remove(id) {
      var o = byId(id);
      if (!o || show.busy) return;
      show.busy = true;
      try {
        objs.splice(objs.indexOf(o), 1);
        delete built[o.id]; delete skey[o.root]; delete skey[o.poly];
        await plugin.build().delete(o.data).commit();     // the whole chain, down to its cartoon (and any surface over it)
        paintLast();
        changed(700);
      } finally { show.busy = false; }
    }
    /** the selected residues, copied into a structure of their own, one per structure they lie in. They are written
        out as mmCIF and read back in, so the copy is an uploaded file in every way: it can be moved, painted, and
        travels in a scene file. It starts where its source sits. Resolves to how many were made. */
    async function extract() {
      if (!plugin) return 0;
      var L = lib.structure.StructureElement.Loci, picks = [], n = 0;
      plugin.managers.structure.selection.entries.forEach(function (e, ref) {
        var src = objs.filter(function (o) { return o.root === ref; })[0];
        if (src && !L.isEmpty(e.selection)) picks.push([src, L.toStructure(e.selection)]);
      });
      for (var i = 0; i < picks.length; i++) {
        var src = picks[i][0], name = src.name + ' selection';
        var text = lib.structure.to_mmCIF('selection', picks[i][1], false);
        if (await openData(name + '.cif', text, { name: name, pose: Array.from(poseOf(src)) })) n++;
      }
      return n;
    }

    /* ---- scene file (developer mode: Export / Import): the structures and the camera; the edits are dev-mode's.
       { structures: [{ name, pose: Mat4, pdb | file + text | base64 }], camera }. A scene of one structure, from before
       there were several ({ structure: { pdb | file … } }), still imports. ---- */
    var plain = function (k, v) { return ArrayBuffer.isView(v) ? Array.from(v) : v; };
    function load(d) {
      var extra = { name: typeof d.name === 'string' ? d.name : undefined, pose: validPose(d.pose) };
      if (d.pdb) {
        if (!PDB_ID.test(d.pdb)) throw new Error('"' + d.pdb + '" is not a PDB ID');
        idInput.value = String(d.pdb).toUpperCase();
        return fetchId(d.pdb, extra);
      }
      if (d.file && (d.text != null || d.base64)) {
        return openData(String(d.file).replace(/^.*[\\/]/, ''), d.base64 ? unb64(d.base64) : String(d.text), extra);
      }
      throw new Error('the scene has a structure this viewer cannot read');
    }
    host.__pvScene = {
      get: function () {
        var cam = plugin && focus && plugin.canvas3d.camera.getSnapshot();
        return {
          structures: objs.length ? objs.map(function (o) {
            return Object.assign({ name: o.name, pose: Array.from(poseOf(o)) }, o.desc);
          }) : undefined,
          camera: cam ? JSON.parse(JSON.stringify(cam, plain)) : undefined
        };
      },
      set: async function (s) {
        var list = s.structures || (s.structure && [s.structure]);
        if (list) {
          if (!Array.isArray(list)) throw new Error('the scene has no list of structures');
          await boot();
          await plugin.clear();
          objs = []; forget(built); forget(skey);
          restoring = true;                                   // one save at the end, not one per structure
          try {
            for (var i = 0; i < list.length; i++) {
              if (!list[i] || typeof list[i] !== 'object') throw new Error('the scene has a structure this viewer cannot read');
              if (!(await load(list[i]))) throw new Error('a structure of this scene could not be loaded');
            }
          } finally { restoring = false; save(); }
        }
        if (s.camera && plugin && focus) {
          plugin.canvas3d.camera.setState(s.camera, 0);
          plugin.canvas3d.requestDraw();
        }
      }
    };

    /* ---- what developer mode's list of structures works with (static/dev/dev-mode.js); pv:objects says it changed ---- */
    host.__pvObjects = {
      list: function () { return objs.map(function (o) { return { id: o.id, name: o.name }; }); },
      rename: function (id, name) { var o = byId(id); if (o) { o.name = name; save(); } },
      remove: remove,
      extract: extract,
      pose: function (id, m) { var o = byId(id); if (o) { setPose(o, m); changed(700); } }
    };

    bar.addEventListener('submit', function (e) {
      e.preventDefault();
      var v = idInput.value.trim();
      if (!PDB_ID.test(v)) return say('A PDB ID is four characters: a digit, then three letters or digits.', true);
      fetchId(v);
    });
    fileInput.addEventListener('change', function () {
      var fs = Array.from(fileInput.files);
      fileInput.value = '';                                   // so picking the same file again still fires
      fs.reduce(function (p, f) { return p.then(function () { return openFile(f); }); }, Promise.resolve());
    });

    var last = [], started = false;
    try {                                                     // a list; before there were several, one ID
      var mem = localStorage.getItem(LAST);
      last = !mem ? [] : mem.charAt(0) === '[' ? JSON.parse(mem) : [{ pdb: mem }];
    } catch (e) { /* none */ }
    last = Array.isArray(last) ? last.filter(function (l) { return l && PDB_ID.test(l.pdb); }) : [];
    if (last.length) idInput.value = String(last[last.length - 1].pdb).toUpperCase();
    function start() {
      if (started) return;
      started = true;
      boot().then(async function () {
        var ok = true;
        restoring = true;                                     // one failed fetch must not drop the others from the memory
        for (var i = 0; i < last.length; i++) {
          ok = (await fetchId(last[i].pdb, { name: last[i].name, pose: validPose(last[i].pose) })) && ok;
        }
        restoring = false;
        if (ok) save();
      }, function (err) {
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
