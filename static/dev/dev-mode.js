/* Developer mode: edit a viewer in place.
 *
 * Works on the one [data-pv-key] viewer of index.html, whichever script mounted it: the four viewers share their markup
 * (gv-viewer.css) and each hands over its Mol* plugin as host.__pv once it is ready. Everything edited is one JSON
 * object in localStorage ("pv-config:<key>"), applied on every load: with developer mode off the page shows the
 * edited viewer and nothing else. The [data-pv-dev] button switches the editing tools on:
 *   - text box and legend: move (drag the tab), resize (drag the corner), edit the text of the current step;
 *   - notes: boxes tied to a residue, placed in 3D, with a line to it;
 *   - which of the viewer's buttons are shown;
 *   - colour and display mode of the selection (click residues in 3D, or drag over the sequence view);
 *   - the structures of the scene, in the open viewer (host.__pvObjects): select, align, rename, delete one, or make
 *     a new one from the selection;
 *   - the sequence view, where the proteins are named;
 *   - export / import of the whole scene as one JSON file (see SCENE below).
 * Text is written in a small markup, see render(): "# " label, "## " title, "[pine] " legend key, **bold**, *italic*.
 */
(function start() {
  'use strict';

  // the viewer script mounts first (index.html loads it before this one), at DOMContentLoaded at the latest
  if (document.readyState === 'loading') return document.addEventListener('DOMContentLoaded', start);
  var host = document.querySelector('[data-pv-key]');
  if (!host) return;
  var stage = host.querySelector('.gvv-stage');
  var STORE = 'pv-config:' + host.getAttribute('data-pv-key');
  var cfg = { box: {}, text: {}, hide: [], paint: [], mode: {}, notes: [], names: {}, seq: false };
  try { Object.assign(cfg, JSON.parse(localStorage.getItem(STORE) || '{}')); } catch (e) { console.warn('[dev-mode]', e); }
  function save() { localStorage.setItem(STORE, JSON.stringify(cfg)); }

  var dev = false, plugin = null, pv, S, V3, M4;
  var objs = host.__pvObjects;                          // the open viewer's list of structures; the others have none

  /* A scene file is one JSON object that rebuilds the viewer on any device: no path, no URL of ours.
       { format: SCENE, version: 1, viewer: <?v= key>, config: <cfg above>,
         structures: [{ name, pose, pdb | file + text | base64 }], camera: <Mol* snapshot> }
     config is this file's; structures and camera come from the viewer if it has any of its own (host.__pvScene, the
     "open" viewer). The other viewers carry their structure in their folder, so their scene is the config alone. */
  var SCENE = 'protein-viewer-scene';

  /* every colour the site has (tokens.css), under the palette's names: nothing else can be picked */
  var PALETTE = { cayenne: '--clr-ring-core', tangerine: '--clr-title', peach: '--clr-title-shadow', frost: '--clr-frosted-blue',
                  pine: '--clr-pine-blue', cyan: '--clr-light-cyan', grey: '--clr-muted', ink: '--clr-ink', white: '--clr-body' };
  var MORE = {                                          // the "more colors" button: more swatches for the selection (not usable in legend keys)
    Teal: 'F0FBFE DDF8F6 CAF8FB AFE2E3 8CDEE5 83CECE 49D2DF 5FB7B9 248F9A 3C8285',
    Peach: 'F9E5DA F8D8C6 F7CCA7 F7C6AC F0BB96',
    Orange: 'F5AB7F F59C69 FB946D F3894F DC6B3C F56D2D DD501B CA3A14',
    'Red and pink': 'FDF0EE F5CEC1 F4ABAC F07E80 EB595C',
    Green: 'E6F5E8 B7E1BE 72A052'
  };
  /* colours for "by chain": blue and orange first, only a few clearly different steps of each (picked for distance in Lab),
     then red / pink / green, then the steps that were left out. ponytail: more chains than this repeat the list */
  var CHAIN = ('248F9A DD501B AFE2E3 F5AB7F 49D2DF F3894F 72A052 EB595C F4ABAC B7E1BE F07E80 F7CCA7 5FB7B9 FB946D DC6B3C ' +
               'F56D2D 8CDEE5 3C8285 CA3A14 F59C69 F0BB96 83CECE F7C6AC').split(' ').map(function (h) { return parseInt(h, 16); });
  var MODES = {                                         // the representations the viewers already draw with
    cartoon: ['Cartoon', { sizeFactor: 0.3 }],
    'gaussian-surface': ['Surface', { quality: 'custom', resolution: 1.6, radiusOffset: 0.8, smoothness: 1.6 }],
    'molecular-surface': ['Molecular surface', { probeRadius: 1.4 }],
    'cartoon-surface': ['Cartoon + surface', null],     // not a Mol* type: a cartoon with a second, see-through representation over it
    spacefill: ['Spheres', { sizeFactor: 1 }]
  };
  var LOOK = { celShaded: true, material: { metalness: 0, roughness: 1, bumpiness: 0 } };   // as in the viewers

  var shades = {};
  /** a custom property as 0xRRGGBB. ponytail: the tokens are all #RRGGBB; resolve through a canvas, as the viewers do, if one ever is not */
  function hex(prop) {
    if (shades[prop] == null) shades[prop] = parseInt(getComputedStyle(host).getPropertyValue(prop).trim().slice(1), 16);
    return shades[prop];
  }
  function css(c) { return '#' + ('00000' + c.toString(16)).slice(-6); }
  function mix(a, b, k) {
    var ch = function (s) { return Math.round((a >> s & 255) + ((b >> s & 255) - (a >> s & 255)) * k) << s; };
    return ch(16) | ch(8) | ch(0);
  }
  function el(tag, cls, html) {
    var e = document.createElement(tag);
    if (cls) e.className = cls;
    if (html != null) e.innerHTML = html;
    return e;
  }
  function hash(s) {
    for (var h = 5381, i = 0; i < s.length; i++) h = (h * 33 ^ s.charCodeAt(i)) >>> 0;
    return h.toString(36);
  }

  /* ---------------------------------------------------------------- markup -- */
  function esc(s) { return s.replace(/[&<>"]/g, function (c) { return '&#' + c.charCodeAt(0) + ';'; }); }
  function inline(s) {
    return esc(s).replace(/\*\*(.+?)\*\*/g, '<strong>$1</strong>').replace(/\*(.+?)\*/g, '<em>$1</em>');
  }
  function key(colour, kind, label) {
    return '<span class="gvv-swatch' + (kind ? ' gvv-swatch--' + kind : '') + '"' +
      (colour ? ' style="--gvv-swatch: var(' + (PALETTE[colour] || '--' + colour) + ')"' : '') + '></span>' +
      '<span class="gvv-key-label">' + inline(label) + '</span>';
  }
  /** one line, one element:  # label (blue bar) · ## title (orange) · [colour] legend key · [scale] left | right · text */
  function render(text, tag) {
    return text.split('\n').map(function (line) {
      var m, cls = 'gvv-text', body;
      if (!line.trim()) return '';
      if ((m = /^\[scale\]\s*(.*?)\s*\|\s*(.*)$/.exec(line))) {
        cls = 'gvv-key gvv-key--scale';
        body = '<span class="gvv-key-label">' + inline(m[1]) + '</span>' + key('', 'scale', m[2]);
      } else if ((m = /^\[([\w-]+)(?::(\w+))?\]\s*(.*)$/.exec(line))) { cls = 'gvv-key'; body = key(m[1], m[2], m[3]); }
      else if ((m = /^##\s*(.*)$/.exec(line))) { cls = 'gvv-title'; body = inline(m[1]); }
      else if ((m = /^#\s*(.*)$/.exec(line))) { cls = 'gvv-eyebrow'; body = inline(m[1]); }
      else body = inline(line);
      return '<' + tag + ' class="' + cls + '">' + body + '</' + tag + '>';
    }).join('');
  }
  /** what a viewer painted, as markup: the text an edit starts from */
  function toMarkup(box) {
    return Array.prototype.map.call(box.children, function (c) {
      var t = c.textContent.trim(), has = function (n) { return c.classList.contains(n); };
      if (has('gvv-eyebrow')) return '# ' + t;
      if (has('gvv-title')) return '## ' + t;
      if (has('gvv-key--scale')) return '[scale] ' + c.firstChild.textContent + ' | ' + c.lastChild.textContent;
      if (has('gvv-key')) {
        var sw = c.querySelector('.gvv-swatch');
        var v = /var\(--([\w-]+)\)/.exec(sw.getAttribute('style') || ''), kind = /gvv-swatch--(\w+)/.exec(sw.className);
        return '[' + (v ? v[1] : 'pine') + (kind ? ':' + kind[1] : '') + '] ' + t;
      }
      return t;
    }).join('\n');
  }

  /* ----------------------------------------------- the text box and the legend --
     A viewer repaints them on every step. What it painted is the box's `orig`, and its hash the key the edited
     text of that step is stored under, so no viewer has to say which step it is on. */
  var items = [], boxes = {}, picked = null;

  function place(box) {
    var g = cfg.box[box.id];
    box.el.classList.toggle('pv-placed', !!g);
    ['x', 'y', 'w', 'h'].forEach(function (k) {
      if (g && g[k] != null) box.el.style.setProperty('--pv-' + k, g[k] + '%'); else box.el.style.removeProperty('--pv-' + k);
    });
  }
  function geometry(e) {                                  // where a box sits in its offset parent, px
    var p = e.offsetParent, r = e.getBoundingClientRect(), q = p.getBoundingClientRect();
    return { x: r.left - q.left - p.clientLeft, y: r.top - q.top - p.clientTop, w: r.width, h: r.height };
  }
  function setBox(box, x, y, w, h) {                      // stored in % of the parent, so it follows the figure's size
    var p = box.el.offsetParent, pc = function (v, of) { return Math.round(v / of * 1000) / 10; };
    cfg.box[box.id] = { x: pc(x, p.clientWidth), y: pc(y, p.clientHeight), w: pc(w, p.clientWidth),
                        h: h == null ? null : pc(h, p.clientHeight) };
    place(box);
  }

  function track(id, label, e) {
    if (!e || !e.isConnected) return;
    var box = boxes[id] = { id: id, label: label, el: e };
    var mo = new MutationObserver(function () { box.orig = null; box.draw(); if (picked === box) pick(box); placeNotes(); });
    box.draw = function () {
      if (box.orig == null) { box.orig = e.innerHTML; box.key = hash(box.orig); box.def = toMarkup(e); }
      var t = cfg.text[box.key];
      e.innerHTML = t == null ? box.orig : render(t, e.tagName === 'UL' ? 'li' : 'span');
      mo.takeRecords();                                   // our own change is not a repaint
    };
    box.text = function () { var t = cfg.text[box.key]; return t == null ? box.def : t; };
    box.setText = function (t) {
      if (t === box.def) delete cfg.text[box.key]; else cfg.text[box.key] = t;
      box.draw();
    };
    box.grab = function () { return Object.assign(geometry(e), { tall: (cfg.box[id] || {}).h != null }); };
    box.move = function (s, dx, dy) { setBox(box, s.x + dx, s.y + dy, s.w, s.tall ? s.h : null); };
    box.size = function (s, dx, dy) { setBox(box, s.x, s.y, Math.max(60, s.w + dx), Math.max(24, s.h + dy)); };
    mo.observe(e, { childList: true, subtree: true, characterData: true });
    box.draw();
    place(box);
    items.push(box);
  }
  track('card', 'Text box', host.querySelector('.gvv-card'));
  track('legend', 'Legend', host.querySelector('.gvv-legend'));
  function stepKey() { return boxes.card ? boxes.card.key : null; }

  /* ------------------------------------------------------------ the buttons -- */
  function applyHide() {
    Array.prototype.forEach.call(host.querySelectorAll('.gvv-tool, .gvv-chips'), function (b) {
      b.hidden = cfg.hide.indexOf(b.getAttribute('data-tool') || 'chips') >= 0;
    });
    items.forEach(function (b) { if (!b.note) b.el.hidden = cfg.hide.indexOf(b.id) >= 0; });   // the text box and the legend
  }
  applyHide();

  /* ------------------------------------------------- frames (developer mode) --
     The boxes lie behind the scene, where no pointer reaches them: each gets a frame above it, with a tab to
     drag and pick it by and a corner to resize it. */
  var frames = el('div', 'pv-frames');
  host.appendChild(frames);

  function drag(e, start, fn) {
    var t = e.currentTarget, x = e.clientX, y = e.clientY;
    e.preventDefault();
    t.setPointerCapture(e.pointerId);
    t.onpointermove = function (m) { fn(start, m.clientX - x, m.clientY - y); };
    t.onpointerup = function () { t.onpointermove = t.onpointerup = null; save(); };
  }
  function frameOf(it) {
    if (it.frame) return it.frame;
    var f = it.frame = el('div', 'pv-frame', '<button type="button" class="pv-grip"></button><span class="pv-corner"></span>');
    f.firstChild.textContent = it.label;
    f.firstChild.title = 'Drag to move, click to edit the text';
    f.firstChild.addEventListener('pointerdown', function (e) { pick(it); drag(e, it.grab(), it.move); });
    f.lastChild.addEventListener('pointerdown', function (e) { pick(it); drag(e, it.grab(), it.size); });
    frames.appendChild(f);
    return f;
  }
  var looping = false;
  function syncFrames() {
    looping = dev;
    if (!dev) return;
    var hr = host.getBoundingClientRect();
    items.forEach(function (it) {
      var f = frameOf(it), r = it.el.getBoundingClientRect();
      var pos = 'left:' + (r.left - hr.left) + 'px;top:' + (r.top - hr.top) + 'px;width:' + r.width + 'px;height:' + r.height + 'px';
      f.hidden = !it.el.offsetParent;
      if (f.style.cssText !== pos) f.style.cssText = pos;
      f.classList.toggle('is-picked', it === picked);
    });
    requestAnimationFrame(syncFrames);
  }

  /* ------------------------------------------------------------------ panel -- */
  var panel = el('aside', 'pv-panel',
    '<h2>Developer mode</h2>' +
    '<section><h3>Text</h3>' +
      '<p class="pv-help">Drag a box by its tab, resize it at its corner. Click the tab to edit its text for this step.</p>' +
      '<label class="pv-label" for="pv-text">Editing: <strong data-el="which">nothing picked</strong></label>' +
      '<textarea id="pv-text" data-el="text" rows="8" spellcheck="false" disabled></textarea>' +
      '<p class="pv-help"><code># line</code> label with the blue bar · <code>## line</code> orange title · ' +
        '<code>**bold**</code> · <code>*italic*</code> · <code>[pine] label</code> legend key, with ' +
        Object.keys(PALETTE).join(', ') + '</p>' +
      '<label class="pv-check" data-el="scope" hidden><input type="checkbox" data-el="every"> Show on every step</label>' +
      '<div class="pv-row"><button type="button" data-act="resetBox">Reset this box</button>' +
        '<button type="button" data-act="dropNote" data-el="drop" title="A deleted text box or legend comes back under Show" hidden>Delete</button></div>' +
    '</section>' +
    '<section><h3>Show</h3><div data-el="buttons"></div></section>' +
    '<section data-el="objs" hidden><h3>Structures</h3>' +
      '<p class="pv-help">Everything in the scene. Select one whole, align it to another or to the current selection, ' +
        'rename it, delete it. Fetch or upload more above the viewer.</p>' +
      '<div data-el="list"></div>' +
      '<div class="pv-row"><button type="button" data-act="extract" ' +
        'title="Copy the selected residues into a structure of their own">New structure from selection</button></div>' +
    '</section>' +
    '<section><h3>Selection</h3>' +
      '<p class="pv-help" data-el="sel">The 3D model is not loaded yet.</p>' +
      '<div class="pv-row">' + Object.keys(PALETTE).map(function (n) {
        return '<button type="button" class="pv-swatch" data-act="colour" data-v="' + n + '" title="' + n + '" aria-label="' + n +
               '" style="background: var(' + PALETTE[n] + ')"></button>';
      }).join('') + '<button type="button" data-act="more" aria-pressed="false">More colors</button></div>' +
      '<div data-el="more" hidden>' + Object.keys(MORE).map(function (g) {
        return '<p class="pv-help">' + g + '</p><div class="pv-row">' + MORE[g].split(' ').map(function (h) {
          return '<button type="button" class="pv-swatch" data-act="colour" data-v="#' + h + '" title="#' + h + '" aria-label="#' + h +
                 '" style="background: #' + h + '"></button>';
        }).join('') + '</div>';
      }).join('') + '</div>' +
      '<div class="pv-row"><button type="button" data-act="scheme" data-v="ss">Secondary structure</button>' +
        '<button type="button" data-act="scheme" data-v="plddt" title="Reads the B-factor column as pLDDT">pLDDT</button>' +
        '<button type="button" data-act="scheme" data-v="hydro">Hydropathy</button>' +
        '<button type="button" data-act="scheme" data-v="chain" title="Every chain its own colour">By chain</button>' +
        '<button type="button" data-act="scheme" data-v="sequence" title="Chains with the same sequence share a colour">By sequence</button>' +
        '<button type="button" data-act="uncolour">Original colours</button></div>' +
      '<div class="pv-row">' + Object.keys(MODES).map(function (m) {
        return '<button type="button" data-act="mode" data-v="' + m + '">' + MODES[m][0] + '</button>';
      }).join('') + '<button type="button" data-act="mode" data-v="">Original display</button></div>' +
      '<div class="pv-row"><button type="button" data-act="note">Add note to selection</button>' +
        '<button type="button" data-act="unselect">Clear selection</button></div>' +
    '</section>' +
    '<section><h3>Sequence</h3><div class="pv-row"><button type="button" data-act="seq" aria-pressed="false">Sequence view</button></div></section>' +
    '<section><h3>Scene</h3>' +
      '<p class="pv-help">One JSON file with the structure (an uploaded file travels inside it), colours, text boxes, notes and camera. ' +
        'Import it on any device, in this viewer.</p>' +
      '<div class="pv-row"><button type="button" data-act="export">Export scene</button>' +
      '<button type="button" data-act="import">Import scene</button><button type="button" data-act="wipe">Reset all</button></div></section>');
  document.body.appendChild(panel);
  function $(name) { return panel.querySelector('[data-el="' + name + '"]'); }

  var ticks = {};                                         // name -> its checkbox, so a Delete can untick it
  var toggles = Array.prototype.map.call(host.querySelectorAll('.gvv-tool, .gvv-chips'), function (b) {
    var name = b.getAttribute('data-tool') || 'chips';
    return [name, name === 'chips' ? 'Option chips (bottom left)' : b.getAttribute('aria-label')];
  }).concat(['card', 'legend'].filter(function (id) { return boxes[id]; }).map(function (id) { return [id, boxes[id].label]; }));
  toggles.forEach(function (t) {
    var name = t[0], l = el('label', 'pv-check', '<input type="checkbox"> ');
    ticks[name] = l.firstChild;
    l.append(t[1]);
    l.firstChild.checked = cfg.hide.indexOf(name) < 0;
    l.firstChild.addEventListener('change', function () {
      cfg.hide = cfg.hide.filter(function (n) { return n !== name; }).concat(this.checked ? [] : [name]);
      applyHide(); save();
    });
    $('buttons').appendChild(l);
  });

  function pick(it) {
    picked = it;
    $('which').textContent = it ? it.label : 'nothing picked';
    $('text').disabled = !it;
    $('text').value = it ? it.text() : '';
    $('scope').hidden = !(it && it.note);
    $('drop').hidden = !it;
    if (it && it.note) $('every').checked = it.note.step == null;
  }
  $('text').addEventListener('input', function () { if (picked) { picked.setText(this.value); placeNotes(); save(); } });
  $('every').addEventListener('change', function () { picked.note.step = this.checked ? null : stepKey(); placeNotes(); save(); });

  var ACT = {
    resetBox: function () {
      if (!picked) return;
      if (picked.note) { delete picked.note.w; delete picked.note.h; picked.note.off = null; placeNotes(); }
      else { delete cfg.box[picked.id]; delete cfg.text[picked.key]; picked.draw(); place(picked); pick(picked); }
    },
    dropNote: function () {
      if (!picked.note) {                                 // the text box and the legend are repainted by the viewer: hide, not remove
        var t = ticks[picked.id];
        t.checked = false; t.dispatchEvent(new Event('change'));
        return pick(null);
      }
      cfg.notes.splice(cfg.notes.indexOf(picked.note), 1);
      items.splice(items.indexOf(picked), 1);
      picked.el.remove(); if (picked.frame) picked.frame.remove();
      pick(null); placeNotes();
    },
    colour: function (v) { paint({ color: PALETTE[v] ? hex(PALETTE[v]) : parseInt(v.slice(1), 16) }); },
    more: function () {
      var m = $('more'), b = panel.querySelector('[data-act="more"]');
      b.setAttribute('aria-pressed', m.hidden ? 'true' : 'false');
      m.hidden = !m.hidden;
    },
    scheme: function (v) { paint({ scheme: v }); },
    uncolour: function () { paint({ clear: true }); },
    mode: function (v) { setMode(v); },
    note: function () { addNote(); },
    objSelect: function (id) {
      var r = shown().filter(function (x) { return x.name === id; })[0];
      if (r) plugin.managers.structure.selection.fromLoci('set', S.StructureElement.Loci.all(r.structure));
    },
    objDelete: function (id) { objs.remove(id); },
    extract: function () {
      if (!plugin) return;
      objs.extract().then(function (n) {
        if (!n) alert('Select residues first: click them in the viewer, or use Select on a structure.');
      }).catch(function (e) { alert('Could not make a structure: ' + e.message); });
    },
    unselect: function () { if (plugin) plugin.managers.structure.selection.clear(); },
    seq: function () { cfg.seq = !cfg.seq; drawSeq(); },
    export: function () {
      var k = host.getAttribute('data-pv-key'), own = host.__pvScene && host.__pvScene.get();
      var a = el('a');
      a.href = URL.createObjectURL(new Blob([JSON.stringify(
        Object.assign({ format: SCENE, version: 1, viewer: k }, own, { config: cfg }), null, 1)], { type: 'application/json' }));
      a.download = k + '-scene.json';
      a.click();
    },
    import: function () {
      var f = el('input');
      f.type = 'file'; f.accept = 'application/json,.json';
      f.onchange = function () {
        f.files[0].text().then(importScene).catch(function (e) { alert('Could not import: ' + e.message); });
      };
      f.click();
    },
    wipe: function () {
      if (!confirm('Throw away every edit of this viewer?')) return;
      save = function () {};                              // nothing may write the edits back before the page is gone
      localStorage.removeItem(STORE);
      location.reload();
    }
  };
  /** a scene file, or a bare config from an earlier export. The structure is loaded first, then the edits go over it */
  async function importScene(text) {
    var s = JSON.parse(text), k = host.getAttribute('data-pv-key');
    if (!s || typeof s !== 'object' || Array.isArray(s)) throw new Error('not a scene file');
    if (s.format !== SCENE) s = { config: s };            // before scene files, an export was the config alone
    if (s.viewer && s.viewer !== k) {
      throw new Error('this scene was made in the "' + s.viewer + '" viewer. Open ?v=' + s.viewer + ' and import it there.');
    }
    if (host.__pvScene && (s.structure || s.structures || s.camera)) await host.__pvScene.set(s);
    setConfig(s.config);
    save();
  }
  /** replace every edit in place (cfg keeps its identity: closures hold it) and redraw what depends on it */
  function setConfig(next) {
    var fresh = { box: {}, text: {}, hide: [], paint: [], mode: {}, notes: [], names: {}, seq: false };
    next = next && typeof next === 'object' ? next : {};
    Object.keys(fresh).forEach(function (k) {             // a key of the wrong kind is dropped, not trusted
      var v = next[k], ok = v != null && typeof v === typeof fresh[k] && Array.isArray(v) === Array.isArray(fresh[k]);
      cfg[k] = ok ? v : fresh[k];
    });
    cfg.notes = cfg.notes.filter(function (n) { return n && typeof n.text === 'string'; });
    pick(null);
    for (var i = items.length; i--;) {                    // the old notes go, the new ones are made
      if (!items[i].note) continue;
      items[i].el.remove();
      if (items[i].frame) items[i].frame.remove();
      items.splice(i, 1);
    }
    cfg.notes.forEach(noteItem);
    items.forEach(function (b) { if (!b.note) { b.draw(); place(b); } });
    Object.keys(ticks).forEach(function (n) { ticks[n].checked = cfg.hide.indexOf(n) < 0; });
    applyHide();
    version++; single = null;                             // paint and sequence view are rebuilt from cfg
    if (plugin) sync();
  }
  panel.addEventListener('click', function (e) {
    var b = e.target.closest('[data-act]');
    if (!b) return;
    ACT[b.getAttribute('data-act')](b.getAttribute('data-v'));
    save();
  });

  /* ------------------------------------------------------------ the 3D part -- */
  /* what the viewer has built: [{ name, ref, cell, sref, skey, structure, on }] */
  var idx = null, skeys = {}, order = null;
  /** a name for a structure that survives a reload: its assembly id, or the file it was read from */
  function structKey(c) {
    var ref = c.transform.ref, cells = plugin.state.data.cells, p = c.transform.params;
    if (pv.skey && pv.skey[ref]) return pv.skey[ref];     // the open viewer names its structures itself: files can share a name
    if (skeys[ref]) return skeys[ref];
    if (p && p.type && p.type.params && p.type.params.id) return (skeys[ref] = p.type.params.id);
    for (; c; c = c.transform.parent === c.transform.ref ? null : cells.get(c.transform.parent)) {
      var u = c.transform.params && c.transform.params.url;
      if (u) return (skeys[ref] = String(u.url || u).replace(/^.*\//, ''));
    }
    return '';
  }
  function index() {
    if (idx) return idx;
    var cells = plugin.state.data.cells;
    idx = []; order = null;
    Object.keys(pv.built).forEach(function (name) {       // a value may nest one level: built[tag] = { closed, open, ghost }
      var v = pv.built[name], refs = typeof v === 'string' ? [[name, v]] : Object.keys(v).map(function (k) { return [name + '_' + k, v[k]]; });
      refs.forEach(function (nr) {
        var c = cells.get(nr[1]), sc = c && cells.get(c.transform.parent);
        if (!c || !c.obj || !sc || !sc.obj || !sc.obj.data.units) return;       // gone, or not a structure (the density)
        idx.push({ name: nr[0], ref: nr[1], cell: c, sref: sc.transform.ref, skey: structKey(sc), structure: sc.obj.data,
                   on: !c.state.isHidden });
      });
    });
    return idx;
  }
  function shown() { return index().filter(function (r) { return r.on; }); }

  /* residues of a unit, in order: copies of a chain share one table */
  var ONE = { ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C', GLN: 'Q', GLU: 'E', GLY: 'G', HIS: 'H', ILE: 'I', LEU: 'L',
              LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P', SER: 'S', THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V' };
  var tables = new WeakMap();
  function table(u) {
    var t = tables.get(u.elements);
    if (t) return t;
    var atoms = u.model.atomicHierarchy.atoms, els = u.elements, last = -1;
    t = { seq: '', start: [], res: [], comp: [], ord: new Map() };   // start: first atom (index in the unit) of each residue
    for (var i = 0; i < els.length; i++) {
      var r = u.residueIndex[els[i]];
      if (r === last) continue;
      last = r;
      t.ord.set(r, t.res.length);
      t.start.push(i); t.res.push(r); t.comp.push(atoms.label_comp_id.value(els[i]));
      t.seq += ONE[t.comp[t.comp.length - 1]] || 'X';
    }
    t.start.push(els.length);
    tables.set(els, t);
    return t;
  }

  /* ---- colours: a list of paint jobs, each a target and a colour, a scheme, or "clear" ----
     target { seq, from, to }: those residues of every copy of that protein (from the sequence view)
     target { at: [{ s, u, r }] }: residues r of unit u of structure s (picked in 3D)
     Each representation gets the result as one Mol* overpaint node, rebuilt from the whole list. */
  var KD = { ALA: 1.8, ARG: -4.5, ASN: -3.5, ASP: -3.5, CYS: 2.5, GLN: -3.5, GLU: -3.5, GLY: -0.4, HIS: -3.2, ILE: 4.5, LEU: 3.8,
             LYS: -3.9, MET: 1.9, PHE: 2.8, PRO: -1.6, SER: -0.8, THR: -0.7, TRP: -0.9, TYR: -1.3, VAL: 4.2 };
  var ssProp;
  var SCHEME = {
    ss: function (r, u, t, i) {
      var all = ssProp && ssProp.get(r.structure).value, e = all && all.get(u.invariantId);
      var f = e ? e.type[e.getIndex(t.res[i])] : 0;         // SecondaryStructureType flags: 2 helix, 4 strand
      return hex(f & 2 ? '--gv-helix' : f & 4 ? '--gv-strand' : '--gv-loop');
    },
    hydro: function (r, u, t, i) {                        // Kyte–Doolittle on the viewers' own scale
      var k = ((KD[t.comp[i]] || 0) + 4.5) / 9;
      return k < .5 ? mix(hex('--gv-philic'), hex('--gv-neutral'), k * 2) : mix(hex('--gv-neutral'), hex('--gv-phobic'), k * 2 - 1);
    },
    chain: function (r, u) { return chainColour('chain', u); },
    sequence: function (r, u) { return chainColour('sequence', u); },
    plddt: function (r, u, t, i) {                        // the B-factor column, in AlphaFold's four bands
      var b = u.model.atomicConformation.B_iso_or_equiv.value(u.elements[t.start[i]]);
      return hex(PALETTE[b > 90 ? 'pine' : b > 70 ? 'frost' : b > 50 ? 'peach' : 'cayenne']);
    }
  };
  function chainKey(by, u) {
    return by === 'sequence' ? table(u).seq : u.model.atomicHierarchy.chains.label_asym_id.value(u.chainIndex[u.elements[0]]);
  }
  /** the colour of a unit's chain. Chains are numbered in the order they were built in, hidden steps included, so
      showing another step does not recolour; by 'sequence' every distinct sequence is one number, by 'chain' every label */
  function chainColour(by, u) {
    var reps = index();
    if (!order) {
      order = { chain: new Map(), sequence: new Map() };
      reps.forEach(function (r) {
        r.structure.units.forEach(function (v) {
          if (v.kind !== 0) return;
          Object.keys(order).forEach(function (b) {
            var k = chainKey(b, v);
            if (!order[b].has(k)) order[b].set(k, order[b].size);
          });
        });
      });
    }
    return CHAIN[order[by].get(chainKey(by, u)) % CHAIN.length];
  }
  async function needSS(structure) {
    if (!cfg.paint.some(function (o) { return o.scheme === 'ss'; })) return;
    ssProp = ssProp || plugin.customStructureProperties.get('molstar_computed_secondary_structure');
    if (!ssProp.get(structure).value) {
      await ssProp.attach({ runtime: { shouldUpdate: false, isSynchronous: true, update: function () {} },
                            assetManager: plugin.managers.asset }, structure);
    }
  }

  var version = 0, single = null, atMaps = new WeakMap();
  function atMap(op) {                                    // 's|u' -> residues, for a job picked in 3D
    var m = atMaps.get(op);
    if (!m) { atMaps.set(op, m = new Map()); op.at.forEach(function (a) { m.set(a.s + '|' + a.u, a.r); }); }
    return m;
  }
  /** the colour of each residue of a unit after every job (null: as the viewer drew it); memo: per table, for a structure */
  function resolve(r, u, memo) {
    var t = table(u), id = r.skey + '|' + u.id;
    if (!single) { single = new Set(); cfg.paint.forEach(function (o) { if (o.at) atMap(o).forEach(function (v, k) { single.add(k); }); }); }
    var shared = !single.has(id);
    if (shared && memo.has(t)) return memo.get(t);
    var out = [];
    for (var i = 0; i < t.res.length; i++) out.push(null);
    cfg.paint.forEach(function (o) {
      var set = function (i) { if (i < out.length) out[i] = o.clear ? null : o.scheme ? SCHEME[o.scheme](r, u, t, i) : o.color; };
      if (o.at) (atMap(o).get(id) || []).forEach(set);
      else if (o.seq === t.seq) for (var i = o.from; i <= o.to; i++) set(i);
    });
    if (shared) memo.set(t, out);
    return out;
  }
  function atomsOf(t, from, to) {                         // unit indices of residues from..to: a sorted Int32Array is a Mol* OrderedSet
    var a = new Int32Array(t.start[to + 1] - t.start[from]);
    for (var i = 0; i < a.length; i++) a[i] = t.start[from] + i;
    return a;
  }

  var nodes = {}, painted = {};
  async function applyPaint(r) {
    if (painted[r.ref] === version) return;
    painted[r.ref] = version;
    if (!cfg.paint.length && !nodes[r.ref]) return;
    await needSS(r.structure);
    var memo = new Map(), groups = new Map(), by = new Map(), L = S.StructureElement.Loci, layers = [];
    r.structure.units.forEach(function (u) {
      if (u.kind !== 0) return;
      var t = table(u), cols = resolve(r, u, memo), g = groups.get(cols);
      if (!g) {                                           // colour -> atoms, once per colouring
        var lists = new Map();
        cols.forEach(function (c, i) {
          if (c == null) return;
          var l = lists.get(c) || lists.set(c, []).get(c);
          for (var k = t.start[i]; k < t.start[i + 1]; k++) l.push(k);
        });
        groups.set(cols, g = new Map());
        lists.forEach(function (l, c) { g.set(c, Int32Array.from(l)); });
      }
      g.forEach(function (ix, c) { (by.get(c) || by.set(c, []).get(c)).push({ unit: u, indices: ix }); });
    });
    by.forEach(function (els, c) {
      layers.push({ bundle: S.StructureElement.Bundle.fromLoci(L(r.structure, els)), color: c, clear: false });
    });
    var b = plugin.build();
    [r.ref, skins[r.ref]].forEach(function (ref) {        // the cartoon and, in "Cartoon + surface", the surface over it
      if (!ref || !plugin.state.data.cells.has(ref)) return;
      if (nodes[ref]) b.to(nodes[ref]).update({ layers: layers });
      else nodes[ref] = b.to(ref).apply(pv.lib.plugin.StateTransforms.Representation.OverpaintStructureRepresentation3DFromBundle,
                                        { layers: layers }, { tags: 'pv-paint' }).ref;
    });
    await b.commit();
  }

  /* ---- display mode: per component, by the name its viewer built it under ---- */
  var was = {}, now = {}, skins = {};                     // skins: component's ref -> the see-through surface drawn over its cartoon
  var SKIN = 'cartoon-surface', SKIN_ALPHA = 0.4;
  async function applyMode(r) {
    var want = cfg.mode[r.name] || null, cells = plugin.state.data.cells;
    if ((now[r.ref] || null) === want) return;
    var old = r.cell.transform.params.type;
    if (!was[r.ref]) was[r.ref] = old;
    var kind = want === SKIN ? 'cartoon' : want;
    var type = !want ? was[r.ref]
      : kind === was[r.ref].name ? was[r.ref]             // the viewer's own cartoon, as it drew it
      : { name: kind, params: Object.assign({}, plugin.representation.structure.registry.get(kind).defaultValues, LOOK, MODES[kind][1]) };
    now[r.ref] = want;
    var at = r.cell.obj.data.repr.state.transform.slice();   // where its viewer has put it: a new type starts at the origin again
    await plugin.build().to(r.ref).update(function (p) {
      p.type = { name: type.name, params: Object.assign({}, type.params, { alpha: old.params.alpha }) };   // a fade keeps running
      return p;
    }).commit();
    var repr = r.cell.obj.data.repr;
    repr.setState({ transform: at });
    plugin.canvas3d.update(repr);
    if (skins[r.ref] && (want !== SKIN || !cells.has(skins[r.ref]))) {
      if (cells.has(skins[r.ref])) await plugin.build().delete(skins[r.ref]).commit();
      skins[r.ref] = null;
    }
    if (want === SKIN && !skins[r.ref]) {                 // the surface takes the cartoon's colour theme, so paint applies to both
      var theme = cells.get(r.ref).transform.params.colorTheme;
      var s = await plugin.builders.structure.representation.addRepresentation(r.sref, { type: 'gaussian-surface',
        typeParams: Object.assign({}, MODES['gaussian-surface'][1], LOOK, { alpha: SKIN_ALPHA }) });
      await plugin.build().to(s.ref).update(function (p) { p.colorTheme = theme; return p; }).commit();
      skins[r.ref] = s.ref;
      delete painted[r.ref];                              // paint the new surface too
    }
  }
  function skinFollow() {                                 // ... and so does its pose: the viewers move parts with the representation's transform
    var cells = plugin.state.data.cells;
    Object.keys(skins).forEach(function (ref) {
      var a = skins[ref] && cells.get(ref), b = skins[ref] && cells.get(skins[ref]);
      if (!a || !b || !a.obj || !b.obj) return;
      var m = a.obj.data.repr.state.transform, n = b.obj.data.repr.state.transform;
      if (m.some(function (v, i) { return v !== n[i]; })) {
        b.obj.data.repr.setState({ transform: m.slice() });
        plugin.canvas3d.requestDraw();
      }
    });
  }
  function skinVisibility(r) {                            // the viewers show and hide the cartoon; the surface follows it
    var s = skins[r.ref], c = s && plugin.state.data.cells.get(s);
    if (c && !!c.state.isHidden !== !r.on) plugin.state.data.updateCellState(s, { isHidden: !r.on });
  }
  function setMode(type) {
    if (!plugin) return;
    var hit = sel.n ? sel.at.map(function (a) { return a.s; }) : null;
    shown().forEach(function (r) {
      if (hit && hit.indexOf(r.skey) < 0) return;
      if (type) cfg.mode[r.name] = type; else delete cfg.mode[r.name];
    });
    sync();
  }

  var queue = Promise.resolve(), timer;
  function touch() { idx = null; clearTimeout(timer); timer = setTimeout(sync, 120); }
  function sync() {
    queue = queue.then(async function () {
      idx = null;
      for (var reps = index(), i = 0; i < reps.length; i++) { await applyMode(reps[i]); await applyPaint(reps[i]); skinVisibility(reps[i]); }
      drawSeq(); placeNotes();
    }).catch(function (e) { console.error('[dev-mode]', e); });
  }

  /* ---- selection: Mol*'s own. Clicks in 3D (selection mode) and drags over the sequence both end up in it ---- */
  var sel = { n: 0, at: [], bySeq: {} };
  var seqSel = null;                                      // [{ seq, from, to }] while the selection came from the sequence view only
  var selfSel = false, selTimer;
  function readSelection() {
    var L = S.StructureElement.Loci, cells = plugin.state.data.cells;
    sel = { n: 0, at: [], bySeq: {} };
    plugin.managers.structure.selection.entries.forEach(function (e, ref) {
      var c = cells.get(ref), per = new Map();
      if (!c || L.isEmpty(e.selection)) return;
      L.forEachLocation(e.selection, function (loc) {
        if (loc.unit.kind !== 0) return;
        var set = per.get(loc.unit) || per.set(loc.unit, new Set()).get(loc.unit);
        set.add(table(loc.unit).ord.get(loc.unit.residueIndex[loc.element]));
      });
      per.forEach(function (set, u) {
        var seq = table(u).seq, all = sel.bySeq[seq] || (sel.bySeq[seq] = new Set());
        set.forEach(function (i) { all.add(i); });
        sel.n += set.size;
        sel.at.push({ s: structKey(c), u: u.id, r: Array.from(set) });
      });
    });
    $('sel').textContent = (sel.n ? sel.n + ' residue' + (sel.n > 1 ? 's' : '') + ' selected.' : 'Nothing selected: colours and display apply to everything shown.') +
      ' Click residues in the viewer (Shift for a range) or drag over the sequence.';
    markSeq();
  }
  function selectSeq(seq, from, to, add) {
    var m = plugin.managers.structure.selection, done = {}, x = { seq: seq, from: from, to: to };
    seqSel = !add ? [x] : seqSel ? seqSel.concat([x]) : sel.n ? null : [x];
    selfSel = true;
    if (!add) m.clear();
    shown().forEach(function (r) {
      if (done[r.sref]) return;
      done[r.sref] = true;
      var memo = new Map(), els = [];
      r.structure.units.forEach(function (u) {
        var t = u.kind === 0 && table(u);
        if (!t || t.seq !== seq) return;
        els.push({ unit: u, indices: memo.get(t) || memo.set(t, atomsOf(t, from, to)).get(t) });
      });
      if (els.length) m.fromLoci('add', S.StructureElement.Loci(r.structure, els));
    });
    selfSel = false;
  }
  function paint(what) {
    if (!plugin) return;
    var targets = seqSel || (sel.n ? [{ at: sel.at }] : rows().map(function (w) { return { seq: w.seq, from: 0, to: w.seq.length - 1 }; }));
    targets.forEach(function (t) {                        // a job that is painted over completely is dropped
      if (t.seq != null) cfg.paint = cfg.paint.filter(function (o) { return !(o.seq === t.seq && o.from >= t.from && o.to <= t.to); });
      cfg.paint.push(Object.assign({}, what, t));
    });
    if (cfg.paint.every(function (o) { return o.clear; })) cfg.paint = [];
    version++; single = null;
    sync();
  }

  /* ---- sequence view: one row per protein (= distinct sequence) on screen ---- */
  var seqEl = el('div', 'pv-seq'), seqRows = [], seqSig = '';
  seqEl.hidden = true;
  host.after(seqEl);
  function rows() {
    var seen = {}, out = [];
    shown().forEach(function (r) {
      var done = new Set();
      r.structure.units.forEach(function (u) {
        if (u.kind !== 0 || done.has(u.elements)) return;
        done.add(u.elements);
        var t = table(u), w = seen[t.seq];
        if (!w) out.push(w = seen[t.seq] = { seq: t.seq, chains: [], rep: r, unit: u });
        var ch = u.model.atomicHierarchy.chains.label_asym_id.value(u.chainIndex[u.elements[0]]);
        if (w.chains.indexOf(ch) < 0) w.chains.push(ch);
      });
    });
    return out;
  }
  function drawSeq() {
    panel.querySelector('[data-act="seq"]').setAttribute('aria-pressed', cfg.seq ? 'true' : 'false');
    seqEl.hidden = !(dev && cfg.seq && plugin);
    if (seqEl.hidden) return;
    var rs = rows(), sig = version + '|' + rs.map(function (w) { return w.seq + w.chains; }).join('|');
    if (sig !== seqSig) {
      seqSig = sig; seqRows = rs;
      seqEl.innerHTML = rs.length ? '' : '<p class="pv-help">No protein on screen.</p>';
      rs.forEach(function (w, n) {
        var cols = resolve(w.rep, w.unit, new Map()), res = w.unit.model.atomicHierarchy.residues, t = table(w.unit);
        var row = el('div', 'pv-seq-row',
          '<input class="pv-seq-name" aria-label="Protein name"><span class="pv-help">' + w.seq.length + ' residues · chain ' +
          esc(w.chains.join(', ')) + '</span><div class="pv-seq-res" data-row="' + n + '">' + w.seq.split('').map(function (ch, i) {
            var c = cols[i], dark = c != null && (c >> 16 & 255) * .3 + (c >> 8 & 255) * .59 + (c & 255) * .11 < 140;
            return '<span data-i="' + i + '"' + ((i + 1) % 10 ? '' : ' data-n="' + (i + 1) + '"') +
              ' title="' + t.comp[i] + ' ' + res.auth_seq_id.value(t.res[i]) + '"' +
              (c == null ? '' : ' style="background:' + css(c) + (dark ? ';color:var(--clr-on-dark)' : '') + '"') + '>' + ch + '</span>';
          }).join('') + '</div>');
        row.firstChild.value = cfg.names[w.seq] || 'Chain ' + w.chains.join(', ');
        row.firstChild.addEventListener('input', function () { cfg.names[w.seq] = this.value; save(); });
        seqEl.appendChild(row);
      });
    }
    markSeq();
  }
  function markSeq() {
    Array.prototype.forEach.call(seqEl.querySelectorAll('.pv-seq-res'), function (row) {
      var set = sel.bySeq[seqRows[row.getAttribute('data-row')].seq];
      Array.prototype.forEach.call(row.children, function (c, i) {
        c.classList.remove('is-pre');
        c.classList.toggle('is-sel', !!set && set.has(i));
      });
    });
  }
  seqEl.addEventListener('pointerdown', function (e) {    // drag over residues: a range; Shift adds it to the selection
    var s = e.target.closest('[data-i]');
    if (!s) return;
    e.preventDefault();
    var row = s.parentNode, a = +s.getAttribute('data-i'), b = a;
    function preview() {
      Array.prototype.forEach.call(row.children, function (c, i) { c.classList.toggle('is-pre', i >= Math.min(a, b) && i <= Math.max(a, b)); });
    }
    function over(ev) {
      var o = ev.target.closest('[data-i]');
      if (o && o.parentNode === row) { b = +o.getAttribute('data-i'); preview(); }
    }
    function up(ev) {
      row.removeEventListener('pointerover', over);
      window.removeEventListener('pointerup', up);
      selectSeq(seqRows[row.getAttribute('data-row')].seq, Math.min(a, b), Math.max(a, b), ev.shiftKey);
    }
    row.addEventListener('pointerover', over);
    window.addEventListener('pointerup', up);
    preview();
  });

  /* ---- structures (the open viewer): the list, and aligning one to another ---- */
  function drawObjs() {
    var l = objs.list(), box = $('list');
    $('objs').hidden = false;
    box.innerHTML = l.length ? '' : '<p class="pv-help">No structure yet.</p>';
    l.forEach(function (o) {
      var id = esc(o.id), row = el('div', 'pv-obj',
        '<input class="pv-obj-name" aria-label="Structure name">' +
        '<div class="pv-row"><button type="button" data-act="objSelect" data-v="' + id + '">Select</button>' +
        '<select data-align="' + id + '" aria-label="Align to"><option value="">Align to…</option>' +
          '<option value="sel">the selection</option>' + l.filter(function (x) { return x.id !== o.id; }).map(function (x) {
            return '<option value="o:' + esc(x.id) + '">' + esc(x.name) + '</option>';
          }).join('') + '</select>' +
        '<button type="button" data-act="objDelete" data-v="' + id + '">Delete</button></div>');
      row.firstChild.value = o.name;
      row.firstChild.addEventListener('input', function () {
        objs.rename(o.id, this.value);
        Array.prototype.forEach.call(box.querySelectorAll('option'), function (op) {   // the other rows name it too
          if (op.value === 'o:' + o.id) op.textContent = row.firstChild.value;
        });
      });
      box.appendChild(row);
    });
  }
  $('list').addEventListener('change', function (e) {
    var s = e.target, to = s.value;
    if (!to || !s.hasAttribute('data-align')) return;
    s.value = '';
    align(s.getAttribute('data-align'), to).catch(function (err) { alert('Could not align: ' + err.message); });
  });

  /** the residues of a shown structure, with its pose on: their one-letter codes and the position of each one's CA
      (C4' in nucleic acids; its first atom otherwise). keep: unit id -> the residues (as numbered in table()) to take */
  function trace(r, keep) {
    var P = r.cell.obj.data.repr.state.transform, out = { seq: '', xyz: [] };
    r.structure.units.forEach(function (u) {
      var only = keep && keep.get(u.id);
      if (u.kind !== 0 || (keep && !only)) return;
      var t = table(u), names = u.model.atomicHierarchy.atoms.label_atom_id;
      for (var i = 0; i < t.res.length; i++) {
        if (only && !only.has(i)) continue;
        var k = t.start[i];
        for (var j = t.start[i]; j < t.start[i + 1]; j++) {
          var n = names.value(u.elements[j]);
          if (n === 'CA' || n === "C4'") { k = j; break; }
        }
        out.seq += t.seq[i];
        out.xyz.push(V3.transformMat4(V3(), u.conformation.position(u.elements[k], V3()), P));
      }
    });
    return out;
  }
  /** Smith-Waterman on two one-letter sequences: [[i, j], …] of the residues it pairs up. ponytail: a matrix of
      n * m bytes, so a few thousand residues each is the end of it; take a chain, or the selection, first */
  function match(a, b) {
    if (a === b) return Array.from(a, function (c, i) { return [i, i]; });          // the same construct, another pose
    var n = a.length, m = b.length, W = m + 1;
    if (n * m > 4e7) throw new Error('too many residues to pair up: select a chain of each first');
    var tb = new Uint8Array((n + 1) * W), up = new Int32Array(W), now = new Int32Array(W), best = 0, bi = 0, bj = 0;
    for (var i = 1; i <= n; i++) {
      for (var j = 1; j <= m; j++) {
        var d = up[j - 1] + (a[i - 1] === b[j - 1] ? 2 : -1), v = up[j] - 2, h = now[j - 1] - 2, s = Math.max(0, d, v, h);
        now[j] = s;
        tb[i * W + j] = s === 0 ? 0 : s === d ? 1 : s === v ? 2 : 3;                // 1 pair, 2 skip in a, 3 skip in b
        if (s > best) { best = s; bi = i; bj = j; }
      }
      var swap = up; up = now; now = swap;
    }
    var out = [];
    for (var x = bi, y = bj, k; (k = tb[x * W + y]);) {
      if (k === 1) { out.push([x - 1, y - 1]); x--; y--; } else if (k === 2) x--; else y--;
    }
    return out.reverse();
  }
  /** the rigid motion (a Mat4) that takes points a onto points b with the least squared distance: Horn's method,
      the quaternion is the top eigenvector of a 4 x 4 matrix of the correlations, found by Jacobi sweeps */
  function superpose(a, b) {
    var n = a.length, ca = V3(), cb = V3(), S = [0, 0, 0, 0, 0, 0, 0, 0, 0], i, j, k;
    for (i = 0; i < n; i++) { V3.add(ca, ca, a[i]); V3.add(cb, cb, b[i]); }
    V3.scale(ca, ca, 1 / n); V3.scale(cb, cb, 1 / n);
    for (i = 0; i < n; i++) {
      for (j = 0; j < 3; j++) for (k = 0; k < 3; k++) S[3 * j + k] += (a[i][j] - ca[j]) * (b[i][k] - cb[k]);
    }
    var N = [[S[0] + S[4] + S[8], S[5] - S[7], S[6] - S[2], S[1] - S[3]],
             [S[5] - S[7], S[0] - S[4] - S[8], S[1] + S[3], S[6] + S[2]],
             [S[6] - S[2], S[1] + S[3], -S[0] + S[4] - S[8], S[5] + S[7]],
             [S[1] - S[3], S[6] + S[2], S[5] + S[7], -S[0] - S[4] + S[8]]];
    var Q = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
    for (var sweep = 0; sweep < 50; sweep++) {
      var moved = false;
      for (var p = 0; p < 3; p++) for (var q = p + 1; q < 4; q++) {
        if (Math.abs(N[p][q]) < 1e-12) continue;
        moved = true;
        var th = (N[q][q] - N[p][p]) / (2 * N[p][q]), t = (th < 0 ? -1 : 1) / (Math.abs(th) + Math.sqrt(th * th + 1));
        var c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (k = 0; k < 4; k++) {                                                   // N <- J' N J, Q <- Q J
          var x = N[k][p], y = N[k][q];
          N[k][p] = c * x - s * y; N[k][q] = s * x + c * y;
        }
        for (k = 0; k < 4; k++) {
          var x2 = N[p][k], y2 = N[q][k];
          N[p][k] = c * x2 - s * y2; N[q][k] = s * x2 + c * y2;
          var x3 = Q[k][p], y3 = Q[k][q];
          Q[k][p] = c * x3 - s * y3; Q[k][q] = s * x3 + c * y3;
        }
      }
      if (!moved) break;
    }
    var top = 0;
    for (i = 1; i < 4; i++) if (N[i][i] > N[top][top]) top = i;
    var w = Q[0][top], qx = Q[1][top], qy = Q[2][top], qz = Q[3][top];
    var R = [[1 - 2 * (qy * qy + qz * qz), 2 * (qx * qy - w * qz), 2 * (qx * qz + w * qy)],
             [2 * (qx * qy + w * qz), 1 - 2 * (qx * qx + qz * qz), 2 * (qy * qz - w * qx)],
             [2 * (qx * qz - w * qy), 2 * (qy * qz + w * qx), 1 - 2 * (qx * qx + qy * qy)]];
    var m = M4.identity();
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) m[4 * j + i] = R[i][j];
      m[12 + i] = cb[i] - (R[i][0] * ca[0] + R[i][1] * ca[1] + R[i][2] * ca[2]);
    }
    return m;
  }
  /** move structure `id` onto another ("o:<id>") or onto the selected residues ("sel"): its residues are paired with
      theirs by sequence, then the pose is the superposition of the pairs */
  async function align(id, to) {
    var me = shown().filter(function (r) { return r.name === id; })[0];
    if (!me) throw new Error('that structure is not shown');
    var target = { seq: '', xyz: [] };
    shown().forEach(function (r) {
      if (r === me || (to !== 'sel' && 'o:' + r.name !== to)) return;
      var keep = null;
      if (to === 'sel') {                                 // only what is selected of it
        keep = new Map();
        sel.at.forEach(function (a) { if (a.s === r.skey) keep.set(a.u, new Set(a.r)); });
        if (!keep.size) return;
      }
      var t = trace(r, keep);
      target.seq += t.seq; target.xyz = target.xyz.concat(t.xyz);
    });
    if (!target.seq) throw new Error(to === 'sel' ? 'select residues of another structure first.' : 'there is nothing to align to.');
    var mine = trace(me), pairs = match(mine.seq, target.seq);
    if (pairs.length < 3) throw new Error('fewer than 3 residues of the two match.');
    var m = superpose(pairs.map(function (p) { return mine.xyz[p[0]]; }), pairs.map(function (p) { return target.xyz[p[1]]; }));
    objs.pose(id, M4.mul(M4(), m, me.cell.obj.data.repr.state.transform));
  }

  /* ---- notes: a box at a 3D point (its residue + an offset), behind the scene like every box; the line to its
     residue is drawn above the scene, or it would end at the protein's outline ----
     { id, s, u, a: atom, off: [x, y, z] Å, text, step: the step it belongs to (null: every step), w, h px } */
  var noteLayer = el('div', 'pv-notes'), leaders = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  leaders.setAttribute('class', 'pv-leaders');
  if (stage) { stage.prepend(noteLayer); stage.appendChild(leaders); }

  function noteItem(n) {
    var e = el('div', 'pv-note'), it = { id: n.id, label: 'Note', el: e, note: n };
    it.text = function () { return n.text; };
    it.setText = function (t) { n.text = t; e.innerHTML = render(t, 'span'); };
    it.grab = function () { return { x: it.x, y: it.y, w: e.offsetWidth, h: e.offsetHeight }; };
    it.move = function (s, dx, dy) { n.off = Array.from(V3.sub(V3(), unproject(s.x + dx, s.y + dy, it.z), it.at)); placeNotes(); };
    it.size = function (s, dx, dy) { n.w = Math.max(60, s.w + dx); n.h = Math.max(24, s.h + dy); placeNotes(); };
    it.setText(n.text);
    e.hidden = true;
    noteLayer.appendChild(e);
    items.push(it);
    return it;
  }
  cfg.notes.forEach(noteItem);

  function unproject(x, y, z) {
    return V3.transformMat4(V3(), V3.create(x / stage.clientWidth * 2 - 1, 1 - y / stage.clientHeight * 2, z),
                            plugin.canvas3d.camera.inverseProjectionView);
  }
  function anchor(n) {                                    // where its residue is drawn right now, or null if it is not
    var r = shown().filter(function (x) { return x.skey === n.s; })[0], u = r && r.structure.unitMap.get(n.u);
    if (!u) return null;
    var v = u.conformation.position(n.a, V3());
    return V3.transformMat4(v, v, r.cell.obj.data.repr.state.transform);   // the viewers move parts by this matrix
  }
  function placeNotes() {
    if (!plugin) return;
    var w = stage.clientWidth, h = stage.clientHeight, pvm = plugin.canvas3d.camera.projectionView, lines = '';
    var px = function (v) { var q = V3.transformMat4(V3(), v, pvm); return [(q[0] * .5 + .5) * w, (1 - (q[1] * .5 + .5)) * h, q[2]]; };
    items.forEach(function (it) {
      var n = it.note, e = it.el;
      if (!n) return;
      var at = (n.step == null || n.step === stepKey()) ? anchor(n) : null, a = at && px(at);
      e.hidden = !a || Math.abs(a[2]) > 1;
      if (e.hidden) return;
      if (!n.off) n.off = Array.from(V3.sub(V3(), unproject(a[0] + 60, a[1] - 90, a[2]), at));    // a new note: up and to the right
      var p = px(V3.add(V3(), at, n.off));
      it.at = at; it.x = p[0]; it.y = p[1]; it.z = p[2];     // where it is drawn, for a drag to start from
      e.style.cssText = 'left:' + p[0] + 'px;top:' + p[1] + 'px' + (n.w ? ';width:' + n.w + 'px;max-width:none' : '') + (n.h ? ';height:' + n.h + 'px' : '');
      // above the canvas only where nothing is drawn at the box's centre, or what is there lies behind it (one answer per box)
      var hit = plugin.canvas3d.identify([p[0] + e.offsetWidth / 2, p[1] + e.offsetHeight / 2]);
      e.classList.toggle('is-front', !hit || !hit.position || p[2] <= V3.transformMat4(V3(), hit.position, pvm)[2]);
      var x = Math.min(Math.max(a[0], p[0]), p[0] + e.offsetWidth), y = Math.min(Math.max(a[1], p[1]), p[1] + e.offsetHeight);
      lines += '<line x1="' + x + '" y1="' + y + '" x2="' + a[0] + '" y2="' + a[1] + '"/><circle cx="' + a[0] + '" cy="' + a[1] + '" r="3"/>';
    });
    if (leaders.innerHTML !== lines) leaders.innerHTML = lines;
  }
  function addNote() {
    var a = sel.at[0], r = a && shown().filter(function (x) { return x.skey === a.s; })[0];
    if (!r) return alert('Select a residue first: click it in the viewer.');
    var u = r.structure.unitMap.get(a.u), t = table(u), i = a.r[0], h = u.model.atomicHierarchy, atom = u.elements[t.start[i]];
    for (var k = t.start[i]; k < t.start[i + 1]; k++) if (h.atoms.label_atom_id.value(u.elements[k]) === 'CA') atom = u.elements[k];
    var n = { id: Date.now().toString(36), s: a.s, u: a.u, a: atom, off: null, step: stepKey(),
              text: '## ' + t.comp[i] + ' ' + h.residues.auth_seq_id.value(t.res[i]) + '\nWhat this residue does.' };
    cfg.notes.push(n);
    pick(noteItem(n));
    placeNotes();
    $('text').focus();
  }

  /* ------------------------------------------------------------- on and off -- */
  var toggle = document.querySelector('[data-pv-dev]');
  function setDev(on) {
    dev = on;
    localStorage.setItem('pv-dev', on ? '1' : '');
    document.body.classList.toggle('pv-dev', on);
    if (toggle) toggle.setAttribute('aria-pressed', on ? 'true' : 'false');
    if (plugin) {
      plugin.selectionMode = on;
      if (on) readSelection(); else plugin.managers.structure.selection.clear();
      drawSeq();
    }
    if (on && !looping) syncFrames();
  }
  if (toggle) toggle.addEventListener('click', function () { setDev(!dev); });

  function ready() {
    pv = host.__pv; plugin = pv.plugin; S = pv.lib.structure; V3 = pv.lib.math.LinearAlgebra.Vec3;
    M4 = pv.lib.math.LinearAlgebra.Mat4;
    var ev = plugin.state.data.events;
    ev.cell.stateUpdated.subscribe(touch);               // a part was shown or hidden
    ev.object.created.subscribe(touch);                  // ... or built
    plugin.managers.structure.selection.events.changed.subscribe(function () {
      if (!selfSel) seqSel = null;
      clearTimeout(selTimer); selTimer = setTimeout(readSelection, 30);
    });
    plugin.canvas3d.didDraw.subscribe(placeNotes);
    plugin.canvas3d.didDraw.subscribe(skinFollow);
    setDev(dev);
    sync();
  }
  if (objs) { host.addEventListener('pv:objects', drawObjs); drawObjs(); }
  setDev(localStorage.getItem('pv-dev') === '1');
  if (host.__pv) ready(); else host.addEventListener('pv:ready', ready);
})();
