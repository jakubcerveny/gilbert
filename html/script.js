// SPDX-License-Identifier: BSD-2-Clause
// Copyright (c) 2018 Jakub Červený

var info = {
  "W" : -1,
  "H": -1,
  "default": { "w": 28, "h": 18 },
  "T": { "x": 10, "y": 10 },
  "S": { "x": 10, "y": 10 },
  "two": null,
  "line": [],
  "two": null,
  "ctx": null,
  "container" : null,
  "canvas": null
};

function update_wh(w,h) {
  let ui_width  = document.getElementById("ui_width");
  let ui_height = document.getElementById("ui_height");

  ui_width.value  = w;
  ui_height.value = h

  info.W = w;
  info.H = h;

  draw_curve();
}


function update_preset() {
  let ele = document.getElementById("ui_preset");
  let v = ele.value;

  if ( v.match( /^\d+x\d+$/ ) ){
    let tok = v.split("x");
    let _w = parseInt(tok[0]);
    let _h = parseInt(tok[1]);

    if (_w < 1) { _w = 1; }
    if (_h < 1) { _h = 1; }

    update_wh(_w,_h);
    return;
  }

  let ui_w = document.getElementById("ui_width");
  let ui_h = document.getElementById("ui_height");

  let _w = parseInt(ui_w.value);
  let _h = parseInt(ui_h.value);

  if (_w < 1) { _w = 1; }
  if (_h < 1) { _h = 1; }

  update_wh(_w,_h);

}

function update_num() {
  let ui_w = document.getElementById("ui_width");
  let ui_h = document.getElementById("ui_height");

  let _w = parseInt(ui_w.value);
  let _h = parseInt(ui_h.value);

  if (_w < 1) { _w = 1; }
  if (_h < 1) { _h = 1; }

  update_wh(_w,_h);
}


function draw_curve() {
  let two = info.two;
  let W = info.W;
  let H = info.H;

  let S = info.S;
  let T = info.T;

  two.clear();

  for (let idx=1; idx<(W*H); idx++) {
    let q = gilbert.d2xy( idx-1, W, H );
    let p = gilbert.d2xy( idx, W, H );

    q.x *= S.x; q.y *= S.y;
    p.x *= S.x; p.y *= S.y;

    q.x += T.x; q.y += T.y;
    p.x += T.x; p.y += T.y;

    let _line = two.makeLine(q.x, q.y, p.x, p.y);

    let hue = Math.floor(360*idx / (W*H)).toString();
    let lit = '55%';
    let crm = '0.35';
    let clr = 'oklch(' + [ lit,crm,hue ].join(" ") + ')';

    _line.stroke = clr;

    info.line.push( _line );
  }

  two.update();
  info.svg = info.container.innerHTML;
}


function init() {
  let W = info.default.w,
      H = info.default.h;

  info.container = document.getElementById("gilbert_container");
  info.two = new Two().appendTo(info.container);

  update_wh(W,H);
}

init();
