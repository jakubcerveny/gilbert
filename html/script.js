// SPDX-License-Identifier: BSD-2-Clause
// Copyright (c) 2018 Jakub Červený

var info = {
  "two": null,
  "line": [],
  "two": null,
  "ctx": null,
  "canvas": null
};

(()=>{

  let W = 15, H = 13;
  let s = { "x": 10, "y": 10 };
  let T = { "x": 10, "y": 10 };

  let container = document.getElementById("gilbert_container");
  let two = new Two().appendTo(container);

  info.two = two;
  
  for (let idx=1; idx<(W*H); idx++) {
    let q = gilbert.d2xy( idx-1, W, H );
    let p = gilbert.d2xy( idx, W, H );

    q.x *= s.x; q.y *= s.y;
    p.x *= s.x; p.y *= s.y;

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
  info.svg = container.innerHTML;

})();
