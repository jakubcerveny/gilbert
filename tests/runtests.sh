#!/bin/bash
#
# SPDX-License-Identifier: BSD-2-Clause
# Copyright (c) 2018 Jakub Červený
#


ln -f -s ../gilbert2d.py .
ln -f -s ../gilbert3d.py .
ln -f -s ../gilbert.js .
ln -f -s ../gilbert .

gilbert_cmp2 () {
  local x=$1
  local y=$2

  echo -n "(python) xy2d[$x,$y]: "
  diff \
    <( ./gilbert2d.py $x $y 2> /dev/null ) \
    <( ./gilbert2d.py --op xy2d $x $y | sort -n | cut -f2- -d' ' 2> /dev/null ) > /dev/null
  if [[ $? != 0 ]] ; then echo "FAIL" ; else echo "pass" ; fi

  echo -n "(python) d2xy[$x,$y]: "
  diff \
    <( ./gilbert2d.py $x $y 2> /dev/null ) \
    <( ./gilbert2d.py --op d2xy $x $y 2> /dev/null ) > /dev/null
  if [[ $? != 0 ]] ; then echo "FAIL" ; else echo "pass" ; fi

  ###

  echo -n "(js) xy2d[$x,$y]: "
  diff \
    <( ./gilbert2d.py $x $y 2> /dev/null ) \
    <( node -e 'require("./gilbert.js").main(["gilbert.js","xy2d",'$x','$y']);' | sort -n | cut -f2- -d' ' 2> /dev/null ) > /dev/null
  if [[ $? != 0 ]] ; then echo "FAIL" ; else echo "pass" ; fi

  echo -n "(js) d2xy[$x,$y]: "
  diff \
    <( ./gilbert2d.py $x $y 2> /dev/null ) \
    <( node -e 'require("./gilbert.js").main(["gilbert.js","d2xy",'$x','$y']);' 2> /dev/null ) > /dev/null
  if [[ $? != 0 ]] ; then echo "FAIL" ; else echo "pass" ; fi

  ###

  echo -n "(c) xy2d[$x,$y]: "
  diff \
    <( ./gilbert2d.py $x $y 2> /dev/null ) \
    <( ./gilbert xy2d $x $y | sort -n | cut -f2- -d' ' 2> /dev/null ) > /dev/null
  if [[ $? != 0 ]] ; then echo "FAIL" ; else echo "pass" ; fi

  echo -n "(c) d2xy[$x,$y]: "
  diff \
    <( ./gilbert2d.py $x $y 2> /dev/null ) \
    <( ./gilbert d2xy $x $y 2> /dev/null ) > /dev/null
  if [[ $? != 0 ]] ; then echo "FAIL" ; else echo "pass" ; fi

}

gilbert_cmp3 () {
  local x=$1
  local y=$2
  local z=$3

  echo -n "(python) xyz2d[$x,$y,$z]: "
  diff \
    <( ./gilbert3d.py $x $y $z 2> /dev/null ) \
    <( ./gilbert3d.py --op xyz2d $x $y $z | sort -n | cut -f2- -d' ' 2> /dev/null ) > /dev/null
  if [[ $? != 0 ]] ; then echo "FAIL" ; else echo "pass" ; fi

  echo -n "(python) d2xyz[$x,$y,$z]: "
  diff \
    <( ./gilbert3d.py $x $y $z 2> /dev/null ) \
    <( ./gilbert3d.py --op d2xyz $x $y $z 2> /dev/null ) > /dev/null
  if [[ $? != 0 ]] ; then echo "FAIL" ; else echo "pass" ; fi

  ###

  echo -n "(js) xyz2d[$x,$y,$z]: "
  diff \
    <( ./gilbert3d.py $x $y $z 2> /dev/null ) \
    <( node -e 'require("./gilbert.js").main(["gilbert.js","xyz2d",'$x','$y','$z']);' | sort -n | cut -f2- -d' ' 2> /dev/null ) > /dev/null
  if [[ $? != 0 ]] ; then echo "FAIL" ; else echo "pass" ; fi

  echo -n "(js) d2xyz[$x,$y,$z]: "
  diff \
    <( ./gilbert3d.py $x $y $z 2> /dev/null ) \
    <( node -e 'require("./gilbert.js").main(["gilbert.js","d2xyz",'$x','$y','$z']);' 2> /dev/null ) > /dev/null
  if [[ $? != 0 ]] ; then echo "FAIL" ; else echo "pass" ; fi

  ###

  echo -n "(c) xyz2d[$x,$y,$z]: "
  diff \
    <( ./gilbert3d.py $x $y $z 2> /dev/null ) \
    <( ./gilbert xyz2d $x $y $z | sort -n | cut -f2- -d' ' 2> /dev/null ) > /dev/null
  if [[ $? != 0 ]] ; then echo "FAIL" ; else echo "pass" ; fi

  echo -n "(c) d2xyz[$x,$y,$z]: "
  diff \
    <( ./gilbert3d.py $x $y $z 2> /dev/null ) \
    <( ./gilbert d2xyz $x $y $z 2> /dev/null ) > /dev/null
  if [[ $? != 0 ]] ; then echo "FAIL" ; else echo "pass" ; fi

}

x=100 ; y=63
gilbert_cmp2 $x $y

x=8 ; y=6 ; z=4
gilbert_cmp3 $x $y $z

x=40 ; y=30
gilbert_cmp2 $x $y

x=40 ; y=30 ; z=20
gilbert_cmp3 $x $y $z

x=20 ; y=12 ; z=2
gilbert_cmp3 $x $y $z

x=15 ; y=12
gilbert_cmp2 $x $y

x=7 ; y=6 ; z=4
gilbert_cmp3 $x $y $z

