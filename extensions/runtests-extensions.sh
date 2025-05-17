#!/bin/bash
#

N=13

COMPARE_TESTS=1

C_TESTS=1
PYTHON_TESTS=1
JAVASCRIPT_TESTS=1
VERBOSE_LEVEL=1
ERR=0

CHECK_OPT=""
if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
  CHECK_OPT='verbose'
fi

gilbert2dpp_cmp () {
  local bin_a=$1
  local bin_b=$2
  local x=$3
  local y=$4

  echo -n "#($bin_a) ($bin_b) (w:$x,h:$y):"
  diff \
    <( $bin_a $x $y ) \
    <( $bin_b $x $y )
  if [[ $? != 0 ]] ; then echo "FAIL" ; else echo "pass" ; fi

}

gilbert3dpp_cmp () {
  local bin_a=$1
  local bin_b=$2
  local x=$3
  local y=$4
  local z=$5

  echo -n "#($bin_a) ($bin_b) (w:$x,h:$y,d:$z):"
  diff \
    <( $bin_a $x $y $z ) \
    <( $bin_b $x $y $z )
  if [[ $? != 0 ]] ; then echo "FAIL" ; else echo "pass" ; fi

}



gilbert2dpp_cmp_js_py_c () {
  local x=$1
  local y=$2

  bin_base="node ./gilbert3dpp.js xy"
  bin_cmp="./gilbert3dpp-c d2xy"

  gilbert2dpp_cmp "$bin_base" "$bin_cmp" $x $y

  bin_cmp="./gilbert3dpp.py -a xy"
  gilbert2dpp_cmp "$bin_base" "$bin_cmp" $x $y
}

gilbert3dpp_cmp_js_py_c () {
  local x=$1
  local y=$2
  local z=$3

  bin_base="node ./gilbert3dpp.js xyz"

  bin_cmp="./gilbert3dpp-c d2xyz"
  gilbert3dpp_cmp "$bin_base" "$bin_cmp" $x $y $z

  bin_cmp="./gilbert3dpp.py -a xyz"
  gilbert3dpp_cmp "$bin_base" "$bin_cmp" $x $y $z
}

if [[ $COMPARE_TESTS ]] ; then

  x=10 ; y=2
  gilbert2dpp_cmp_js_py_c $x $y

  x=40 ; y=30 ; z=20
  gilbert3dpp_cmp_js_py_c $x $y $z

  x=41 ; y=30 ; z=20
  gilbert3dpp_cmp_js_py_c $x $y $z

  x=41 ; y=31 ; z=20
  gilbert3dpp_cmp_js_py_c $x $y $z

  x=41 ; y=31 ; z=21
  gilbert3dpp_cmp_js_py_c $x $y $z

fi


if [[ $C_TESTS ]] ; then

  if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
    echo "# C"
  fi

  ## d2xy (sync)
  ##
  if [[ "$VERBOSE_LEVEL" -gt 0 ]] ; then
    echo -n "# C: gilbert2d d2xy (sync): "
    if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
      echo ""
    fi
  fi

  for y in `seq 1 $N` ; do
    for x in `seq 1 $N` ; do

      if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
        echo "## $x $y"
      fi

      a0=`echo "$x%2" | bc`
      if [[ "$a0" = 1 ]] && ( [[ "$y" = 2 ]] || [[ "$z" == 1 ]] ) ; then
        if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
          echo "# skipping (known bad configuration: $x $y $z)"
        fi
        continue
      fi
      if [[ "$a0" = 1 ]] && ( [[ "$y" = 1 ]] || [[ "$z" == 2 ]] ) ; then
        if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
          echo "# skipping (known bad configuration: $x $y $z)"
        fi
        continue
      fi

      ./gilbert3dpp-c d2xy $x $y | \
        grep -v '#' | \
        ./gilbert-check-curve $CHECK_OPT

      if [[ "$?" != 0 ]] ; then
        ERR=1
        break
      fi


    done
    if [[ "$ERR" != 0 ]] ; then break ; fi
  done

  if [[ "$ERR" != 0 ]] ; then
    echo FAIL
    exit -1
  fi
  echo pass

  ## xy2d (sync)
  ##
  if [[ "$VERBOSE_LEVEL" -gt 0 ]] ; then
    echo -n "# C: gilbert2d xy2d (sync): "
    if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
      echo ""
    fi
  fi


  for y in `seq 1 $N` ; do
    for x in `seq 1 $N` ; do

      if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
        echo "## $x $y"
      fi

      a0=`echo "$x%2" | bc`
      if [[ "$a0" = 1 ]] && ( [[ "$y" = 2 ]] || [[ "$z" == 1 ]] ) ; then
        if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
          echo "# skipping (known bad configuration: $x $y $z)"
        fi
        continue
      fi
      if [[ "$a0" = 1 ]] && ( [[ "$y" = 1 ]] || [[ "$z" == 2 ]] ) ; then
        if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
          echo "# skipping (known bad configuration: $x $y $z)"
        fi
        continue
      fi

      ./gilbert3dpp-c xy2d $x $y | \
        grep -v '#' | \
        sort -n | cut -f2- -d' ' | \
        ./gilbert-check-curve $CHECK_OPT

      if [[ "$?" != 0 ]] ; then
        ERR=1
        break
      fi


    done
    if [[ "$ERR" != 0 ]] ; then break ; fi
  done

  if [[ "$ERR" != 0 ]] ; then
    echo FAIL
    exit -1
  fi
  echo pass



  ## d2xy (sync)
  ##
  if [[ "$VERBOSE_LEVEL" -gt 0 ]] ; then
    echo -n "# C: gilbert3d d2xyz (sync): "
    if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
      echo ""
    fi
  fi

  for z in `seq 1 $N` ; do
    for y in `seq 1 $N` ; do
      for x in `seq 1 $N` ; do

        if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
          echo "## $x $y $z"
        fi

        a0=`echo "$x%2" | bc`
        if [[ "$a0" = 1 ]] && ( [[ "$y" = 2 ]] || [[ "$z" == 1 ]] ) ; then
          if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
            echo "# skipping (known bad configuration: $x $y $z)"
          fi
          continue
        fi
        if [[ "$a0" = 1 ]] && ( [[ "$y" = 1 ]] || [[ "$z" == 2 ]] ) ; then
          if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
            echo "# skipping (known bad configuration: $x $y $z)"
          fi
          continue
        fi

        ./gilbert3dpp-c d2xyz $x $y $z | \
          grep -v '#' | \
          ./gilbert-check-curve $CHECK_OPT

        if [[ "$?" != 0 ]] ; then

          ERR=1
          break
        fi


      done
      if [[ "$ERR" != 0 ]] ; then break ; fi
    done
    if [[ "$ERR" != 0 ]] ; then break ; fi
  done

  if [[ "$ERR" != 0 ]] ; then

    echo FAIL

    echo "# FAIL C g3d d2xyz $x $y $z"

    exit -1
  fi
  echo pass

  ## xyz2d (sync)
  ##
  if [[ "$VERBOSE_LEVEL" -gt 0 ]] ; then
    echo -n "# C: gilbert3d xyz2d (sync): "
    if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
      echo ""
    fi
  fi


  for z in `seq 1 $N` ; do
    for y in `seq 1 $N` ; do
      for x in `seq 1 $N` ; do

        if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
          echo "## $x $y $z"
        fi

        a0=`echo "$x%2" | bc`
        if [[ "$a0" = 1 ]] && ( [[ "$y" = 2 ]] || [[ "$z" == 1 ]] ) ; then
          if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
            echo "# skipping (known bad configuration: $x $y $z)"
          fi
          continue
        fi
        if [[ "$a0" = 1 ]] && ( [[ "$y" = 1 ]] || [[ "$z" == 2 ]] ) ; then
          if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
            echo "# skipping (known bad configuration: $x $y $z)"
          fi
          continue
        fi

        ./gilbert3dpp-c xyz2d $x $y $z | \
          grep -v '#' | \
          sort -n | cut -f2- -d' ' | \
          ./gilbert-check-curve $CHECK_OPT

        if [[ "$?" != 0 ]] ; then
          ERR=1
          break
        fi


      done
      if [[ "$ERR" != 0 ]] ; then break ; fi
    done
    if [[ "$ERR" != 0 ]] ; then break ; fi
  done

  if [[ "$ERR" != 0 ]] ; then
    echo FAIL

    echo "# FAIL C g3d xyz2d $x $y $z"

    exit -1
  fi
  echo pass

fi

####
####
####

if [[ $PYTHON_TESTS ]] ; then

  if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
    echo "# Python"
  fi

  ## async
  ##
  if [[ "$VERBOSE_LEVEL" -gt 0 ]] ; then
    echo -n "#Py: gilbert2d (async): "
    if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
      echo ""
    fi
  fi

  for y in `seq 1 $N` ; do
    for x in `seq 1 $N` ; do

      if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
        echo "## $x $y"
      fi

      a0=`echo "$x%2" | bc`
      if [[ "$a0" = 1 ]] && ( [[ "$y" = 2 ]] || [[ "$z" == 1 ]] ) ; then
        if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
          echo "# skipping (known bad configuration: $x $y $z)"
        fi
        continue
      fi
      if [[ "$a0" = 1 ]] && ( [[ "$y" = 1 ]] || [[ "$z" == 2 ]] ) ; then
        if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
          echo "# skipping (known bad configuration: $x $y $z)"
        fi
        continue
      fi

      ./gilbert3dpp.py -a xy $x $y | \
        grep -v '#' | \
        ./gilbert-check-curve $CHECK_OPT

      if [[ "$?" != 0 ]] ; then
        ERR=1
        break
      fi


    done
    if [[ "$ERR" != 0 ]] ; then break ; fi
  done

  if [[ "$ERR" != 0 ]] ; then
    echo FAIL
    exit -1
  fi
  echo pass


  ## d2xy
  ##
  if [[ "$VERBOSE_LEVEL" -gt 0 ]] ; then
    echo -n "#Py: gilbert2d d2xy (sync): "
    if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
      echo ""
    fi
  fi

  for y in `seq 1 $N` ; do
    for x in `seq 1 $N` ; do

      if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
        echo "## $x $y"
      fi

      a0=`echo "$x%2" | bc`
      if [[ "$a0" = 1 ]] && ( [[ "$y" = 2 ]] || [[ "$z" == 1 ]] ) ; then
        if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
          echo "# skipping (known bad configuration: $x $y $z)"
        fi
        continue
      fi
      if [[ "$a0" = 1 ]] && ( [[ "$y" = 1 ]] || [[ "$z" == 2 ]] ) ; then
        if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
          echo "# skipping (known bad configuration: $x $y $z)"
        fi
        continue
      fi

      ./gilbert3dpp.py -a d2xy $x $y | \
        grep -v '#' | \
        ./gilbert-check-curve $CHECK_OPT

      if [[ "$?" != 0 ]] ; then
        ERR=1
        break
      fi


    done
    if [[ "$ERR" != 0 ]] ; then break ; fi
  done

  if [[ "$ERR" != 0 ]] ; then
    echo FAIL
    exit -1
  fi
  echo pass

  ## xy2d
  ##
  if [[ "$VERBOSE_LEVEL" -gt 0 ]] ; then
    echo -n "#Py: gilbert2d xy2d (sync): "
    if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
      echo ""
    fi
  fi

  for y in `seq 1 $N` ; do
    for x in `seq 1 $N` ; do

      if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
        echo "## $x $y $z"
      fi


      a0=`echo "$x%2" | bc`
      if [[ "$a0" = 1 ]] && ( [[ "$y" = 2 ]] || [[ "$z" == 1 ]] ) ; then
        if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
          echo "# skipping (known bad configuration: $x $y $z)"
        fi
        continue
      fi
      if [[ "$a0" = 1 ]] && ( [[ "$y" = 1 ]] || [[ "$z" == 2 ]] ) ; then
        if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
          echo "# skipping (known bad configuration: $x $y $z)"
        fi
        continue
      fi

      ./gilbert3dpp.py -a xy2d $x $y | \
        grep -v '#' | \
        sort -n | cut -f2- -d' ' | \
        ./gilbert-check-curve $CHECK_OPT

      if [[ "$?" != 0 ]] ; then
        ERR=1
        break
      fi

    done
    if [[ "$ERR" != 0 ]] ; then break ; fi
  done

  if [[ "$ERR" != 0 ]] ; then
    echo FAIL
    exit -1
  fi
  echo pass


  ## async
  ##
  if [[ "$VERBOSE_LEVEL" -gt 0 ]] ; then
    echo -n "#Py: gilbert3d (async): "
    if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
      echo ""
    fi
  fi

  for z in `seq 1 $N` ; do
    for y in `seq 1 $N` ; do
      for x in `seq 1 $N` ; do

        if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
          echo "## $x $y $z"
        fi


        a0=`echo "$x%2" | bc`
        if [[ "$a0" = 1 ]] && ( [[ "$y" = 2 ]] || [[ "$z" == 1 ]] ) ; then
          if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
            echo "# skipping (known bad configuration: $x $y $z)"
          fi
          continue
        fi
        if [[ "$a0" = 1 ]] && ( [[ "$y" = 1 ]] || [[ "$z" == 2 ]] ) ; then
          if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
            echo "# skipping (known bad configuration: $x $y $z)"
          fi
          continue
        fi

        ./gilbert3dpp.py -a xyz $x $y $z | \
          grep -v '#' | \
          ./gilbert-check-curve $CHECK_OPT

        if [[ "$?" != 0 ]] ; then
          ERR=1
          break
        fi


      done
      if [[ "$ERR" != 0 ]] ; then break ; fi
    done
    if [[ "$ERR" != 0 ]] ; then break ; fi
  done

  if [[ "$ERR" != 0 ]] ; then
    echo FAIL
    exit -1
  fi
  echo pass

  ## d2xyz
  ##
  ## xyz2d
  ##
  if [[ "$VERBOSE_LEVEL" -gt 0 ]] ; then
    echo -n "#Py: gilbert3d d2xyz (sync): "
    if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
      echo ""
    fi
  fi

  for z in `seq 1 $N` ; do
    for y in `seq 1 $N` ; do
      for x in `seq 1 $N` ; do

        if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
          echo "## $x $y $z"
        fi


        a0=`echo "$x%2" | bc`
        if [[ "$a0" = 1 ]] && ( [[ "$y" = 2 ]] || [[ "$z" == 1 ]] ) ; then
          if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
            echo "# skipping (known bad configuration: $x $y $z)"
          fi
          continue
        fi
        if [[ "$a0" = 1 ]] && ( [[ "$y" = 1 ]] || [[ "$z" == 2 ]] ) ; then
          if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
            echo "# skipping (known bad configuration: $x $y $z)"
          fi
          continue
        fi

        ./gilbert3dpp.py -a d2xyz $x $y $z | \
          grep -v '#' | \
          ./gilbert-check-curve $CHECK_OPT

        if [[ "$?" != 0 ]] ; then
          ERR=1
          break
        fi


      done
      if [[ "$ERR" != 0 ]] ; then break ; fi
    done
    if [[ "$ERR" != 0 ]] ; then break ; fi
  done

  if [[ "$ERR" != 0 ]] ; then
    echo FAIL
    exit -1
  fi
  echo pass


  ## xyz2d
  ##
  if [[ "$VERBOSE_LEVEL" -gt 0 ]] ; then
    echo -n "#Py: gilbert3d xyz2d (sync): "
    if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
      echo ""
    fi
  fi

  for z in `seq 1 $N` ; do
    for y in `seq 1 $N` ; do
      for x in `seq 1 $N` ; do

        if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
          echo "## $x $y $z"
        fi


        a0=`echo "$x%2" | bc`
        if [[ "$a0" = 1 ]] && ( [[ "$y" = 2 ]] || [[ "$z" == 1 ]] ) ; then
          if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
            echo "# skipping (known bad configuration: $x $y $z)"
          fi
          continue
        fi
        if [[ "$a0" = 1 ]] && ( [[ "$y" = 1 ]] || [[ "$z" == 2 ]] ) ; then
          if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
            echo "# skipping (known bad configuration: $x $y $z)"
          fi
          continue
        fi

        ./gilbert3dpp.py -a xyz2d $x $y $z | \
          grep -v '#' | \
          sort -n | cut -f2- -d' ' | \
          ./gilbert-check-curve $CHECK_OPT

        if [[ "$?" != 0 ]] ; then
          ERR=1
          break
        fi

      done
      if [[ "$ERR" != 0 ]] ; then break ; fi
    done
    if [[ "$ERR" != 0 ]] ; then break ; fi
  done

  if [[ "$ERR" != 0 ]] ; then
    echo FAIL
    exit -1
  fi
  echo pass

fi

####
####
####

if [[ $JAVASCRIPT_TESTS ]] ; then

  if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
    echo "# JavaScript"
  fi

  if [[ "$VERBOSE_LEVEL" -gt 0 ]] ; then
    echo -n "#JS: gilbert2d xy (async): "
    if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
      echo ""
    fi
  fi


  for y in `seq 1 $N` ; do
    for x in `seq 1 $N` ; do

      if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
        echo "## $x $y"
      fi

      a0=`echo "$x%2" | bc`
      if [[ "$a0" = 1 ]] && [[ "$y" = 2 ]] ; then
        if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
          echo "# skipping (known bad configuration: $x $y)"
        fi
        continue
      fi

      node ./gilbert3dpp.js xy $x $y  | \
        grep -v '#' | \
        ./gilbert-check-curve
      if [[ "$?" != 0 ]] ; then
        ERR=1
        break
      fi
    done
    if [[ "$ERR" != 0 ]] ; then break ; fi
  done

  if [[ "$ERR" != 0 ]] ; then
    echo FAIL
    exit -1
  fi
  echo pass



  if [[ "$VERBOSE_LEVEL" -gt 0 ]] ; then
    echo -n "#JS: gilbert2d d2xy (sync): "
    if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
      echo ""
    fi
  fi

  for y in `seq 1 $N` ; do
    for x in `seq 1 $N` ; do
      if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
        echo "## $x $y"
      fi


      a0=`echo "$x%2" | bc`
      if [[ "$a0" = 1 ]] && [[ "$y" = 2 ]] ; then
        if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
          echo "# skipping (known bad configuration: $x $y)"
        fi
        continue
      fi

      node ./gilbert3dpp.js d2xy $x $y | \
        grep -v '#' | \
        ./gilbert-check-curve
      if [[ "$?" != 0 ]] ; then
        ERR=1
        break
      fi
    done
    if [[ "$ERR" != 0 ]] ; then break ; fi
  done

  if [[ "$ERR" != 0 ]] ; then
    echo FAIL
    exit -1
  fi

  echo pass


  if [[ "$VERBOSE_LEVEL" -gt 0 ]] ; then
    echo -n "#JS: gilbert2d xy2d (sync): "
    if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
      echo ""
    fi
  fi


  for y in `seq 1 $N` ; do
    for x in `seq 1 $N` ; do
      if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
        echo "## $x $y"
      fi


      a0=`echo "$x%2" | bc`
      if [[ "$a0" = 1 ]] && [[ "$y" = 2 ]] ; then
        if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
          echo "# skipping (known bad configuration: $x $y)"
        fi
        continue
      fi

      node ./gilbert3dpp.js xy2d $x $y | \
        grep -v '#' | \
        sort -n | \
        cut -f2- -d' ' | ./gilbert-check-curve

      if [[ "$?" != 0 ]] ; then
        ERR=1
        break
      fi

    done
    if [[ "$ERR" != 0 ]] ; then break ; fi
  done

  if [[ "$ERR" != 0 ]] ; then
    echo FAIL
    exit -1
  fi
  echo pass


  if [[ "$VERBOSE_LEVEL" -gt 0 ]] ; then
    echo -n "#JS: gilbert3d xyz (async): "
    if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
      echo ""
    fi
  fi


  ## JS version
  ##
  for z in `seq 1 $N` ; do
    for y in `seq 1 $N` ; do
      for x in `seq 1 $N` ; do

        if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
          echo "## $x $y $z"
        fi


        a0=`echo "$x%2" | bc`
        if [[ "$a0" = 1 ]] && ( [[ "$y" = 2 ]] || [[ "$z" == 1 ]] ) ; then
          if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
            echo "# skipping (known bad configuration: $x $y $z)"
          fi
          continue
        fi
        if [[ "$a0" = 1 ]] && ( [[ "$y" = 1 ]] || [[ "$z" == 2 ]] ) ; then
          if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
            echo "# skipping (known bad configuration: $x $y $z)"
          fi
          continue
        fi

        node ./gilbert3dpp.js xyz $x $y $z | \
          grep -v '#' | \
          ./gilbert-check-curve

        if [[ "$?" != 0 ]] ; then
          ERR=1
          break
        fi


      done
      if [[ "$ERR" != 0 ]] ; then break ; fi
    done
    if [[ "$ERR" != 0 ]] ; then break ; fi
  done

  if [[ "$ERR" != 0 ]] ; then
    echo FAIL
    exit -1
  fi
  echo pass


  if [[ "$VERBOSE_LEVEL" -gt 0 ]] ; then
    echo -n "#JS: gilbert3d d2xyz (sync): "
    if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
      echo ""
    fi
  fi


  for z in `seq 1 $N` ; do
    for y in `seq 1 $N` ; do
      for x in `seq 1 $N` ; do
        if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
          echo "## $x $y $z"
        fi


        a0=`echo "$x%2" | bc`
        if [[ "$a0" = 1 ]] && ( [[ "$y" = 2 ]] || [[ "$z" == 1 ]] ) ; then
          if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
            echo "# skipping (known bad configuration: $x $y $z)"
          fi
          continue
        fi
        if [[ "$a0" = 1 ]] && ( [[ "$y" = 1 ]] || [[ "$z" == 2 ]] ) ; then
          if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
            echo "# skipping (known bad configuration: $x $y $z)"
          fi
          continue
        fi


        node ./gilbert3dpp.js d2xyz $x $y $z | \
          grep -v '#' | \
          ./gilbert-check-curve

        if [[ "$?" != 0 ]] ; then
          ERR=1
          break
        fi


      done
      if [[ "$ERR" != 0 ]] ; then break ; fi
    done
    if [[ "$ERR" != 0 ]] ; then break ; fi
  done

  if [[ "$ERR" != 0 ]] ; then
    echo FAIL
    exit -1
  fi
  echo pass


  if [[ "$VERBOSE_LEVEL" -gt 0 ]] ; then
    echo -n "#JS: gilbert3d xyz2d (sync): "
    if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
      echo ""
    fi
  fi


  for z in `seq 1 $N` ; do
    for y in `seq 1 $N` ; do
      for x in `seq 1 $N` ; do
        if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
          echo "## $x $y $z"
        fi


        a0=`echo "$x%2" | bc`
        if [[ "$a0" = 1 ]] && ( [[ "$y" = 2 ]] || [[ "$z" == 1 ]] ) ; then
          if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
            echo "# skipping (known bad configuration: $x $y $z)"
          fi
          continue
        fi
        if [[ "$a0" = 1 ]] && ( [[ "$y" = 1 ]] || [[ "$z" == 2 ]] ) ; then
          if [[ "$VERBOSE_LEVEL" -gt 1 ]] ; then
            echo "# skipping (known bad configuration: $x $y $z)"
          fi
          continue
        fi

        node ./gilbert3dpp.js xyz2d $x $y $z | \
          grep -v '#' | \
          sort -n | \
          cut -f2- -d' ' | ./gilbert-check-curve

        if [[ "$?" != 0 ]] ; then
          ERR=1
          break
        fi


      done
      if [[ "$ERR" != 0 ]] ; then break ; fi
    done
    if [[ "$ERR" != 0 ]] ; then break ; fi
  done

  if [[ "$ERR" != 0 ]] ; then
    echo FAIL
    exit -1
  fi
  echo pass

fi

exit 0

