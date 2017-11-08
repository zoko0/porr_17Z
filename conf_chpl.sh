#!/bin/bash
export CHPL_HOME=~/Pulpit/magisterka_sem3/porr/projekt/chapel-master
export CHPL_HOST_PLATFORM=`$CHPL_HOME/util/chplenv/chpl_platform.py`
export PATH="$PATH":"$CHPL_HOME/bin/$CHPL_HOST_PLATFORM"
export MANPATH="$MANPATH":"$CHPL_HOME"/man
