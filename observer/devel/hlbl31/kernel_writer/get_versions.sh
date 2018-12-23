#!/bin/bash


TOPLEVEL=$(git rev-parse --show-toplevel)


OBSERV_REL=$(git describe --abbrev=4 --dirty --always --tags)
cd $TOPLEVEL
cd qdpxx/
QDP_REL=$(git describe --abbrev=4 --dirty --always --tags)
cd $TOPLEVEL
cd qmp/
QMP_REL=$(git describe --abbrev=4 --dirty --always --tags)

echo "-DOBSERV_VERS=\"$OBSERV_REL\" -DQDP_VERS=\"$QDP_REL\" -DQMP_VERS=\"$QMP_REL\""
