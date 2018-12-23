
#ifndef OBSERVER_VERSION_H
#define OBSERVER_VERSION_H

#define STR2(s) #s
#define STR(s) STR2(s)
#ifndef OBSERV_VERS
#error "Version string needs to be provided for OBSERVER!\n Use the script get_versions.sh in the top level!"
#endif
#define observer_RELEASE STR(OBSERV_VERS)
#ifndef QDP_VERS
#error "Version string needs to be provided for QDP! \
Use the script get_versions.sh in the top level! \
VERSIONINFO=$(shell ../../../get_versions.sh) \
 \
CFLAGS =  ... $(VERSIONINFO)"


#endif
#define qdp_RELEASE STR(QDP_VERS)


#ifndef QMP_VERS
#error "Version string needs to be provided for QMP!\n Use the script get_versions.sh in the top level!"
#endif
#define qmp_RELEASE STR(QMP_VERS)

#endif
