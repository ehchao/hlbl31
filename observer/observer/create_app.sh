#!/bin/bash



if [ "$#" -ne 1 ] ; then
  echo "Usage: $0 APPNAME" >&2
  exit 1
fi

APPNAME="$1"

if [ -d "apps/$1" ] || [ -d "devel/$1" ] ; then
  echo "App/Devel directory with that name already exists"
  exit 1
fi


mkdir apps/"$APPNAME"

touch apps/"$APPNAME"/calculate_"$APPNAME".cc  
touch apps/"$APPNAME"/read_"$APPNAME"_parms.cc  
touch apps/"$APPNAME"/write_"$APPNAME"_data.cc

mkdir devel/"$APPNAME"

cat << EOF > "devel/${APPNAME}/${APPNAME}.cc" 
#define MAIN_PROGRAM
#include "observer.h"

using namespace QDP;

int main(int argc, char *argv[])
{

  int cfno;
  int status[3];
  StopWatch swatch;


  /* Initialize Observer Code */
  observer_init(argc, argv);
  multi1d < LatticePropagatorD > quark_props;
  print_info();


  for (cfno = observer.first_cnfg; cfno <= observer.last_cnfg;
       cfno += observer.step_cnfg)
  {
    start_time(compute_config);
    import_config(cfno);
    chs_ubnd(-1);
    start_time(Deflation_Setup);
    dfl_modes(status);
    stop_time(Deflation_Setup);
    for (int j = 0; j < observer.no_src_pos; j++)
    {
      calculate_propagators(observer.src_pos[j], quark_props);
      /* App goes here */

    }
    recompute_propagators();
    stop_time(compute_config);
  }

  QDP_finalize();
}

EOF



cp devel/template/Makefile devel/${APPNAME}/
command=$(echo "sed -i 's/TEMPLATENAME/${APPNAME}/g' devel/${APPNAME}/Makefile")
echo $command
eval $command

#mkdir devel/"$APPNAME"

