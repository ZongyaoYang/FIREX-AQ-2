#!/bin/bash -l
date
echo "Host is " `hostname`
pushd ~/FIREX-AQ

source ../loadPython
python3 firex-plot.py year=today species=isopr,O3_SFC,PM2_5_DRY_SFC,co_fire,nh3,co_anth,c2h6,PM2_5_DRY source=aq-watch domain=d01 height=-1.0		# real
#python3 firex-plot.py year=2025 month=5 day=13 species=isopr,O3_SFC,PM2_5_DRY_SFC,co_fire,nh3,co_anth,c2h6,PM2_5_DRY source=aq-watch domain=d01 height=-1.0
#python3 firex-plot.py year=2025 month=5 day=14 species=isopr,O3_SFC,PM2_5_DRY_SFC,co_fire,nh3,co_anth,c2h6,PM2_5_DRY source=aq-watch domain=d01 height=-1.0
#python3 firex-plot.py year=2025 month=5 day=15 species=isopr,O3_SFC,PM2_5_DRY_SFC,co_fire,nh3,co_anth,c2h6,PM2_5_DRY source=aq-watch domain=d01 height=-1.0
#python3 firex-plot.py year=2025 month=5 day=7 hour=20 species=isopr source=aq-watch domain=d01,d02 height=-1.0
#python3 firex-plot.py year=2024 month=1 day=2 species=O3_SFC source=aq-watch domain=d01,d02 height=-1.0
#python3 firex-plot.py year=2021 month=2 species=O3_SFC source=aq-watch domain=d01,d02 height=-1.0
#python3 firex-plot.py year=2021 month=3 species=O3_SFC source=aq-watch domain=d01,d02 height=-1.0
#python3 firex-plot.py year=2021 month=6 day=17 hour=11 species=nh3,c2h6 source=aq-watch domain=d02
#python3 firex-plot.py year=2021 month=6 day=17 hour=11 species=co_anth,co_fire source=aq-watch domain=d02
#python3 firex-plot.py year=2021 month=6 day=17 hour=11 species=PM2_5_DRY source=aq-watch domain=d02
#python3 firex-plot.py year=today species=nh3,c2h6,co_anth,co_fire,PM2_5_DRY source=aq-watch domain=d02
#python3 firex-plot.py year=2023 month=12 species=vent_rate,PBLH,RH source=aq-watch domain=d02,d03	# backfill

# ./copyToModeling1.sh

popd
date

