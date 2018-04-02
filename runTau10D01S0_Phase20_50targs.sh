# bash script to average over targets
# assumes you are in the top level of BEAM diectory
# before running you should compile 
# set up code (make_particles_wakes.f90 or mpw_iquad_xy.f90) 
# and Monte Carlo code (mc.F90)
# see the README.txt for more information
alpha=20
shaparam=100
prefile=mc-beamS0Phase
# LitGeom suffix means both sun and observer are on the same side of ring plane
postfile=LitGeom.in
optdepth=tau10
density=D01
binParm=nbins40
#set home directory to current diectory
home1=$(pwd)
echo ${home1}
#set number of targets
MAXTARGETS=50
results=${home1}/BEAM_${optdepth}${density}S0Phase${alpha}SHA${shaparam}Rev3_50targs # specific target params
echo ${results}
mkdir -p ${results}
cd ${results}
cp ${home1}/input/if_sum.out .
cp ${home1}/input/flux_sum.out .
#loop over targets
for iter in `seq 1 $MAXTARGETS`; do 
echo $iter
setup=$home1/bin/wakex
part_input=${home1}/input/part-no_wakes${density}${optdepth}_${binParm}.in
echo ${part_input}
${setup} < ${part_input} > part_output
target_dir=${results}/targetnum_${iter}/ 
echo ${target_dir}
mkdir -p ${target_dir}
cp nowakes_${density}_${optdepth}_${binParm}.partdat ${target_dir}
cp part_output ${target_dir}
# these all have sunlat=26.7, obslat=41.7 
beam_input=${home1}/input/${prefile}${alpha}sha${shaparam}${optdepth}${density}${postfile}
echo ${beam_input}
${home1}/bin/mcx < ${beam_input} > beam_output
cp if.out ${target_dir}
cp flux.out ${target_dir}
cp beam_output ${target_dir}
cp nphotpass.out ${target_dir}
python ${home1}/python/SumFlux.py flux.out flux_sum.out if.out if_sum.out $iter $MAXTARGETS
mv flux_new.out flux_sum.out
mv if_new.out if_sum.out
done # end of loop  

