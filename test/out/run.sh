ccol=~/src/prog/cCOLUMBUS
echo molint
${ccol}/molint < int.in > int.out
echo cdip
${ccol}/cdip < int.in > cdip.out
echo pkfs
${ccol}/pkfl < dum.in > pkfl.out
echo cscf
${ccol}/cscf < scf.in > scf.out
echo CISD_ground
${ccol}/molci < ci1.in > ci1.out
mv CSF CSFG
mv CIVEC CIVECG
echo CISD_EE
${ccol}/molci < ci2.in > ci2.out
